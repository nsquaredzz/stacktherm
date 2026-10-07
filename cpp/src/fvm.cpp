// Cell-centred finite-volume discretisation of steady heat conduction.
//
// Each face conductance is the series sum of the two half-cell resistances and
// the interface (Kapitza / barrier) resistance on that face, which makes
// material jumps and boundary resistance exact for face-aligned interfaces:
//
//     g = A_face / (d_1 / (2 k_1) + R_int + d_2 / (2 k_2))
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace st {

void Grid::allocate() {
    size_t m = n();
    for (Vec* v : {&kx, &ky, &kz, &rx, &ry, &rz, &q, &rhoc}) v->assign(m, 0.0);
    mat.assign(m, 0);
    if (layer_of_k.empty()) layer_of_k.assign(nz(), 0);
}

std::vector<std::pair<int, int>> Grid::layer_slices() const {
    std::vector<std::pair<int, int>> out;
    int start = 0;
    for (int k = 1; k <= nz(); ++k) {
        if (k == nz() || layer_of_k[k] != layer_of_k[k - 1]) {
            out.emplace_back(start, k);
            start = k;
        }
    }
    return out;
}

// Edges on [0, t] whose cells grow geometrically away from both ends, so thick
// uniform layers stay cheap while the zone next to a structured neighbour is
// resolved.
Vec graded_edges(double t, double h_bot, double h_top, double ratio, double h_max) {
    if (h_max < 0) h_max = t;
    h_bot = std::min(h_bot, t);
    h_top = std::min(h_top, t);
    Vec lo, hi;
    double pos_lo = 0, pos_hi = 0, a = h_bot, b = h_top;
    while (pos_lo + pos_hi + std::min(a, b) < t) {
        if (a <= b) { lo.push_back(a); pos_lo += a; a = std::min(a * ratio, h_max); }
        else { hi.push_back(b); pos_hi += b; b = std::min(b * ratio, h_max); }
    }
    Vec sizes(lo);
    sizes.insert(sizes.end(), hi.rbegin(), hi.rend());
    if (sizes.empty()) sizes.push_back(t);
    else {
        double sum = std::accumulate(sizes.begin(), sizes.end(), 0.0), gap = t - sum;
        double largest = *std::max_element(sizes.begin(), sizes.end());
        if (gap > 0.5 * largest || sizes.size() < 2) sizes.insert(sizes.begin() + lo.size(), gap);
    }
    double sum = std::accumulate(sizes.begin(), sizes.end(), 0.0);
    Vec edges(sizes.size() + 1, 0.0);
    for (size_t i = 0; i < sizes.size(); ++i) edges[i + 1] = edges[i] + sizes[i] * (t / sum);
    edges.back() = t;
    return edges;
}

void Csr::mul(const Vec& x, Vec& y) const {
    y.assign(n, 0.0);
    for (int r = 0; r < n; ++r) {
        double s = 0.0;
        for (ptrdiff_t p = ptr[r]; p < ptr[r + 1]; ++p) s += val[p] * x[col[p]];
        y[r] = s;
    }
}

System assemble(const Grid& g, const BC& bottom, const BC& top, bool periodic_x, bool periodic_y,
                double jump_x, double jump_y, bool sources) {
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    const size_t n = g.n(), nxy = size_t(nx) * ny;
    System s;
    s.grid = &g;
    s.bottom = bottom;
    s.top = top;
    s.jump[0] = jump_x;
    s.jump[1] = jump_y;
    s.gx.assign(n, 0.0);
    s.gy.assign(n, 0.0);
    s.gz.assign(nz > 1 ? n - nxy : 0, 0.0);
    s.gbot.assign(nxy, 0.0);
    s.gtop.assign(nxy, 0.0);

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = g.idx(i, j, k);
                int ip = (i + 1) % nx, jp = (j + 1) % ny;
                if (nx > 1 && (periodic_x || i < nx - 1)) {
                    size_t e = g.idx(ip, j, k);
                    s.gx[c] = g.dy(j) * g.dz(k) /
                              (g.dx(i) / (2 * g.kx[c]) + g.dx(ip) / (2 * g.kx[e]) + g.rx[c]);
                }
                if (ny > 1 && (periodic_y || j < ny - 1)) {
                    size_t e = g.idx(i, jp, k);
                    s.gy[c] = g.dx(i) * g.dz(k) /
                              (g.dy(j) / (2 * g.ky[c]) + g.dy(jp) / (2 * g.ky[e]) + g.ry[c]);
                }
                if (k < nz - 1) {
                    size_t e = c + nxy;
                    s.gz[c] = g.dx(i) * g.dy(j) /
                              (g.dz(k) / (2 * g.kz[c]) + g.dz(k + 1) / (2 * g.kz[e]) + g.rz[c]);
                }
            }

    Vec diag(n, 0.0);
    s.b.assign(n, 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = g.idx(i, j, k);
                size_t w = g.idx((i + nx - 1) % nx, j, k), so = g.idx(i, (j + ny - 1) % ny, k);
                double d = s.gx[c] + s.gx[w] + s.gy[c] + s.gy[so];
                if (k < nz - 1) d += s.gz[c];
                if (k > 0) d += s.gz[c - nxy];
                diag[c] = d;
                if (sources) s.b[c] = g.q[c] * g.dx(i) * g.dy(j) * g.dz(k);
            }
    if (jump_x != 0.0)
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
                size_t last = g.idx(nx - 1, j, k), first = g.idx(0, j, k);
                s.b[last] += s.gx[last] * jump_x;
                s.b[first] -= s.gx[last] * jump_x;
            }
    if (jump_y != 0.0)
        for (int k = 0; k < nz; ++k)
            for (int i = 0; i < nx; ++i) {
                size_t last = g.idx(i, ny - 1, k), first = g.idx(i, 0, k);
                s.b[last] += s.gy[last] * jump_y;
                s.b[first] -= s.gy[last] * jump_y;
            }

    auto boundary = [&](const BC& bc, int k, Vec& gb) {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = g.idx(i, j, k);
                double area = g.dx(i) * g.dy(j);
                if (bc.fixes_temperature()) {
                    double gg = area / (g.dz(k) / (2 * g.kz[c]) + bc.resistance());
                    gb[size_t(j) * nx + i] = gg;
                    diag[c] += gg;
                    s.b[c] += gg * bc.T;
                } else if (bc.kind == BC::Flux) {
                    s.b[c] += bc.q * area;
                }
            }
    };
    boundary(bottom, 0, s.gbot);
    boundary(top, nz - 1, s.gtop);

    s.pinned = !(bottom.fixes_temperature() || top.fixes_temperature());
    if (s.pinned) {
        // Pure Neumann / periodic problem: the solution is defined up to a
        // constant.  A diagonal shift on one cell selects T[0] = 0 without
        // perturbing the fluxes, because the right-hand side sums to zero.
        diag[0] += std::accumulate(diag.begin(), diag.end(), 0.0) / double(n);
    }

    Csr& A = s.A;
    A.n = int(n);
    A.ptr.assign(n + 1, 0);
    A.col.reserve(7 * n);
    A.val.reserve(7 * n);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = g.idx(i, j, k);
                std::pair<ptrdiff_t, double> e[7];
                int m = 0;
                auto add = [&](size_t column, double value) {
                    if (value == 0.0) return;
                    for (int t = 0; t < m; ++t)
                        if (e[t].first == ptrdiff_t(column)) { e[t].second += value; return; }
                    e[m++] = {ptrdiff_t(column), value};
                };
                add(c, diag[c]);
                size_t w = g.idx((i + nx - 1) % nx, j, k), ea = g.idx((i + 1) % nx, j, k);
                size_t so = g.idx(i, (j + ny - 1) % ny, k), no = g.idx(i, (j + 1) % ny, k);
                if (ea != c) add(ea, -s.gx[c]);
                if (w != c) add(w, -s.gx[w]);
                if (no != c) add(no, -s.gy[c]);
                if (so != c) add(so, -s.gy[so]);
                if (k < nz - 1) add(c + nxy, -s.gz[c]);
                if (k > 0) add(c - nxy, -s.gz[c - nxy]);
                std::sort(e, e + m);
                for (int t = 0; t < m; ++t) { A.col.push_back(e[t].first); A.val.push_back(e[t].second); }
                A.ptr[c + 1] = ptrdiff_t(A.col.size());
            }
    return s;
}

void System::face_flows(const Vec& T, Vec& fx, Vec& fy, Vec& fz) const {
    const Grid& g = *grid;
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    const size_t nxy = size_t(nx) * ny;
    fx.assign(g.n(), 0.0);
    fy.assign(g.n(), 0.0);
    fz.assign(gz.size(), 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = g.idx(i, j, k);
                double dtx = T[c] - T[g.idx((i + 1) % nx, j, k)] - (i == nx - 1 ? jump[0] : 0.0);
                double dty = T[c] - T[g.idx(i, (j + 1) % ny, k)] - (j == ny - 1 ? jump[1] : 0.0);
                fx[c] = gx[c] * dtx;
                fy[c] = gy[c] * dty;
                if (k < nz - 1) fz[c] = gz[c] * (T[c] - T[c + nxy]);
            }
}

void System::boundary_flows(const Vec& T, Vec& qbot, Vec& qtop) const {
    const Grid& g = *grid;
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    qbot.assign(size_t(nx) * ny, 0.0);
    qtop.assign(size_t(nx) * ny, 0.0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            size_t col = size_t(j) * nx + i;
            double area = g.dx(i) * g.dy(j);
            if (bottom.fixes_temperature()) qbot[col] = gbot[col] * (T[g.idx(i, j, 0)] - bottom.T);
            else if (bottom.kind == BC::Flux) qbot[col] = -bottom.q * area;
            if (top.fixes_temperature()) qtop[col] = gtop[col] * (T[g.idx(i, j, nz - 1)] - top.T);
            else if (top.kind == BC::Flux) qtop[col] = -top.q * area;
        }
}

void System::cell_flux(const Vec& T, Vec& qx, Vec& qy, Vec& qz) const {
    const Grid& g = *grid;
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    const size_t nxy = size_t(nx) * ny;
    Vec fx, fy, fz, qb, qt;
    face_flows(T, fx, fy, fz);
    boundary_flows(T, qb, qt);
    qx.assign(g.n(), 0.0);
    qy.assign(g.n(), 0.0);
    qz.assign(g.n(), 0.0);
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = g.idx(i, j, k), col = size_t(j) * nx + i;
                qx[c] = 0.5 * (fx[c] + fx[g.idx((i + nx - 1) % nx, j, k)]) / (g.dy(j) * g.dz(k));
                qy[c] = 0.5 * (fy[c] + fy[g.idx(i, (j + ny - 1) % ny, k)]) / (g.dx(i) * g.dz(k));
                double below = k == 0 ? -qb[col] : fz[c - nxy];
                double above = k == nz - 1 ? qt[col] : fz[c];
                qz[c] = 0.5 * (below + above) / (g.dx(i) * g.dy(j));
            }
}

}  // namespace st
