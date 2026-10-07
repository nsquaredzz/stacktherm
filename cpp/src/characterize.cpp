// Cell problems that turn a voxel grid into per-layer effective properties.
//
// through-plane: unit temperature drop bottom to top, no sources.  Plane-averaged
//   face temperatures at every layer boundary split the total resistance into
//   one term per layer and one per inter-layer interface.
// in-plane x and y: unit mean gradient with a periodic fluctuation, adiabatic
//   top and bottom.  The volume-averaged flux in each layer gives its in-situ
//   lateral conductivity.
// operating: real boundary conditions and heat sources, for the temperature field.
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace st {

namespace {

// Area-averaged temperature just below (low) and just above (high) every z face.
// Index f is the face under z cell f.
void plane_temperatures(const System& s, const Vec& T, Vec& t_low, Vec& t_high) {
    const Grid& g = *s.grid;
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    const size_t nxy = size_t(nx) * ny;
    Vec fx, fy, fz, qb, qt;
    s.face_flows(T, fx, fy, fz);
    s.boundary_flows(T, qb, qt);
    const double total = g.lx() * g.ly();
    t_low.assign(nz + 1, 0.0);
    t_high.assign(nz + 1, 0.0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double area = g.dx(i) * g.dy(j), w = area / total;
            const size_t col = size_t(j) * nx + i;
            for (int k = 0; k < nz - 1; ++k) {
                size_t c = g.idx(i, j, k);
                double flux = fz[c] / area;
                t_low[k + 1] += w * (T[c] - flux * g.dz(k) / (2 * g.kz[c]));
                t_high[k + 1] += w * (T[c + nxy] + flux * g.dz(k + 1) / (2 * g.kz[c + nxy]));
            }
            size_t c0 = g.idx(i, j, 0), c1 = g.idx(i, j, nz - 1);
            t_high[0] += w * (T[c0] - (qb[col] / area) * g.dz(0) / (2 * g.kz[c0]));
            t_low[nz] += w * (T[c1] - (qt[col] / area) * g.dz(nz - 1) / (2 * g.kz[c1]));
        }
    t_low[0] = t_high[0];
    t_high[nz] = t_low[nz];
}

}  // namespace

ThroughPlane through_plane(const Grid& g, double tol) {
    System s = assemble(g, BC::dirichlet(1.0), BC::dirichlet(0.0), true, true, 0.0, 0.0, false);
    ThroughPlane out;
    out.T = LinearSolver(s.A, tol).solve(s.b, &out.info);
    Vec qb, qt, t_low, t_high;
    s.boundary_flows(out.T, qb, qt);
    double flow = 0.0;
    for (double v : qt) flow += v;
    const double flux = flow / (g.lx() * g.ly());
    plane_temperatures(s, out.T, t_low, t_high);
    for (auto [k0, k1] : g.layer_slices()) {
        out.r_layer.push_back((t_high[k0] - t_low[k1]) / flux);
        out.r_int.push_back(k1 < g.nz() ? (t_low[k1] - t_high[k1]) / flux : 0.0);
    }
    out.r_total = 1.0 / flux;
    return out;
}

Characterization characterize(const Grid& g, double tol) {
    const int nx = g.nx(), ny = g.ny();
    const auto slices = g.layer_slices();
    const size_t nl = slices.size();
    Characterization ch;
    ThroughPlane tp = through_plane(g, tol);
    ch.solves["through_plane"] = tp.info;

    Vec volume(nl, 0.0), kxl(nl, 0.0), kyl(nl, 0.0), kzmix(nl, 0.0), rhoc(nl, 0.0);
    for (size_t l = 0; l < nl; ++l)
        for (int k = slices[l].first; k < slices[l].second; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    size_t c = g.idx(i, j, k);
                    double v = g.dx(i) * g.dy(j) * g.dz(k);
                    volume[l] += v;
                    kzmix[l] += g.kz[c] * v;
                    rhoc[l] += g.rhoc[c] * v;
                    if (nx == 1) kxl[l] += g.kx[c] * v;      // no structure: arithmetic mean is exact
                    if (ny == 1) kyl[l] += g.ky[c] * v;
                }

    std::unique_ptr<LinearSolver> solver;
    for (int axis = 0; axis < 2; ++axis) {
        if ((axis == 0 ? nx : ny) == 1) continue;
        System s = assemble(g, BC(), BC(), true, true, axis == 0 ? g.lx() : 0.0, axis == 1 ? g.ly() : 0.0,
                            false);                           // mean gradient of 1 K/m
        if (!solver) solver.reset(new LinearSolver(s.A, tol));   // the matrix is the same for both axes
        SolveInfo info;
        Vec T = solver->solve(s.b, &info);
        ch.solves[axis == 0 ? "in_plane_x" : "in_plane_y"] = info;
        Vec fx, fy, fz;
        s.face_flows(T, fx, fy, fz);
        for (size_t l = 0; l < nl; ++l) {
            double moment = 0.0;
            for (int k = slices[l].first; k < slices[l].second; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i) {
                        size_t c = g.idx(i, j, k);
                        if (axis == 0) moment += fx[c] * 0.5 * (g.dx(i) + g.dx((i + 1) % nx));
                        else moment += fy[c] * 0.5 * (g.dy(j) + g.dy((j + 1) % ny));
                    }
            (axis == 0 ? kxl : kyl)[l] = -moment;
        }
    }

    double sum_kx = 0, sum_ky = 0;
    for (size_t l = 0; l < nl; ++l) {
        LayerProps p;
        p.name = l < g.layer_names.size() ? g.layer_names[l] : "layer" + std::to_string(l);
        p.thickness = g.ze[slices[l].second] - g.ze[slices[l].first];
        p.kx = kxl[l] / volume[l];
        p.ky = kyl[l] / volume[l];
        p.R = tp.r_layer[l];
        p.kz = p.thickness / p.R;
        p.R_above = tp.r_int[l];
        p.rhoc = rhoc[l] / volume[l];
        p.kz_mix = kzmix[l] / volume[l];
        sum_kx += p.kx * p.thickness;
        sum_ky += p.ky * p.thickness;
        ch.layers.push_back(p);
    }
    ch.thickness = g.lz();
    ch.kx = sum_kx / ch.thickness;
    ch.ky = sum_ky / ch.thickness;
    ch.R_total = tp.r_total;
    ch.kz = ch.thickness / ch.R_total;
    return ch;
}

std::string Characterization::table() const {
    char buf[256];
    std::string out;
    std::snprintf(buf, sizeof buf, "%-24s%10s%9s%9s%9s%9s%13s%13s\n", "layer", "t [nm]", "kx", "ky", "kz",
                  "kz mix", "R [mm2K/W]", "R_int above");
    out += buf;
    for (const auto& p : layers) {
        std::snprintf(buf, sizeof buf, "%-24s%10.1f%9.3g%9.3g%9.3g%9.3g%13.4g%13.4g\n", p.name.c_str(),
                      p.thickness * 1e9, p.kx, p.ky, p.kz, p.kz_mix, p.R * 1e6, p.R_above * 1e6);
        out += buf;
    }
    std::snprintf(buf, sizeof buf, "%-24s%10.1f%9.3g%9.3g%9.3g%9s%13.4g", "TOTAL", thickness * 1e9, kx, ky, kz,
                  "", R_total * 1e6);
    return out + buf;
}

Operating operate(const Grid& g, const BC& bottom, const BC& top, double tol) {
    if (!(bottom.fixes_temperature() || top.fixes_temperature()))
        throw std::invalid_argument("the operating solve needs a temperature reference: make the top or "
                                    "bottom boundary 'dirichlet' or 'convection'");
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    System s = assemble(g, bottom, top);
    Operating op;
    op.T = LinearSolver(s.A, tol).solve(s.b, &op.info);
    Vec qb, qt, t_low, t_high;
    s.boundary_flows(op.T, qb, qt);
    for (double v : qb) op.q_bottom += v;
    for (double v : qt) op.q_top += v;
    const double area = g.lx() * g.ly();
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) op.power += g.q[g.idx(i, j, k)] * g.dx(i) * g.dy(j) * g.dz(k);
    for (const BC* bc : {&bottom, &top})
        if (bc->kind == BC::Flux) op.power += bc->q * area;
    double out = (bottom.fixes_temperature() ? op.q_bottom : 0.0) + (top.fixes_temperature() ? op.q_top : 0.0);
    op.imbalance = std::abs(op.power - out) / std::max(std::abs(op.power), 1e-300);

    size_t hot = size_t(std::max_element(op.T.begin(), op.T.end()) - op.T.begin());
    op.T_max = op.T[hot];
    int hk = int(hot / (size_t(nx) * ny)), hj = int((hot / nx) % ny), hi = int(hot % nx);
    op.hotspot = {0.5 * (g.xe[hi] + g.xe[hi + 1]), 0.5 * (g.ye[hj] + g.ye[hj + 1]),
                  0.5 * (g.ze[hk] + g.ze[hk + 1])};
    plane_temperatures(s, op.T, t_low, t_high);
    op.T_face_bottom = t_high[0];
    op.T_face_top = t_low[nz];
    for (auto [k0, k1] : g.layer_slices()) {
        double sum = 0, vol = 0, lo = kInf, hi2 = -kInf;
        for (int k = k0; k < k1; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    double t = op.T[g.idx(i, j, k)], v = g.dx(i) * g.dy(j) * g.dz(k);
                    sum += t * v;
                    vol += v;
                    lo = std::min(lo, t);
                    hi2 = std::max(hi2, t);
                }
        op.layer_T_mean.push_back(sum / vol);
        op.layer_T_min.push_back(lo);
        op.layer_T_max.push_back(hi2);
    }
    return op;
}

}  // namespace st
