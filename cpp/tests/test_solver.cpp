// Solver verification against closed-form solutions.
#include "stacktherm/stacktherm.hpp"

#include <cstdio>
#include <cstdlib>
#include <random>

using namespace st;

static int failures = 0;

static void expect_close(const char* what, double got, double want, double rel) {
    bool ok = std::abs(got - want) <= rel * std::abs(want) + 1e-300;
    std::printf("%s %-58s got %.10g  want %.10g\n", ok ? "ok  " : "FAIL", what, got, want);
    if (!ok) ++failures;
}
static void expect(const char* what, bool ok) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok) ++failures;
}

static Vec uniform(double length, int n) {
    Vec e(n + 1);
    for (int i = 0; i <= n; ++i) e[i] = length * i / n;
    return e;
}

int main() {
    {   // three films with a Kapitza resistance between each: R = sum(t/k) + sum(R_int)
        const double t[3] = {50e-9, 200e-9, 30e-9}, k[3] = {1.4, 148.0, 0.3}, r_int[2] = {2e-9, 5e-9};
        const int per = 4;
        Grid g;
        g.xe = uniform(1e-7, 3);
        g.ye = uniform(1e-7, 3);
        g.ze = {0.0};
        for (int l = 0; l < 3; ++l)
            for (int c = 0; c < per; ++c) { g.ze.push_back(g.ze.back() + t[l] / per); g.layer_of_k.push_back(l); }
        g.allocate();
        for (int kz = 0; kz < g.nz(); ++kz)
            for (int j = 0; j < 3; ++j)
                for (int i = 0; i < 3; ++i) {
                    size_t c = g.idx(i, j, kz);
                    g.kx[c] = g.ky[c] = g.kz[c] = k[kz / per];
                    if (kz == per - 1) g.rz[c] = r_int[0];
                    if (kz == 2 * per - 1) g.rz[c] = r_int[1];
                }
        Characterization c = characterize(g);
        expect_close("series stack: total resistance", c.R_total,
                     t[0] / k[0] + t[1] / k[1] + t[2] / k[2] + r_int[0] + r_int[1], 1e-9);
        for (int l = 0; l < 3; ++l) expect_close("series stack: layer resistance", c.layers[l].R, t[l] / k[l], 1e-9);
        expect_close("series stack: interface 0", c.layers[0].R_above, r_int[0], 1e-8);
        expect_close("series stack: interface 1", c.layers[1].R_above, r_int[1], 1e-8);
        expect_close("series stack: in-plane k of layer 1", c.layers[1].kx, k[1], 1e-9);
    }
    {   // laminate: lines along x.  kxx arithmetic, kyy harmonic incl. two interfaces per period
        const double pitch = 100e-9, width = 40e-9, km = 300.0, kd = 0.3, r_int = 1e-8;
        const int ny = 50;
        Grid g;
        g.xe = uniform(pitch, 2);
        g.ye = uniform(pitch, ny);
        g.ze = uniform(50e-9, 4);
        g.allocate();
        std::vector<bool> metal(ny);
        for (int j = 0; j < ny; ++j) { double yc = 0.5 * (g.ye[j] + g.ye[j + 1]); metal[j] = yc > 30e-9 && yc < 70e-9; }
        for (int kz = 0; kz < 4; ++kz)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < 2; ++i) {
                    size_t c = g.idx(i, j, kz);
                    g.kx[c] = g.ky[c] = g.kz[c] = metal[j] ? km : kd;
                    if (metal[j] != metal[(j + 1) % ny]) g.ry[c] = r_int;
                }
        Characterization c = characterize(g);
        const double f = width / pitch;
        expect_close("laminate: kx (arithmetic mean)", c.kx, f * km + (1 - f) * kd, 1e-9);
        expect_close("laminate: kz (arithmetic mean)", c.kz, f * km + (1 - f) * kd, 1e-9);
        expect_close("laminate: ky (harmonic mean with interfaces)", c.ky,
                     pitch / (width / km + (pitch - width) / kd + 2 * r_int), 1e-9);
    }
    {   // uniform source, adiabatic bottom, convection on top: exact 1D solution and energy balance
        const double t[2] = {100e-9, 400e-9}, k[2] = {10.0, 2.0}, h = 2e5, flux = 1e6;
        const int per = 20;
        Grid g;
        g.xe = uniform(1e-7, 2);
        g.ye = uniform(1e-7, 2);
        g.ze = {0.0};
        for (int l = 0; l < 2; ++l)
            for (int c = 0; c < per; ++c) { g.ze.push_back(g.ze.back() + t[l] / per); g.layer_of_k.push_back(l); }
        g.allocate();
        for (int kz = 0; kz < g.nz(); ++kz)
            for (int j = 0; j < 2; ++j)
                for (int i = 0; i < 2; ++i) {
                    size_t c = g.idx(i, j, kz);
                    g.kx[c] = g.ky[c] = g.kz[c] = k[kz / per];
                    if (kz < per) g.q[c] = flux / t[0];
                }
        Operating op = operate(g, BC(), BC::convection(h, 300.0));
        expect("operating: energy balance closes", op.imbalance < 1e-8);
        expect_close("operating: hottest temperature", op.T_max,
                     300.0 + flux / h + flux * t[1] / k[1] + flux * t[0] / (2 * k[0]), 2e-3);
    }
    {   // random heterogeneous cell with 50:1 aspect ratio: residual and symmetry of the result
        const int n = 24;
        std::mt19937 rng(0);
        std::uniform_real_distribution<double> u(-1.0, 2.5);
        Grid g;
        g.xe = uniform(1e-6, n);
        g.ye = uniform(1e-6, n);
        g.ze = uniform(2e-8, n);
        g.allocate();
        for (size_t c = 0; c < g.n(); ++c) g.kx[c] = g.ky[c] = g.kz[c] = std::pow(10.0, u(rng));
        System s = assemble(g, BC::dirichlet(1.0), BC::dirichlet(0.0), true, true, 0, 0, false);
        SolveInfo info;
        Vec T = LinearSolver(s.A).solve(s.b, &info), r;
        s.A.mul(T, r);
        double rn = 0, bn = 0;
        for (size_t c = 0; c < r.size(); ++c) { rn += (r[c] - s.b[c]) * (r[c] - s.b[c]); bn += s.b[c] * s.b[c]; }
        std::printf("     heterogeneous solve: %d iterations, setup %.3f s, solve %.3f s\n", info.iterations,
                    info.setup_s, info.solve_s);
        expect("heterogeneous solve: residual below 1e-9", std::sqrt(rn / bn) < 1e-9);
        Vec qb, qt;
        s.boundary_flows(T, qb, qt);
        double in = 0, out = 0;
        for (double v : qb) in -= v;
        for (double v : qt) out += v;
        expect_close("heterogeneous solve: flow in equals flow out", in, out, 1e-7);
    }
    expect_close("graded edges end at the thickness", graded_edges(1e-5, 1e-8, 2e-8, 1.3).back(), 1e-5, 1e-15);
    expect_close("E3(0.3) exponential integral", expint(3, 0.3), 0.300041826564, 1e-10);
    expect_close("E5(2.5) exponential integral", expint(5, 2.5), 0.0119073798263, 1e-10);
    expect_close("Voigt profile reduces to a Gaussian", voigt(0.4, 1.0, 1e-12),
                 std::exp(-0.08) / std::sqrt(2 * kPi), 1e-9);
    expect_close("Voigt profile at the centre (sigma = gamma = 1)", voigt(0.0, 1.0, 1.0), 0.20870928052, 1e-10);
    expect_close("Voigt profile off centre", voigt(1.3, 0.7, 0.4), 0.13890764492, 1e-10);
    expect_close("unit parser: 42 nm", parse_quantity("42 nm"), 42e-9, 1e-12);
    expect_close("unit parser: conductance to resistance", parse_interface("150 MW/m2K"), 1.0 / 1.5e8, 1e-12);
    std::printf("\n%s\n", failures ? "FAILED" : "all passed");
    return failures ? 1 : 0;
}
