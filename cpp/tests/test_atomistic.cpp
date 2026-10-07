// Atomistic glass model and Wigner conductivity engine.
#include "stacktherm/stacktherm.hpp"

#include <cstdio>
#include <random>

using namespace st;

static int failures = 0;
static void expect_close(const char* what, double got, double want, double rel) {
    bool ok = std::abs(got - want) <= rel * std::abs(want) + 1e-300;
    std::printf("%s %-62s got %.10g  want %.10g\n", ok ? "ok  " : "FAIL", what, got, want);
    if (!ok) ++failures;
}
static void expect(const char* what, bool ok) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok) ++failures;
}

int main() {
    const double CM1 = 2.0 * kPi * 2.99792458e10;
    {   // pair potential: reference values and derivatives by finite differences
        double v, v1, v2, v3;
        pair_potential(1.61, 1, v, v1, v2, &v3);
        expect_close("Si-O at 1.61 A: energy", v, -18.79588621, 1e-9);
        expect_close("Si-O at 1.61 A: first derivative", v1, 9.907819954, 1e-9);
        expect_close("Si-O at 1.61 A: second derivative", v2, 22.61975171, 1e-9);
        expect_close("Si-O at 1.61 A: third derivative", v3, -158.8187448, 1e-9);
        double worst = 0;
        for (int kind = 0; kind < 3; ++kind)
            for (double r = 1.2; r < 10.0; r += 0.37) {
                const double h = 1e-5;
                double a[4], b[4], c[4], t;
                pair_potential(r + h, kind, a[0], a[1], a[2], &a[3]);
                pair_potential(r - h, kind, b[0], b[1], b[2], &b[3]);
                pair_potential(r, kind, c[0], c[1], c[2], &c[3]);
                (void)t;
                for (int d = 1; d < 4; ++d)
                    worst = std::max(worst, std::abs((a[d - 1] - b[d - 1]) / (2 * h) - c[d]) / (std::abs(c[d]) + 1e-4));
            }
        expect("all three derivatives match finite differences", worst < 1e-5);
        pair_potential(0.9, 1, v, v1, v2);
        expect("no Buckingham collapse at short range", v > 50.0);
        pair_potential(10.169, 2, v, v1, v2);
        expect("potential vanishes smoothly at the cut-off", std::abs(v) < 1e-12 && std::abs(v1) < 1e-10);
    }

    Glass start = cristobalite(3);
    Coordination c0 = start.coordination();
    expect("cristobalite start: 648 atoms at 2.20 g/cm3", start.n() == 648 && std::abs(start.density() - 2.20) < 1e-9);
    expect("cristobalite start: perfect network", c0.si_four_fold == 1.0 && c0.o_two_fold == 1.0);

    Glass g = load_glass("a-SiO2-648");
    std::vector<int> pi, pj;
    g.pairs(10.17, pi, pj);
    Vec f;
    g.energy_forces(pi, pj, g.pos, f);
    double fmax = 0;
    for (double v : f) fmax = std::max(fmax, std::abs(v));
    Coordination c = g.coordination();
    expect("stored glass is at a true minimum (forces below 1e-8 eV/A)", fmax < 1e-8);
    expect("stored glass is a four-fold / two-fold network", c.si_four_fold > 0.9 && c.o_two_fold > 0.9);
    expect("stored glass is disordered (Si-O bond spread)", c.si_o_bond_std > 0.01 && c.si_o_bond > 1.58 && c.si_o_bond < 1.66);

    {   // Hessian against finite differences of the forces along a random direction
        const int m = 3 * g.n();
        Vec h, step(m), xp(g.pos), xm(g.pos), fp, fm;
        g.hessian(pi, pj, g.pos, h);
        std::mt19937 rng(0);
        std::normal_distribution<double> normal;
        for (double& s : step) s = normal(rng);
        const double eps = 1e-5;
        for (int i = 0; i < m; ++i) { xp[i] += eps * step[i]; xm[i] -= eps * step[i]; }
        g.energy_forces(pi, pj, xp, fp);
        g.energy_forces(pi, pj, xm, fm);
        double worst = 0, scale = 0;
        for (int r = 0; r < m; ++r) {
            double hv = 0;
            for (int col = 0; col < m; ++col) hv += h[size_t(r) + size_t(m) * col] * step[col];
            worst = std::max(worst, std::abs(-(fp[r] - fm[r]) / (2 * eps) - hv));
            scale = std::max(scale, std::abs(hv));
        }
        expect("Hessian matches finite-difference forces", worst < 1e-5 * scale);
    }

    Modes md = glass_modes(g);
    expect("three translations, every other mode stable", md.omega[2] / CM1 < 0.05 && md.omega[3] / CM1 > 5.0);
    expect_close("lowest vibrational mode (cm-1)", md.omega[3] / CM1, 33.37947683, 1e-6);
    expect_close("highest vibrational mode (cm-1)", md.omega[md.m - 1] / CM1, 1293.093976, 1e-8);

    {   // cubic force constants: analytic tensor against finite differences, and permutation symmetry
        const int m = md.m, s = 700;
        std::vector<long long> coords;
        Vec data, analytic(size_t(m) * m, 0.0), u(m);
        g.third_order(-1.0, coords, data);
        for (int i = 0; i < m; ++i) u[i] = md.vec[size_t(i) + size_t(m) * s] / std::sqrt(g.mass(i / 3));
        for (size_t e = 0; e < data.size(); ++e)
            analytic[size_t(coords[3 * e]) + size_t(m) * coords[3 * e + 1]] += data[e] * u[coords[3 * e + 2]];
        double largest = 0;
        for (double v : u) largest = std::max(largest, std::abs(v));
        const double delta = 0.01 / largest;
        Vec xp(g.pos), xm(g.pos), hp, hm;
        for (int i = 0; i < m; ++i) { xp[i] += delta * u[i]; xm[i] -= delta * u[i]; }
        g.hessian(pi, pj, xp, hp);
        g.hessian(pi, pj, xm, hm);
        double worst = 0, scale = 0;
        for (size_t k = 0; k < analytic.size(); ++k) {
            double numeric = (hp[k] - hm[k]) / (2 * delta);
            worst = std::max(worst, std::abs(numeric - analytic[k]));
            scale = std::max(scale, std::abs(numeric));
        }
        expect("analytic third-order tensor matches finite differences", worst < 1e-3 * scale);
        Vec va = cubic_row(md, 300, pi, pj), vb = cubic_row(md, 900, pi, pj);
        expect_close("cubic constants are permutation symmetric", va[size_t(900) + size_t(m) * 1500],
                     vb[size_t(300) + size_t(m) * 1500], 5e-3);
    }

    Vec diff;
    double k300 = wigner_conductivity(md, 300.0, nullptr, -1.0, &diff), k100 = wigner_conductivity(md, 100.0);
    expect_close("harmonic limit at 300 K (reference)", k300, 1.165481109, 1e-7);
    expect_close("harmonic limit at 100 K (reference)", k100, 0.4339640893, 1e-7);
    expect("glass-like: conductivity rises with temperature", k100 < k300);
    Vec small(md.m, 0.2 * CM1);
    expect_close("a small uniform linewidth changes little", wigner_conductivity(md, 300.0, &small), k300, 0.05);

    std::vector<int> chosen;
    Vec gam;
    linewidths(md, {100.0, 300.0}, 40, 3.0, chosen, gam);
    double m100 = 0, m300 = 0;
    for (size_t i = 0; i < chosen.size(); ++i) { m100 += gam[i]; m300 += gam[chosen.size() + i]; }
    expect("linewidths broaden with temperature", m300 > 1.5 * m100);
    Vec g300 = smooth_linewidths(md, chosen, Vec(gam.begin() + chosen.size(), gam.end()));
    double kw = wigner_conductivity(md, 300.0, &g300);
    expect_close("Wigner conductivity at 300 K against the measured 1.38 W/m/K", kw, 1.38, 0.10);
    double lo = wigner_conductivity(md, 300.0, &g300, 0.5), hi = wigner_conductivity(md, 300.0, &g300, 8.0);
    expect("anharmonic linewidths remove the dependence on the regularisation", hi / lo - 1.0 < 0.05);
    std::printf("\n%s\n", failures ? "FAILED" : "all passed");
    return failures ? 1 : 0;
}
