// Atomistic model of amorphous silica for the Wigner conductivity engine.
//
// A periodic cell of a few hundred atoms is melted and quenched by molecular
// dynamics with the BKS pair potential, relaxed to a force-free minimum, and its
// harmonic and cubic force constants are assembled analytically.
//
// Potential: van Beest, Kramer and van Santen, Phys. Rev. Lett. 64, 1955 (1990),
// with the Coulomb term truncated and screened as in Carre et al., J. Chem. Phys.
// 127, 114512 (2007), and a short-range 30-6 repulsion that removes the
// unphysical collapse of the Buckingham form at small separations.
//
// Units: eV, angstrom, atomic mass units; time in angstrom * sqrt(amu / eV).
#include "lapack.hpp"
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace st {

namespace {
const double KB_EV = 8.617333262e-5, KE = 14.399645, TIME_FS = 10.1805057;
const double MASS[2] = {28.0855, 15.9994}, CHARGE[2] = {2.4, -1.2};
const double RC_COULOMB = 10.17, RC_SHORT = 5.5, GAMMA = 0.5;
// pair type = species_i + species_j: 0 Si-Si, 1 Si-O, 2 O-O    (A eV, b 1/A, C eV A^6, eps eV, sigma A)
const double BUCK[3][5] = {{0.0, 1.0, 0.0, 0.0, 1.0},
                           {18003.7572, 4.87318, 133.5381, 3.0976e-3, 1.313},
                           {1388.7730, 2.76000, 175.0000, 1.0514e-3, 1.779}};

// exp(-gamma^2 / (r - rc)^2) and its first three derivatives
void smooth(double r, double rc, double g[4]) {
    double u = std::min(r - rc, -1e-9), a = 2.0 * GAMMA * GAMMA / (u * u * u);
    double a1 = -6.0 * GAMMA * GAMMA / std::pow(u, 4), a2 = 24.0 * GAMMA * GAMMA / std::pow(u, 5);
    g[0] = std::exp(-GAMMA * GAMMA / (u * u));
    g[1] = g[0] * a;
    g[2] = g[0] * (a * a + a1);
    g[3] = g[0] * (a * a * a + 3.0 * a * a1 + a2);
}

void short_range(double r, int kind, double f[4]) {
    const double A = BUCK[kind][0], b = BUCK[kind][1], C = BUCK[kind][2], eps = BUCK[kind][3], sig = BUCK[kind][4];
    double e = A * std::exp(-b * r), s6 = std::pow(sig / r, 6), s30 = std::pow(sig / r, 30);
    f[0] = e - C / std::pow(r, 6) + 4 * eps * (s30 - s6);
    f[1] = -b * e + 6 * C / std::pow(r, 7) + 4 * eps * (-30 * s30 + 6 * s6) / r;
    f[2] = b * b * e - 42 * C / std::pow(r, 8) + 4 * eps * (930 * s30 - 42 * s6) / (r * r);
    f[3] = -b * b * b * e + 336 * C / std::pow(r, 9) + 4 * eps * (-29760 * s30 + 336 * s6) / (r * r * r);
}

// energy and -V'/r tabulated on a uniform grid in r^2: molecular dynamics spends
// its time here, and a lookup avoids every exp, power and square root
struct Table {
    double lo, scale;
    int n;
    Vec e, g;
    Table() : lo(0.7 * 0.7), n(200000), e(3 * 200000), g(3 * 200000) {
        double hi = (RC_COULOMB + 1.0) * (RC_COULOMB + 1.0);
        scale = (n - 1) / (hi - lo);
        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < n; ++i) {
                double r = std::sqrt(lo + i / scale), v, v1, v2;
                pair_potential(r, k, v, v1, v2);
                e[size_t(k) * n + i] = v;
                g[size_t(k) * n + i] = v1 / r;
            }
    }
};

void minimum_image(const Vec& x, int i, int j, double box, double d[3]) {
    for (int c = 0; c < 3; ++c) {
        d[c] = x[3 * j + c] - x[3 * i + c];
        d[c] -= box * std::nearbyint(d[c] / box);
    }
}
}  // namespace

void pair_potential(double r, int kind, double& v, double& v1, double& v2, double* v3) {
    v = v1 = v2 = 0.0;
    if (v3) *v3 = 0.0;
    if (r >= RC_COULOMB) return;
    const double qq = KE * (kind == 0 ? CHARGE[0] * CHARGE[0] : kind == 1 ? CHARGE[0] * CHARGE[1] : CHARGE[1] * CHARGE[1]);
    const double rc = RC_COULOMB;
    double c = qq * (1.0 / r - 1.0 / rc + (r - rc) / (rc * rc)), c1 = qq * (-1.0 / (r * r) + 1.0 / (rc * rc));
    double c2 = qq * 2.0 / (r * r * r), c3 = -qq * 6.0 / std::pow(r, 4), g[4];
    smooth(r, rc, g);
    v = c * g[0];
    v1 = c1 * g[0] + c * g[1];
    v2 = c2 * g[0] + 2 * c1 * g[1] + c * g[2];
    double t3 = c3 * g[0] + 3 * c2 * g[1] + 3 * c1 * g[2] + c * g[3];
    if (kind && r < RC_SHORT) {
        double f[4], f0[4], s[4];
        short_range(r, kind, f);
        short_range(RC_SHORT, kind, f0);
        f[0] -= f0[0];
        smooth(r, RC_SHORT, s);
        v += f[0] * s[0];
        v1 += f[1] * s[0] + f[0] * s[1];
        v2 += f[2] * s[0] + 2 * f[1] * s[1] + f[0] * s[2];
        t3 += f[3] * s[0] + 3 * f[2] * s[1] + 3 * f[1] * s[2] + f[0] * s[3];
    }
    if (v3) *v3 = t3;
}

double Glass::mass(int i) const { return MASS[species[i]]; }

double Glass::density() const {
    double m = 0;
    for (int i = 0; i < n(); ++i) m += mass(i);
    return m * 1.66053907 / (box * box * box);
}

void Glass::pairs(double cutoff, std::vector<int>& pi, std::vector<int>& pj) const {
    if (2 * cutoff > box) throw std::invalid_argument("cell too small for the minimum-image convention");
    pi.clear();
    pj.clear();
    const double c2 = cutoff * cutoff;
    double d[3];
    for (int i = 0; i < n(); ++i)
        for (int j = i + 1; j < n(); ++j) {
            minimum_image(pos, i, j, box, d);
            if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < c2) { pi.push_back(i); pj.push_back(j); }
        }
}

double Glass::energy_forces(const std::vector<int>& pi, const std::vector<int>& pj, const Vec& x, Vec& forces) const {
    forces.assign(3 * n(), 0.0);
    double energy = 0, d[3], v, v1, v2;
    for (size_t p = 0; p < pi.size(); ++p) {
        int i = pi[p], j = pj[p];
        minimum_image(x, i, j, box, d);
        double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        pair_potential(r, species[i] + species[j], v, v1, v2);
        energy += v;
        for (int c = 0; c < 3; ++c) {      // force on i points along +d when V' > 0
            forces[3 * i + c] += v1 / r * d[c];
            forces[3 * j + c] -= v1 / r * d[c];
        }
    }
    return energy;
}

void Glass::hessian(const std::vector<int>& pi, const std::vector<int>& pj, const Vec& x, Vec& h, Vec* dout,
                    Vec* blocks) const {
    const size_t m = size_t(3) * n();
    h.assign(m * m, 0.0);
    if (dout) dout->resize(3 * pi.size());
    if (blocks) blocks->resize(9 * pi.size());
    double d[3], v, v1, v2;
    for (size_t p = 0; p < pi.size(); ++p) {
        int i = pi[p], j = pj[p];
        minimum_image(x, i, j, box, d);
        double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        pair_potential(r, species[i] + species[j], v, v1, v2);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                double k = (v2 - v1 / r) * d[a] * d[b] / (r * r) + (a == b ? v1 / r : 0.0);
                h[(3 * i + a) + m * (3 * i + b)] += k;
                h[(3 * j + a) + m * (3 * j + b)] += k;
                h[(3 * i + a) + m * (3 * j + b)] -= k;
                h[(3 * j + a) + m * (3 * i + b)] -= k;
                if (blocks) (*blocks)[9 * p + 3 * a + b] = k;
            }
        if (dout) for (int c = 0; c < 3; ++c) (*dout)[3 * p + c] = d[c];
    }
}

// Each pair contributes a block on its two atoms; the energy depends on x_j - x_i,
// so an index on atom i carries a minus sign.
void Glass::third_order(double cutoff, std::vector<long long>& coords, Vec& data) const {
    std::vector<int> pi, pj;
    pairs(cutoff > 0 ? cutoff : RC_COULOMB, pi, pj);
    coords.clear();
    data.clear();
    coords.reserve(3 * 216 * pi.size());
    data.reserve(216 * pi.size());
    double d[3], v, v1, v2, v3;
    for (size_t p = 0; p < pi.size(); ++p) {
        int i = pi[p], j = pj[p];
        minimum_image(pos, i, j, box, d);
        double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), nn[3] = {d[0] / r, d[1] / r, d[2] / r};
        pair_potential(r, species[i] + species[j], v, v1, v2, &v3);
        const double bc = v2 / r - v1 / (r * r);
        const int atom[2] = {i, j};
        const double sign[2] = {-1.0, 1.0};
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                for (int c = 0; c < 2; ++c)
                    for (int al = 0; al < 3; ++al)
                        for (int be = 0; be < 3; ++be)
                            for (int ga = 0; ga < 3; ++ga) {
                                double t = (v3 - 3 * bc) * nn[al] * nn[be] * nn[ga] +
                                           bc * (nn[al] * (be == ga) + nn[be] * (al == ga) + nn[ga] * (al == be));
                                coords.push_back(3 * atom[a] + al);
                                coords.push_back(3 * atom[b] + be);
                                coords.push_back(3 * atom[c] + ga);
                                data.push_back(sign[a] * sign[b] * sign[c] * t);
                            }
    }
}

Coordination Glass::coordination(double cutoff) const {
    std::vector<int> pi, pj, count(n(), 0);
    pairs(cutoff, pi, pj);
    double d[3], sum = 0, sum2 = 0;
    int bonds = 0;
    for (size_t p = 0; p < pi.size(); ++p) {
        if (species[pi[p]] == species[pj[p]]) continue;
        ++count[pi[p]];
        ++count[pj[p]];
        minimum_image(pos, pi[p], pj[p], box, d);
        double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        sum += r;
        sum2 += r * r;
        ++bonds;
    }
    Coordination c;
    int si = 0, ox = 0, si4 = 0, o2 = 0;
    for (int i = 0; i < n(); ++i) {
        if (species[i] == 0) { ++si; si4 += count[i] == 4; }
        else { ++ox; o2 += count[i] == 2; }
    }
    c.si_four_fold = double(si4) / std::max(si, 1);
    c.o_two_fold = double(o2) / std::max(ox, 1);
    if (bonds) {
        c.si_o_bond = sum / bonds;
        c.si_o_bond_std = std::sqrt(std::max(0.0, sum2 / bonds - c.si_o_bond * c.si_o_bond));
    }
    return c;
}

// Idealised beta-cristobalite supercell (24 atoms per cubic cell), rescaled to
// the requested density.  Only a starting point for the melt.
Glass cristobalite(int cells, double density) {
    const double fcc[4][3] = {{0, 0, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0}};
    const double tet[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    std::vector<std::array<double, 3>> basis;
    std::vector<int> kinds;
    for (int s = 0; s < 2; ++s)
        for (auto& f : fcc) { basis.push_back({f[0] + 0.25 * s, f[1] + 0.25 * s, f[2] + 0.25 * s}); kinds.push_back(0); }
    for (auto& t : tet)
        for (auto& f : fcc) { basis.push_back({f[0] + 0.125 * t[0], f[1] + 0.125 * t[1], f[2] + 0.125 * t[2]}); kinds.push_back(1); }
    Glass g;
    for (int a = 0; a < cells; ++a)
        for (int b = 0; b < cells; ++b)
            for (int c = 0; c < cells; ++c)
                for (size_t k = 0; k < basis.size(); ++k) {
                    const double shift[3] = {double(a), double(b), double(c)};
                    for (int x = 0; x < 3; ++x) {
                        double f = basis[k][x] - std::floor(basis[k][x]);
                        g.pos.push_back((f + shift[x]) / cells);
                    }
                    g.species.push_back(kinds[k]);
                }
    double m = 0;
    for (int i = 0; i < g.n(); ++i) m += g.mass(i);
    g.box = std::cbrt(m * 1.66053907 / density);
    for (double& x : g.pos) x *= g.box;
    return g;
}

// velocity-Verlet molecular dynamics with a Berendsen thermostat
Glass melt_quench(const Glass& start, double t_melt, double t_end, double melt_ps, double quench_ps, double hold_ps,
                  double dt_fs, unsigned seed, Log log) {
    static const Table table;
    Glass w = start;
    const int n = w.n();
    const double dt = dt_fs / TIME_FS;
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Vec vel(3 * n), f(3 * n, 0.0);
    double mom[3] = {0, 0, 0}, mtot = 0;
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < 3; ++c) { vel[3 * i + c] = normal(rng) * std::sqrt(KB_EV * t_melt / w.mass(i)); mom[c] += vel[3 * i + c] * w.mass(i); }
        mtot += w.mass(i);
    }
    for (int i = 0; i < n; ++i) for (int c = 0; c < 3; ++c) vel[3 * i + c] -= mom[c] / mtot;
    const long steps[3] = {std::lround(melt_ps * 1000 / dt_fs), std::lround(quench_ps * 1000 / dt_fs),
                           std::lround(hold_ps * 1000 / dt_fs)};
    const long total = steps[0] + steps[1] + steps[2];
    // the list must stay inside the minimum-image radius; a thin skin is rebuilt often
    const double skin = std::min(1.0, 0.5 * w.box - RC_COULOMB - 1e-3);
    const long every = std::max(1L, long(8 * skin));
    std::vector<int> pi, pj;
    double energy = 0;
    auto forces = [&] {
        std::fill(f.begin(), f.end(), 0.0);
        energy = 0;
        double d[3];
        for (size_t p = 0; p < pi.size(); ++p) {
            int i = pi[p], j = pj[p];
            minimum_image(w.pos, i, j, w.box, d);
            double t = std::clamp((d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - table.lo) * table.scale, 0.0, table.n - 1.001);
            size_t k = size_t(t), o = size_t(w.species[i] + w.species[j]) * table.n + k;
            double frac = t - double(k), g = table.g[o] * (1 - frac) + table.g[o + 1] * frac;
            energy += table.e[o] * (1 - frac) + table.e[o + 1] * frac;
            for (int c = 0; c < 3; ++c) { f[3 * i + c] += g * d[c]; f[3 * j + c] -= g * d[c]; }
        }
    };
    w.pairs(RC_COULOMB + skin, pi, pj);
    forces();
    for (long step = 0; step < total; ++step) {
        double target = step < steps[0] ? t_melt
                        : step < steps[0] + steps[1] ? t_melt + (t_end - t_melt) * double(step - steps[0]) / steps[1] : t_end;
        for (int i = 0; i < n; ++i)
            for (int c = 0; c < 3; ++c) {
                vel[3 * i + c] += 0.5 * dt * f[3 * i + c] / w.mass(i);
                w.pos[3 * i + c] += dt * vel[3 * i + c];
            }
        if (step % every == 0) {
            for (double& x : w.pos) x -= w.box * std::floor(x / w.box);
            w.pairs(RC_COULOMB + skin, pi, pj);
        }
        forces();
        double kinetic = 0;
        for (int i = 0; i < n; ++i)
            for (int c = 0; c < 3; ++c) {
                vel[3 * i + c] += 0.5 * dt * f[3 * i + c] / w.mass(i);
                kinetic += w.mass(i) * vel[3 * i + c] * vel[3 * i + c];
            }
        double temp = kinetic / (3 * n * KB_EV), lambda = std::sqrt(1.0 + dt_fs / 200.0 * (target / temp - 1.0));
        for (double& v : vel) v *= lambda;
        if (log && step % 5000 == 0) {
            char buf[160];
            std::snprintf(buf, sizeof buf, "step %6ld/%ld  T %7.0f K (target %5.0f)  E %8.4f eV/atom", step, total, temp,
                          target, energy / n);
            log(buf);
        }
    }
    for (double& x : w.pos) x -= w.box * std::floor(x / w.box);
    return w;
}

// Minimise the energy: damped dynamics (FIRE), then Newton steps on the exact
// Hessian until the forces vanish to rounding error, which a clean harmonic
// spectrum needs.
Glass relax(const Glass& start, Log log) {
    Glass w = start;
    const int n = w.n(), m = 3 * n;
    std::vector<int> pi, pj;
    w.pairs(std::min(RC_COULOMB + 1.5, 0.5 * w.box - 1e-3), pi, pj);
    Vec f, vel(m, 0.0);
    double dt = 0.1, alpha = 0.1;
    int positive = 0;
    w.energy_forces(pi, pj, w.pos, f);
    for (int it = 0; it < 20000; ++it) {
        double fmax = 0, power = 0, fn = 0, vn = 0;
        for (int i = 0; i < m; ++i) { fmax = std::max(fmax, std::abs(f[i])); power += f[i] * vel[i]; fn += f[i] * f[i]; vn += vel[i] * vel[i]; }
        if (fmax < 1e-5) break;
        if (power > 0) {
            double mix = alpha * std::sqrt(vn / std::max(fn, 1e-300));
            for (int i = 0; i < m; ++i) vel[i] = (1 - alpha) * vel[i] + mix * f[i];
            if (++positive > 5) { dt = std::min(dt * 1.1, 0.6); alpha *= 0.99; }
        } else {
            std::fill(vel.begin(), vel.end(), 0.0);
            dt *= 0.5; alpha = 0.1; positive = 0;
        }
        for (int i = 0; i < n; ++i)
            for (int c = 0; c < 3; ++c) { vel[3 * i + c] += dt * f[3 * i + c] / w.mass(i); w.pos[3 * i + c] += dt * vel[3 * i + c]; }
        w.energy_forces(pi, pj, w.pos, f);
    }
    for (int it = 0; it < 8; ++it) {
        w.pairs(RC_COULOMB, pi, pj);
        w.energy_forces(pi, pj, w.pos, f);
        double fmax = 0;
        for (double v : f) fmax = std::max(fmax, std::abs(v));
        if (log) { char buf[80]; std::snprintf(buf, sizeof buf, "newton %d: max force %.2e eV/A", it, fmax); log(buf); }
        if (fmax < 1e-10) break;
        Vec h, eig;
        w.hessian(pi, pj, w.pos, h);
        symmetric_eigen(m, h, eig);
        Vec proj(m, 0.0);                       // V^T f, scaled by 1/eigenvalue on the vibrational modes
        for (int s = 0; s < m; ++s) {
            if (eig[s] <= 1e-6) continue;       // drop the three translations
            double dot = 0;
            for (int i = 0; i < m; ++i) dot += h[size_t(i) + size_t(m) * s] * f[i];
            proj[s] = dot / eig[s];
        }
        for (int i = 0; i < m; ++i) {
            double step = 0;
            for (int s = 0; s < m; ++s) step += h[size_t(i) + size_t(m) * s] * proj[s];
            w.pos[i] += step;
        }
    }
    return w;
}

Glass load_glass(const std::string& name) {
    const std::string path = data_dir() + "/glass/" + name + ".glass";
    std::ifstream f(path);
    if (!f) throw std::invalid_argument("no stored glass model '" + name + "' (" + path + ")");
    Glass g;
    int n;
    f >> n >> g.box;
    std::string line;
    std::getline(f, line);
    std::getline(f, line);
    if (line.size() > 2 && line[0] == '#') g.info = line.substr(2);
    g.pos.resize(3 * n);
    g.species.resize(n);
    for (int i = 0; i < n; ++i) f >> g.species[i] >> g.pos[3 * i] >> g.pos[3 * i + 1] >> g.pos[3 * i + 2];
    if (!f) throw std::runtime_error("glass file '" + path + "' is truncated");
    return g;
}

void save_glass(const Glass& g, const std::string& name) {
    std::ofstream f(data_dir() + "/glass/" + name + ".glass");
    f.precision(17);
    f << g.n() << ' ' << g.box << "\n# " << (g.info.empty() ? "{}" : g.info) << "\n";
    for (int i = 0; i < g.n(); ++i)
        f << g.species[i] << ' ' << g.pos[3 * i] << ' ' << g.pos[3 * i + 1] << ' ' << g.pos[3 * i + 2] << '\n';
}

}  // namespace st
