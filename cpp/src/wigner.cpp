// Wigner thermal conductivity of a disordered solid from its atomistic model.
//
// The Wigner formulation of heat transport (Simoncelli, Marzari and Mauri, Nat.
// Phys. 15, 809 (2019); Phys. Rev. X 12, 041011 (2022)) gives the conductivity of
// crystals and glasses from one expression.  For the vibrational modes s of a
// periodic cell of volume V it reads
//
//     k = 1/V * sum_{s,s'} (w_s + w_s')/4 * (C_s/w_s + C_s'/w_s') * |v_ss'|^2 / 3
//               * pi * L(w_s - w_s'; (G_s + G_s')/2)
//
// with C the modal heat capacity, v the velocity operator, G the anharmonic
// linewidths and L a unit-area Lorentzian of that half-width.  Off-diagonal terms
// are wave-like tunnelling between modes, which is what carries heat in a glass;
// with G -> 0 the expression reduces to the Allen-Feldman theory.
//
// The modes come from the exact Hessian of a pair-potential glass at the zone
// centre, and the linewidths from three-phonon scattering with cubic force
// constants obtained by finite differences of the Hessian along each mode.  A
// Gaussian of width eta is folded into the Lorentzian (a Voigt profile) to
// regularise the discrete spectrum of the finite cell, as in Simoncelli, Mauri
// and Marzari, npj Comput. Mater. 9, 106 (2023).
#include "lapack.hpp"
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace st {

namespace {
const double HBAR = 1.054571817e-34, KB = 1.380649e-23, EV = 1.602176634e-19, AMU = 1.66053907e-27, ANG = 1e-10;
const double CM1 = 2.0 * kPi * 2.99792458e10;        // rad/s per cm^-1
const double D_UNIT = EV / (ANG * ANG * AMU);          // (eV/A^2/amu) -> s^-2
const double RC = 10.17;

double bose(double omega, double T) { return 1.0 / std::expm1(std::clamp(HBAR * omega / (KB * T), 1e-12, 700.0)); }
double heat_capacity(double omega, double T) {
    double x = std::clamp(HBAR * omega / (KB * T), 1e-12, 700.0), em = std::expm1(x);
    return KB * x * x * std::exp(x) / (em * em);
}
Vec inverse_sqrt_mass(const Glass& g) {
    Vec inv(3 * g.n());
    for (int i = 0; i < g.n(); ++i) inv[3 * i] = inv[3 * i + 1] = inv[3 * i + 2] = 1.0 / std::sqrt(g.mass(i));
    return inv;
}
}  // namespace

Modes glass_modes(const Glass& g) {
    Modes md;
    md.glass = &g;
    const int n = g.n(), m = 3 * n;
    md.m = m;
    std::vector<int> pi, pj;
    g.pairs(RC, pi, pj);
    Vec h, d, blocks;
    g.hessian(pi, pj, g.pos, h, &d, &blocks);
    const Vec inv = inverse_sqrt_mass(g);
    for (int c = 0; c < m; ++c)
        for (int r = 0; r < m; ++r) h[size_t(r) + size_t(m) * c] *= inv[r] * inv[c];
    Vec w;
    symmetric_eigen(m, h, w);
    if (w[3] <= 0) throw std::runtime_error("the structure is not at a minimum: negative curvature found");
    md.vec = std::move(h);
    md.omega.resize(m);
    for (int s = 0; s < m; ++s) md.omega[s] = std::sqrt(std::max(w[s], 0.0) * D_UNIT);
    md.volume = g.box * g.box * g.box * ANG * ANG * ANG;

    // velocity operator between every pair of modes: v_ss' = e_s^T G e_s' / (2 sqrt(w_s w_s')),
    // G_ij = d_ij * Phi_ij / sqrt(m_i m_j), antisymmetric in the atom pair
    md.vel2.assign(size_t(m) * m, 0.0);
    Vec gm, tmp, a;
    for (int axis = 0; axis < 3; ++axis) {
        gm.assign(size_t(m) * m, 0.0);
        for (size_t p = 0; p < pi.size(); ++p) {
            const int i = pi[p], j = pj[p];
            const double scale = -d[3 * p + axis] * inv[3 * i] * inv[3 * j];
            for (int al = 0; al < 3; ++al)
                for (int be = 0; be < 3; ++be) {
                    double e = scale * blocks[9 * p + 3 * al + be];
                    gm[size_t(3 * i + al) + size_t(m) * (3 * j + be)] += e;
                    gm[size_t(3 * j + al) + size_t(m) * (3 * i + be)] -= e;
                }
        }
        matmul(m, 'N', gm, 'N', md.vec, tmp);
        matmul(m, 'T', md.vec, 'N', tmp, a);
        for (size_t k = 0; k < a.size(); ++k) md.vel2[k] += a[k] * a[k];
    }
    const double unit = (ANG * D_UNIT) * (ANG * D_UNIT) / 12.0;
    for (int c = 0; c < m; ++c)
        for (int r = 0; r < m; ++r) {
            double& v = md.vel2[size_t(r) + size_t(m) * c];
            v = (r < 3 || c < 3 || r == c) ? 0.0 : v * unit / (md.omega[r] * md.omega[c]);
        }
    return md;
}

// Third derivative of the energy with respect to mode s and every pair of
// modes, by central differences of the Hessian along the mode (SI units).
Vec cubic_row(const Modes& md, int s, const std::vector<int>& pi, const std::vector<int>& pj, double amplitude) {
    const Glass& g = *md.glass;
    const int m = md.m;
    const Vec inv = inverse_sqrt_mass(g);
    Vec direction(m);
    double largest = 0;
    for (int i = 0; i < m; ++i) { direction[i] = md.vec[size_t(i) + size_t(m) * s] * inv[i]; largest = std::max(largest, std::abs(direction[i])); }
    const double delta = amplitude / largest;          // angstrom * sqrt(amu)
    Vec xp(g.pos), xm(g.pos), hp, hm, tmp, out;
    for (int i = 0; i < m; ++i) { xp[i] += delta * direction[i]; xm[i] -= delta * direction[i]; }
    g.hessian(pi, pj, xp, hp);
    g.hessian(pi, pj, xm, hm);
    for (int c = 0; c < m; ++c)
        for (int r = 0; r < m; ++r) {
            size_t k = size_t(r) + size_t(m) * c;
            hp[k] = (hp[k] - hm[k]) * inv[r] * inv[c] / (2.0 * delta);
        }
    matmul(m, 'N', hp, 'N', md.vec, tmp);
    matmul(m, 'T', md.vec, 'N', tmp, out);
    const double unit = D_UNIT / (ANG * std::sqrt(AMU));
    for (double& v : out) v *= unit;
    return out;
}

// Three-phonon linewidths (full width, rad/s) from Fermi's golden rule:
//   G_s = pi hbar / (8 w_s) * sum |V_ss's''|^2 / (w' w'')
//         * [(n' + n'' + 1) delta(w_s - w' - w'') + 2 (n' - n'') delta(w_s + w' - w'')]
void linewidths(const Modes& md, const Vec& temperatures, int sample, double sigma_cm, std::vector<int>& chosen,
                Vec& gamma, Log log) {
    const int m = md.m, active = m - 3;
    chosen.clear();
    for (int c = 0; c < sample; ++c) {
        int k = 3 + int(std::nearbyint(double(c) * (active - 1) / std::max(sample - 1, 1)));
        if (chosen.empty() || chosen.back() != k) chosen.push_back(k);
    }
    std::vector<int> pi, pj;
    md.glass->pairs(RC, pi, pj);
    const double sigma = sigma_cm * CM1, norm = 1.0 / (sigma * std::sqrt(2.0 * kPi));
    const size_t nt = temperatures.size(), nc = chosen.size();
    gamma.assign(nt * nc, 0.0);
    Vec wd(size_t(m) * m), wa(size_t(m) * m), occ(m);
    for (size_t c = 0; c < nc; ++c) {
        const int s = chosen[c];
        Vec v3 = cubic_row(md, s, pi, pj);
        for (int b = 3; b < m; ++b)
            for (int a = 3; a < m; ++a) {
                size_t k = size_t(a) + size_t(m) * b;
                double wgt = v3[k] * v3[k] / (md.omega[a] * md.omega[b]);
                double x1 = (md.omega[s] - md.omega[a] - md.omega[b]) / sigma, x2 = (md.omega[s] + md.omega[a] - md.omega[b]) / sigma;
                wd[k] = std::abs(x1) < 8 ? wgt * std::exp(-0.5 * x1 * x1) * norm : 0.0;
                wa[k] = std::abs(x2) < 8 ? wgt * std::exp(-0.5 * x2 * x2) * norm : 0.0;
            }
        for (size_t t = 0; t < nt; ++t) {
            for (int a = 3; a < m; ++a) occ[a] = bose(md.omega[a], temperatures[t]);
            double total = 0;
            for (int b = 3; b < m; ++b)
                for (int a = 3; a < m; ++a) {
                    size_t k = size_t(a) + size_t(m) * b;
                    total += wd[k] * (occ[a] + occ[b] + 1.0) + 2.0 * wa[k] * (occ[a] - occ[b]);
                }
            gamma[t * nc + c] = kPi * HBAR / (8.0 * md.omega[s]) * total;
        }
        if (log && c % 20 == 0) {
            char buf[160];
            std::snprintf(buf, sizeof buf, "linewidth %zu/%zu: %7.1f cm-1 -> %.2f cm-1 at %.0f K", c + 1, nc,
                          md.omega[s] / CM1, gamma[(nt - 1) * nc + c] / CM1, temperatures.back());
            log(buf);
        }
    }
}

// linewidth of every mode from the sampled ones, by a Gaussian kernel in frequency
Vec smooth_linewidths(const Modes& md, const std::vector<int>& chosen, const Vec& values, double width_cm) {
    Vec out(md.m, 0.0);
    const double width = width_cm * CM1;
    for (int s = 3; s < md.m; ++s) {
        double num = 0, den = 0;
        for (size_t c = 0; c < chosen.size(); ++c) {
            double x = (md.omega[s] - md.omega[chosen[c]]) / width, w = std::exp(-0.5 * x * x);
            num += w * values[c];
            den += w;
        }
        out[s] = num / std::max(den, 1e-300);
    }
    return out;
}

double wigner_conductivity(const Modes& md, double T, const Vec* gamma, double eta_cm, Vec* diffusivity) {
    const int m = md.m;
    const double eta = eta_cm >= 0 ? eta_cm * CM1 : 2.0 * md.level_spacing();
    Vec heat(m, 0.0), diff(m, 0.0);
    for (int s = 3; s < m; ++s) heat[s] = heat_capacity(md.omega[s], T);
    double kappa = 0;
    for (int b = 3; b < m; ++b)
        for (int a = 3; a < m; ++a) {
            if (a == b) continue;
            double v2 = md.vel2[size_t(a) + size_t(m) * b];
            if (v2 == 0.0) continue;
            double half = gamma ? 0.5 * ((*gamma)[a] + (*gamma)[b]) : 0.0;
            double kernel = v2 * kPi * voigt(md.omega[a] - md.omega[b], eta, half);
            diff[a] += (md.omega[a] + md.omega[b]) / (2.0 * md.omega[a]) * kernel;
            kappa += 0.25 * (md.omega[a] + md.omega[b]) * heat[a] / md.omega[a] * kernel;
        }
    if (diffusivity) *diffusivity = diff;
    return 2.0 * kappa / md.volume;
}

std::string wigner_compute(const std::string& name, const Vec& temps, int sample, Log log) {
    Glass glass = load_glass(name);
    Modes md = glass_modes(glass);
    const int m = md.m;
    std::vector<int> chosen;
    Vec gam;
    linewidths(md, temps, sample, 3.0, chosen, gam, log);
    const size_t nt = temps.size(), nc = chosen.size();
    std::ostringstream o;
    o.precision(10);
    o << "{\"name\":\"" << name << "\",\"backend\":\"native\",\"structure\":" << (glass.info.empty() ? "{}" : glass.info)
      << ",\"modes\":" << m - 3 << ",\"level_spacing_cm\":" << md.level_spacing() / CM1 << ",\"lowest_mode_cm\":"
      << md.omega[3] / CM1 << ",\"highest_mode_cm\":" << md.omega[m - 1] / CM1 << ",\"temperatures\":[";
    std::string table;
    char buf[160];
    size_t i300 = 0;
    for (size_t t = 0; t < nt; ++t) if (std::abs(temps[t] - 300.0) < std::abs(temps[i300] - 300.0)) i300 = t;
    Vec g300, d300;
    for (size_t t = 0; t < nt; ++t) {
        Vec g_all = smooth_linewidths(md, chosen, Vec(gam.begin() + t * nc, gam.begin() + (t + 1) * nc)), diff;
        double kw = wigner_conductivity(md, temps[t], &g_all, -1.0, &diff), kaf = wigner_conductivity(md, temps[t]);
        double mean = 0;
        for (int s = 3; s < m; ++s) mean += g_all[s];
        mean /= (m - 3) * CM1;
        o << (t ? "," : "") << "{\"T\":" << temps[t] << ",\"k_wigner\":" << kw << ",\"k_allen_feldman\":" << kaf
          << ",\"mean_linewidth_cm\":" << mean << "}";
        std::snprintf(buf, sizeof buf, "  %7.0f%9.3f%10.3f%8.2f cm-1\n", temps[t], kw, kaf, mean);
        table += buf;
        if (log) log(std::string("T ") + fmt(temps[t], 4) + " K  Wigner " + fmt(kw, 4) + "  harmonic " + fmt(kaf, 4) + " W/m/K");
        if (t == i300) { g300 = g_all; d300 = diff; }
    }
    o << "],\"linewidth_samples\":{\"frequency_cm\":[";
    for (size_t c = 0; c < nc; ++c) o << (c ? "," : "") << md.omega[chosen[c]] / CM1;
    o << "],\"T\":[";
    for (size_t t = 0; t < nt; ++t) o << (t ? "," : "") << temps[t];
    o << "],\"gamma_cm\":[";
    for (size_t t = 0; t < nt; ++t) {
        o << (t ? "," : "") << "[";
        for (size_t c = 0; c < nc; ++c) o << (c ? "," : "") << gam[t * nc + c] / CM1;
        o << "]";
    }
    o << "]},\"diffusivity_300\":{\"frequency_cm\":[";
    Vec centres(65);
    for (int c = 0; c < 65; ++c) centres[c] = 20.0 + (1300.0 - 20.0) * c / 64.0;
    for (int c = 0; c < 65; ++c) o << (c ? "," : "") << centres[c];
    o << "],\"mm2_per_s\":[";
    for (int c = 0; c < 65; ++c) {
        double sum = 0;
        int cnt = 0;
        for (int s = 3; s < m; ++s)
            if (std::abs(md.omega[s] / CM1 - centres[c]) < 10.0) { sum += d300[s]; ++cnt; }
        o << (c ? "," : "");
        if (cnt) o << sum / cnt * 1e6; else o << "null";
    }
    o << "]},\"vdos\":{\"frequency_cm\":[";
    const int bins = 130;
    const double top = 1400.0, width = top / bins;
    std::vector<int> hist(bins, 0);
    for (int s = 3; s < m; ++s) { int b = int(md.omega[s] / CM1 / width); if (b >= 0 && b < bins) ++hist[b]; }
    for (int b = 0; b < bins; ++b) o << (b ? "," : "") << (b + 0.5) * width;
    o << "],\"density\":[";
    for (int b = 0; b < bins; ++b) o << (b ? "," : "") << hist[b] / (double(m - 3) * width);
    o << "]},\"eta_scan_300\":[";
    int first = 1;
    for (double e : {0.5, 1.0, 2.0, 4.0, 8.0}) {
        o << (first ? "" : ",") << "{\"eta_cm\":" << e << ",\"k\":" << wigner_conductivity(md, temps[i300], &g300, e) << "}";
        first = 0;
    }
    o << "]}";
    std::ofstream(data_dir() + "/glass/" + name + ".wigner.json") << o.str();
    return table;
}

std::string wigner_report(const std::string& name) {
    const std::string base = data_dir() + "/glass/" + name;
    if (!std::ifstream(base + ".wigner.json").good())
        throw std::invalid_argument("no stored result for '" + name + "'; run: stacktherm wigner --recompute --model " + name);
    YAML::Node res = YAML::LoadFile(base + ".wigner.json"), st = res["structure"];
    std::map<double, double> kaldo;
    if (std::ifstream(base + ".kaldo.json").good())
        for (const auto& r : YAML::LoadFile(base + ".kaldo.json")["temperatures"]) kaldo[r["T"].as<double>()] = r["k_wigner"].as<double>();
    char buf[400];
    std::string out = "Conductivity of amorphous SiO2 from atoms, model '" + name + "'\n\n";
    if (st["atoms"]) {
        std::snprintf(buf, sizeof buf, "  structure: %d atoms, %.2f g/cm3, quenched at %.1e K/s; Si four-fold %.1f %%, "
                                       "O two-fold %.1f %%, Si-O bond %.3f A\n",
                      st["atoms"].as<int>(), st["density_g_cm3"].as<double>(), st["quench_rate_K_per_s"].as<double>(),
                      st["si_four_fold"].as<double>() * 100, st["o_two_fold"].as<double>() * 100, st["si_o_bond"].as<double>());
        out += buf;
    }
    std::snprintf(buf, sizeof buf, "  spectrum: %d modes from %.0f to %.0f cm-1, mean level spacing %.2f cm-1\n\n"
                                   "  %7s%9s%9s%10s%11s\n", res["modes"].as<int>(), res["lowest_mode_cm"].as<double>(),
                  res["highest_mode_cm"].as<double>(), res["level_spacing_cm"].as<double>(), "T [K]", "native", "kALDo",
                  "harmonic", "linewidth");
    out += buf;
    for (const auto& r : res["temperatures"]) {
        double T = r["T"].as<double>();
        char k2[16];
        if (kaldo.count(T)) std::snprintf(k2, sizeof k2, "%9.3f", kaldo[T]); else std::snprintf(k2, sizeof k2, "%9s", "-");
        std::snprintf(buf, sizeof buf, "  %7.0f%9.3f%s%10.3f%8.2f cm-1\n", T, r["k_wigner"].as<double>(), k2,
                      r["k_allen_feldman"].as<double>(), r["mean_linewidth_cm"].as<double>());
        out += buf;
    }
    out += "\n  'native' is the Wigner expression in this code; 'kALDo' is the quasi-harmonic\n"
           "  Green-Kubo method on the same force constants (Python bridge); 'harmonic' is the\n"
           "  Allen-Feldman limit (no anharmonic linewidths).\n\n  regularisation check at 300 K (Gaussian width in cm-1 -> W/m/K):";
    for (const auto& e : res["eta_scan_300"]) {
        std::snprintf(buf, sizeof buf, " %g -> %.3f,", e["eta_cm"].as<double>(), e["k"].as<double>());
        out += buf;
    }
    out.back() = '\n';
    return out;
}

}  // namespace st
