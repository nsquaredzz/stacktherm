// Material models and the library that resolves a material name plus its
// geometric context into conductivities, and a pair of materials into an
// interface resistance.
//
// Metals: Fuchs-Sondheimer surface scattering and Mayadas-Shatzkes grain-boundary
//   scattering give the resistivity of a wire or film; Wiedemann-Franz converts
//   it, with a Lorenz number pinned to the bulk metal.
//   (Fuchs 1938; Sondheimer 1952; Mayadas & Shatzkes 1970; Steinhoegl et al.
//   2002, 2005; rho*lambda products from Gall, J. Appl. Phys. 119, 085101 (2016))
// Silicon: kinetic-theory integral over a Born-von Karman acoustic dispersion
//   with isotope and Umklapp scattering, boundary suppression applied mode by
//   mode.  (Callaway 1959; Dames & Chen 2004; Yang & Dames 2013)
// Interfaces: measured boundary conductance where the library has one, a stated
//   default otherwise; the diffuse mismatch model is available as a fallback.
//   (Swartz & Pohl, Rev. Mod. Phys. 61, 605 (1989))
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace st {

namespace {
const double HBAR = 1.054571817e-34, KB = 1.380649e-23;

struct Nodes { Vec x, w; };
const Nodes& unit_nodes() {                       // 96 Gauss-Legendre nodes on (0, 1)
    static const Nodes n = [] { Nodes v; gauss_legendre(96, 0.0, 1.0, v.x, v.w); return v; }();
    return n;
}
double cross_plane(double kappa) { return 1.0 / (1.0 + 4.0 / (3.0 * kappa)); }
}  // namespace

// ---------------------------------------------------------------- metals
double fuchs_sondheimer(double kappa, double p) {
    const Nodes& n = unit_nodes();
    double integral = 0.0;
    for (size_t i = 0; i < n.x.size(); ++i) {
        double u = n.x[i], e = std::exp(-kappa / u);
        integral += (u - u * u * u) * (1.0 - e) / (1.0 - p * e) * n.w[i];
    }
    return 1.0 - 1.5 * (1.0 - p) / kappa * integral;
}

double mayadas_shatzkes(double a) {
    if (a < 1e-4) return 1.0 + 1.5 * a;
    return 1.0 / (3.0 * (1.0 / 3.0 - a / 2.0 + a * a - a * a * a * std::log1p(1.0 / a)));
}

double Metal::rho_line(double w, double h, double T) const {
    double lam = mfp(T), alpha = lam / (grain_factor * std::min(w, h)) * R / (1.0 - R);
    double fs = 0.375 * 1.2 * (1.0 - p) * lam * (1.0 / w + 1.0 / h);
    return rho_bulk(T) * (mayadas_shatzkes(alpha) + fs);
}
double Metal::rho_film(double d, double T) const {
    double lam = mfp(T), alpha = lam / (grain_factor * d) * R / (1.0 - R);
    return rho_bulk(T) * (mayadas_shatzkes(alpha) + 1.0 / fuchs_sondheimer(d / lam, p) - 1.0);
}
double Metal::rho_film_cross(double d, double T) const {
    double lam = mfp(T), alpha = lam / (grain_factor * d) * R / (1.0 - R);
    return rho_bulk(T) * (mayadas_shatzkes(alpha) + 1.0 / cross_plane(d / lam) - 1.0);
}
std::pair<double, double> Metal::k_film(double d, double T) const {
    return {lorenz() * T / rho_film(d, T), lorenz() * T / rho_film_cross(d, T)};
}

const std::map<std::string, Metal>& builtin_metals() {
    // bulk resistivity (293 K) and rho*lambda from Gall (2016); k, density, cp are handbook values
    static const std::map<std::string, Metal> m = {
        {"Cu", {"Cu", 1.678e-8, 3.93e-3, 6.70e-16, 401.0, 0.0, 0.3, 1.0, 8960.0, 385.0}},
        {"Al", {"Al", 2.650e-8, 4.29e-3, 5.01e-16, 237.0, 0.0, 0.3, 1.0, 2700.0, 897.0}},
        {"W", {"W", 5.28e-8, 4.5e-3, 8.20e-16, 173.0, 0.0, 0.3, 1.0, 19300.0, 132.0}},
        {"Co", {"Co", 6.2e-8, 6.0e-3, 7.31e-16, 100.0, 0.0, 0.3, 1.0, 8900.0, 421.0}},
        {"Ru", {"Ru", 7.8e-8, 4.1e-3, 5.14e-16, 117.0, 0.0, 0.3, 1.0, 12370.0, 238.0}},
        {"Mo", {"Mo", 5.34e-8, 4.6e-3, 5.99e-16, 138.0, 0.0, 0.3, 1.0, 10280.0, 251.0}}};
    return m;
}

// ---------------------------------------------------------------- silicon
double fs_phonon(double delta) {
    double d = std::max(delta, 1e-12);
    return 1.0 - 3.0 / (8.0 * d) * (1.0 - 4.0 * expint(3, d) + 4.0 * expint(5, d));
}

namespace {
// per-mode bulk contribution dk (W/m/K) and mean free path (m)
void phonon_modes(const PhononCrystal& c, double T, Vec& dk, Vec& mfp) {
    static const Nodes n = [] { Nodes v; gauss_legendre(400, 0.0, 1.0, v.x, v.w); return v; }();
    const double q0 = std::cbrt(6.0 * kPi * kPi * c.n_prim);
    dk.resize(n.x.size());
    mfp.resize(n.x.size());
    for (size_t i = 0; i < n.x.size(); ++i) {
        double q = q0 * n.x[i], wq = q0 * n.w[i], phase = 0.5 * kPi * q / q0;
        double omega = 2.0 * c.v_s * q0 / kPi * std::sin(phase), v = c.v_s * std::cos(phase);
        double y = HBAR * omega / (KB * T), em = std::expm1(y);
        double heat = 3.0 * q * q / (2.0 * kPi * kPi) * KB * y * y * std::exp(y) / (em * em);
        double tau = 1.0 / (c.A * std::pow(omega, 4) + c.P * omega * omega * T * std::exp(-c.C_U / T));
        mfp[i] = v * tau;
        dk[i] = heat * v * mfp[i] / 3.0 * wq;
    }
}
}  // namespace

double PhononCrystal::k_bulk(double T) const {
    Vec dk, mfp;
    phonon_modes(*this, T, dk, mfp);
    double s = 0;
    for (double v : dk) s += v;
    return s;
}
std::pair<double, double> PhononCrystal::k_film(double d, double T) const {
    Vec dk, mfp;
    phonon_modes(*this, T, dk, mfp);
    double in = 0, cross = 0;
    for (size_t i = 0; i < dk.size(); ++i) {
        in += dk[i] * fs_phonon(d / mfp[i]);
        cross += dk[i] / (1.0 + 4.0 * mfp[i] / (3.0 * d));
    }
    return {in, cross};
}
double PhononCrystal::k_wire(double w, double h, double T) const {
    Vec dk, mfp;
    phonon_modes(*this, T, dk, mfp);
    double s = 0, casimir = 1.12 * std::sqrt(w * h);
    for (size_t i = 0; i < dk.size(); ++i) s += dk[i] / (1.0 + mfp[i] / casimir);
    return s;
}

const PhononCrystal& silicon() {
    // v_s: Debye average of the [100] velocities.  A: isotope scattering of
    // natural Si, V0*Gamma/(4 pi v_s^3), not fitted.  P and C_U: fitted to the
    // bulk k(T) curve (within 2 % from 150 K to 600 K).
    static const PhononCrystal si{"Si", 6354.0, 2.5e28, 1.247e-45, 1.518e-19, 137.7, 2329.0, 705.0};
    return si;
}

// ---------------------------------------------------------------- interfaces
namespace {
struct Debye { double vl, vt, n; };
const std::map<std::string, Debye>& debye_table() {
    static const std::map<std::string, Debye> t = {
        {"Si", {8433.0, 5843.0, 5.00e28}}, {"SiO2", {5953.0, 3743.0, 6.62e28}},
        {"Cu", {4760.0, 2325.0, 8.47e28}}, {"Al", {6420.0, 3040.0, 6.02e28}},
        {"W", {5220.0, 2890.0, 6.31e28}}};
    return t;
}
}  // namespace

double dmm_conductance(const std::string& a, const std::string& b, double T) {
    const auto& t = debye_table();
    if (!t.count(a) || !t.count(b)) return 0.0;
    auto branches = [](const Debye& d, double v[3], double w[3]) {
        double qd = std::cbrt(6.0 * kPi * kPi * d.n);
        v[0] = d.vl; v[1] = v[2] = d.vt;
        for (int i = 0; i < 3; ++i) w[i] = v[i] * qd;
    };
    double va[3], wa[3], vb[3], wb[3];
    branches(t.at(a), va, wa);
    branches(t.at(b), vb, wb);
    double wmax = 0;
    for (int i = 0; i < 3; ++i) wmax = std::max({wmax, wa[i], wb[i]});
    const int N = 4000;
    const double dw = wmax / N;
    double sum = 0.0;
    for (int m = 1; m <= N; ++m) {
        double omega = m * dw, sa = 0, sb = 0;
        for (int i = 0; i < 3; ++i) {
            if (omega <= wa[i]) sa += 1.0 / (va[i] * va[i]);
            if (omega <= wb[i]) sb += 1.0 / (vb[i] * vb[i]);
        }
        if (sa + sb <= 0) continue;
        double x = HBAR * omega / (KB * T), em = std::expm1(x);
        sum += HBAR * omega * omega * omega * (x / T * std::exp(x) / (em * em)) * (sb / (sa + sb)) * sa;
    }
    return sum * dw / (8.0 * kPi * kPi);
}

double cahill_pohl(double v_l, double v_t, double n_atoms, double T) {
    double total = 0.0;
    for (double v : {v_l, v_t, v_t}) {
        double theta = v * HBAR / KB * std::cbrt(6.0 * kPi * kPi * n_atoms);
        const int N = 2000;
        double integral = 0.0, prev = 0.0, uprev = 0.0;
        for (int i = 0; i < N; ++i) {
            double x = 1e-6 + (1.0 - 1e-6) * i / (N - 1), u = x * theta / T, em = std::expm1(u);
            double f = u * u * u * std::exp(u) / (em * em);
            if (i) integral += 0.5 * (f + prev) * (u - uprev);
            prev = f;
            uprev = u;
        }
        total += v * (T / theta) * (T / theta) * integral;
    }
    return std::cbrt(kPi / 6.0) * KB * std::pow(n_atoms, 2.0 / 3.0) * total;
}

// ---------------------------------------------------------------- library
namespace {
std::pair<std::string, std::string> pair_key(const std::string& a, const std::string& b) {
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}
double interp(const Vec& x, const Vec& y, double v) {
    if (v <= x.front()) return y.front();
    if (v >= x.back()) return y.back();
    size_t i = size_t(std::upper_bound(x.begin(), x.end(), v) - x.begin());
    return y[i - 1] + (y[i] - y[i - 1]) * (v - x[i - 1]) / (x[i] - x[i - 1]);
}
bool starts_with(const std::string& s, const std::string& p) { return s.compare(0, p.size(), p) == 0; }
}  // namespace

MaterialLibrary::MaterialLibrary() {
    YAML::Node data = YAML::LoadFile(data_dir() + "/library.yaml");
    for (const auto& kv : data["solids"]) {
        const YAML::Node& d = kv.second;
        solids_[kv.first.as<std::string>()] = {d["k"].as<double>(), d["density"] ? d["density"].as<double>() : 0.0,
                                               d["cp"] ? d["cp"].as<double>() : 0.0,
                                               d["source"] ? d["source"].as<std::string>() : ""};
    }
    metals = builtin_metals();
    if (data["metal_scattering"])
        for (const auto& kv : data["metal_scattering"]) {
            Metal& m = metals.at(kv.first.as<std::string>());
            if (kv.second["p"]) m.p = kv.second["p"].as<double>();
            if (kv.second["R"]) m.R = kv.second["R"].as<double>();
            if (kv.second["grain_factor"]) m.grain_factor = kv.second["grain_factor"].as<double>();
        }
    if (data["interfaces"])
        for (const auto& e : data["interfaces"]) {
            auto a = e["pair"][0].as<std::string>(), b = e["pair"][1].as<std::string>();
            measured_[pair_key(a, b)] = {1.0 / parse_quantity(e["tbc"].as<std::string>(), "MW/m2K"),
                                         e["source"] ? e["source"].as<std::string>() : ""};
        }
    default_tbr = 1.0 / parse_quantity(data["default_interface"]["tbc"].as<std::string>(), "MW/m2K");
    if (data["default_interface"]["fallback"]) fallback = data["default_interface"]["fallback"].as<std::string>();
    if (data["wigner"])
        for (const auto& kv : data["wigner"]) {
            const YAML::Node& d = kv.second;
            std::vector<std::string> models;
            if (d["models"]) for (const auto& m : d["models"]) models.push_back(m.as<std::string>());
            else models.push_back(d["model"].as<std::string>());
            add_table(kv.first.as<std::string>(), models, d["base"] ? d["base"].as<std::string>() : "",
                      d["density"] ? d["density"].as<double>() : 0.0, d["cp"] ? d["cp"].as<double>() : 0.0);
        }
}

// Material whose k(T) comes from stored atomistic results.  Several models are
// averaged; for each the kALDo table is preferred, then the native Wigner table.
bool MaterialLibrary::add_table(const std::string& name, const std::vector<std::string>& models,
                                const std::string& base, double density, double cp) {
    std::vector<std::pair<Vec, Vec>> curves;
    std::string used;
    for (const auto& model : models)
        for (const char* backend : {"kaldo", "wigner"}) {
            std::string path = data_dir() + "/glass/" + model + "." + backend + ".json";
            if (!std::ifstream(path).good()) continue;
            YAML::Node rows = YAML::LoadFile(path)["temperatures"];
            Vec t, k;
            for (const auto& r : rows) { t.push_back(r["T"].as<double>()); k.push_back(r["k_wigner"].as<double>()); }
            curves.emplace_back(t, k);
            used += (used.empty() ? "" : ", ") + model + (std::string(backend) == "kaldo" ? " (kALDo)" : " (native)");
            break;
        }
    if (curves.empty()) return false;
    Vec temps;
    for (const auto& c : curves)
        for (double t : c.first) {
            bool inside = true;
            for (const auto& o : curves) inside = inside && t >= o.first.front() && t <= o.first.back();
            if (inside && std::find(temps.begin(), temps.end(), t) == temps.end()) temps.push_back(t);
        }
    std::sort(temps.begin(), temps.end());
    Table tab;
    tab.T = temps;
    for (double t : temps) {
        double s = 0;
        for (const auto& c : curves) s += interp(c.first, c.second, t);
        tab.k.push_back(s / curves.size());
    }
    tab.rhoc = density * cp;
    tab.note = "computed from atoms: " + used;
    tables_[name] = tab;
    if (!base.empty()) base_[name] = base;
    return true;
}

std::string MaterialLibrary::table_note(const std::string& name) const {
    auto it = tables_.find(name);
    return it == tables_.end() ? "" : it->second.note;
}

void MaterialLibrary::add(const std::string& name, const MaterialSpec& s) {
    if (!s.tables.empty()) {
        if (!add_table(name, s.tables, s.base, s.density.value_or(0.0), s.cp.value_or(0.0)))
            throw std::invalid_argument("no stored conductivity table for '" + s.tables.front() +
                                        "'; create one with 'stacktherm wigner'");
        return;
    }
    if (s.base.empty()) {
        if (!s.k) throw std::invalid_argument("material '" + name + "' needs 'k' or 'base'");
        solids_[name] = {*s.k, s.density.value_or(0.0), s.cp.value_or(0.0), s.source.empty() ? "user" : s.source};
        return;
    }
    base_[name] = canonical(s.base);
    if (metals.count(s.base)) {
        Metal m = metals.at(s.base);
        if (s.p) m.p = *s.p;
        if (s.R) m.R = *s.R;
        if (s.grain_factor) m.grain_factor = *s.grain_factor;
        metals[name] = m;
    } else if (s.base == "Si") {
        throw std::invalid_argument("silicon cannot be derived into '" + name + "' yet");
    } else {
        auto it = solids_.find(s.base);
        if (it == solids_.end()) throw std::invalid_argument("unknown base material '" + s.base + "'");
        double phi = s.porosity.value_or(0.0);
        // differential-effective-medium scaling for spherical pores
        double k = s.k ? *s.k : it->second.k * std::pow(1.0 - phi, 1.5);
        solids_[name] = {k, s.density.value_or(it->second.density * (1.0 - phi)), s.cp.value_or(it->second.cp),
                         "derived from " + s.base};
    }
}

void MaterialLibrary::set_interface(const std::string& a, const std::string& b, double resistance) {
    user_tbr_[pair_key(canonical(a), canonical(b))] = resistance;
}

void MaterialLibrary::without_interfaces() {
    measured_.clear();
    user_tbr_.clear();
    default_tbr = 0.0;
    fallback = "assumed";
}

std::string MaterialLibrary::canonical(const std::string& name) const {
    auto it = base_.find(name);
    return it == base_.end() ? name : it->second;
}

double MaterialLibrary::scale_of(const std::string& key) const {
    auto it = scale.find(key);
    return it == scale.end() ? 1.0 : it->second;
}

Props MaterialLibrary::props(const std::string& name, const std::string& shape_in, double d0, double d1, double T,
                             bool size_effect) const {
    const std::string shape = size_effect ? shape_in : "bulk";
    const std::string base = canonical(name);
    char buf[96];
    if (metals.count(name)) {
        const Metal& m = metals.at(name);
        const double f = scale_of("metal_k"), rhoc = m.density * m.cp;
        if (shape == "film") {
            auto [ki, kc] = m.k_film(d0, T);
            std::snprintf(buf, sizeof buf, "FS+MS film, p=%g R=%g", m.p, m.R);
            return {ki * f, ki * f, kc * f, rhoc, base, buf};
        }
        if (shape == "wire") {
            double k = m.k_line(d0, d1, T) * f;
            std::snprintf(buf, sizeof buf, "FS+MS wire, p=%g R=%g", m.p, m.R);
            return {k, k, k, rhoc, base, buf};
        }
        double k = m.k_bulk_T(T) * f;
        return {k, k, k, rhoc, base, "bulk, Wiedemann-Franz"};
    }
    if (name == "Si" || base == "Si") {
        const PhononCrystal& c = silicon();
        const double f = scale_of("silicon_k"), rhoc = c.density * c.cp;
        if (shape == "film") {
            auto [ki, kc] = c.k_film(d0, T);
            return {ki * f, ki * f, kc * f, rhoc, base, "phonon BTE film (FS in-plane, EPRT cross-plane)"};
        }
        if (shape == "wire") {
            double k = c.k_wire(d0, d1, T) * f;
            return {k, k, k, rhoc, base, "phonon BTE wire (Casimir)"};
        }
        double k = c.k_bulk(T) * f;
        return {k, k, k, rhoc, base, "phonon BTE bulk"};
    }
    auto tab = tables_.find(name);
    if (tab != tables_.end()) {
        double k = interp(tab->second.T, tab->second.k, T) * scale_of("computed_k");
        return {k, k, k, tab->second.rhoc, base, tab->second.note};
    }
    auto sol = solids_.find(name);
    if (sol != solids_.end()) {
        const Solid& s = sol->second;
        double k = s.k * scale_of(starts_with(s.source, "estimate") ? "estimated_k" : "tabulated_k");
        return {k, k, k, s.density * s.cp, base, s.source.empty() ? "library" : s.source};
    }
    throw std::invalid_argument("unknown material '" + name + "'; define it under 'materials:' in the stack file");
}

Tbr MaterialLibrary::tbr(const std::string& a_in, const std::string& b_in, double T) const {
    const std::string a = canonical(a_in), b = canonical(b_in);
    if (a == b || starts_with(a, "homog:") || starts_with(b, "homog:")) return {0.0, "same", ""};
    auto key = pair_key(a, b);
    auto u = user_tbr_.find(key);
    if (u != user_tbr_.end()) return {u->second * scale_of("user_tbr"), "user", ""};
    auto m = measured_.find(key);
    if (m != measured_.end()) return {m->second.first * scale_of("measured_tbr"), "measured", m->second.second};
    if (fallback == "dmm") {
        double g = dmm_conductance(a, b, T);
        if (g > 0) return {scale_of("assumed_tbr") / g, "dmm", "diffuse mismatch model"};
    }
    return {default_tbr * scale_of("assumed_tbr"), "assumed", "library default"};
}

}  // namespace st
