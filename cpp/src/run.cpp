// High-level runs (simulate, error budget, mesh convergence, temperature sweep)
// and writers: layer table for package-scale tools (CSV, JSON, Ansys APDL
// material cards), the 3D fields (VTK rectilinear grid) and the HTML report.
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace st {

namespace {

double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

std::string q(const std::string& s) { return "\"" + json_escape(s) + "\""; }
std::string num(double v) { return fmt(v, 12); }

std::string base64(const unsigned char* data, size_t n) {
    static const char* tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve((n + 2) / 3 * 4);
    for (size_t i = 0; i < n; i += 3) {
        unsigned v = data[i] << 16;
        if (i + 1 < n) v |= data[i + 1] << 8;
        if (i + 2 < n) v |= data[i + 2];
        out += tbl[(v >> 18) & 63];
        out += tbl[(v >> 12) & 63];
        out += i + 1 < n ? tbl[(v >> 6) & 63] : '=';
        out += i + 2 < n ? tbl[v & 63] : '=';
    }
    return out;
}

// the top-level model followed by each distinct sub-stack model, with a name prefix
std::vector<std::pair<std::string, const Model*>> scoped(const Model& m) {
    std::vector<std::pair<std::string, const Model*>> out = {{"", &m}};
    std::set<const Model*> seen;
    for (const auto& s : m.sub_models)
        if (seen.insert(s.second.get()).second) out.emplace_back(s.first + "/", s.second.get());
    return out;
}

std::string layer_json(const LayerProps& p) {
    std::ostringstream o;
    o << "{\"layer\":" << q(p.name) << ",\"thickness_m\":" << num(p.thickness) << ",\"kx_W_mK\":" << num(p.kx)
      << ",\"ky_W_mK\":" << num(p.ky) << ",\"kz_W_mK\":" << num(p.kz) << ",\"R_m2K_W\":" << num(p.R)
      << ",\"R_interface_above_m2K_W\":" << num(p.R_above) << ",\"kz_folded_W_mK\":"
      << num(p.thickness / (p.R + p.R_above)) << ",\"kz_rule_of_mixtures_W_mK\":" << num(p.kz_mix)
      << ",\"rhoc_J_m3K\":" << num(p.rhoc) << "}";
    return o.str();
}

std::string layers_json(const Characterization& c) {
    std::string out = "[";
    for (size_t i = 0; i < c.layers.size(); ++i) out += (i ? "," : "") + layer_json(c.layers[i]);
    return out + "]";
}

}  // namespace

// ---------------------------------------------------------------- runs
Result simulate(const Stack& stack, bool operating, double tol) {
    double t0 = now();
    Result r;
    r.stack = stack;
    r.model = stack.build();
    const Grid& g = r.model.grid;
    r.ch = characterize(g, tol);
    bool load = stack.bottom.kind == BC::Flux || stack.top.kind == BC::Flux;
    for (double v : g.q) if (v != 0.0) { load = true; break; }
    if (operating && load) r.op = operate(g, stack.bottom, stack.top, tol);
    r.wall_s = now() - t0;
    return r;
}

std::optional<double> Result::junction_resistance() const {
    if (!op || op->power <= 0) return std::nullopt;
    double sink = kInf;
    for (const BC* bc : {&stack.bottom, &stack.top})
        if (bc->fixes_temperature()) sink = std::min(sink, bc->T);
    const Grid& g = model.grid;
    return (op->T_max - sink) / (op->power / (g.lx() * g.ly()));
}

std::string Result::summary() const {
    const Grid& g = model.grid;
    char buf[400];
    std::string out;
    std::snprintf(buf, sizeof buf, "%s: %d x %d x %d cells (%zu unknowns), %.2f s\nmaterials evaluated at %.1f K\n\n",
                  stack.name.c_str(), g.nx(), g.ny(), g.nz(), g.n(), wall_s, stack.temperature);
    out += buf + ch.table() + "\n";
    if (!model.sub.empty()) {
        out += "\nsub-stacks characterised in their own unit cell:\n";
        for (const auto& s : model.sub) {
            std::snprintf(buf, sizeof buf, "  %s: kx %.3g  ky %.3g  kz %.3g W/m/K, R %.4g mm2K/W\n", s.first.c_str(),
                          s.second.kx, s.second.ky, s.second.kz, s.second.R_total * 1e6);
            out += buf;
        }
    }
    if (op) {
        double flux = op->power / (g.lx() * g.ly());
        std::snprintf(buf, sizeof buf,
                      "\noperating point: %.4g W/cm2 (%.4g uW per cell), energy imbalance %.1e\n"
                      "  T_max %.3f K at x=%.0f nm, y=%.0f nm, z=%.0f nm\n"
                      "  face temperatures: bottom %.3f K, top %.3f K; rise inside the stack %.3f K\n",
                      flux / 1e4, op->power * 1e6, op->imbalance, op->T_max, op->hotspot[0] * 1e9,
                      op->hotspot[1] * 1e9, op->hotspot[2] * 1e9, op->T_face_bottom, op->T_face_top,
                      op->T_max - std::min(op->T_face_bottom, op->T_face_top));
        out += buf;
        if (auto rj = junction_resistance()) {
            std::snprintf(buf, sizeof buf, "  junction-to-coolant resistance %.4g mm2K/W\n", *rj * 1e6);
            out += buf;
        }
    }
    for (const std::string& w : model.warnings) out += "\nwarning: " + w + "\n";
    std::set<std::string> assumed;
    for (const auto& sm : scoped(model))
        for (const auto& i : sm.second->interfaces)
            if (i.kind.compare(0, 7, "assumed") == 0) assumed.insert(i.a + "|" + i.b);
    if (!assumed.empty())
        out += "\nnote: " + std::to_string(assumed.size()) +
               " interface type(s) use the library default resistance; see the interfaces table in the result file\n";
    return out;
}

Budget error_budget(const Stack& stack) {
    // input classes: key, what it covers, and the relative uncertainty assumed for it
    static const struct { const char* key; const char* label; double unc; } classes[] = {
        {"estimated_k", "conductivities marked as estimates in the library", 0.30},
        {"tabulated_k", "tabulated conductivities with a cited source", 0.10},
        {"computed_k", "glass conductivity computed from atoms", 0.10},
        {"metal_k", "size-dependent metal model", 0.15},
        {"silicon_k", "silicon phonon model", 0.10},
        {"assumed_tbr", "interfaces on the library default", 0.50},
        {"measured_tbr", "interfaces with a measured value", 0.20},
        {"user_tbr", "bond interfaces and barriers given in the stack file", 0.30}};
    MaterialLibrary& lib = *stack.lib;
    const auto saved = lib.scale;
    auto resistance = [&] { Model m = stack.build(); return through_plane(m.grid).r_total; };
    Budget b;
    try {
        b.R_total = resistance();
        for (const auto& c : classes) {
            lib.scale = saved;
            lib.scale[c.key] = (saved.count(c.key) ? saved.at(c.key) : 1.0) * 1.1;
            double e = (resistance() / b.R_total - 1.0) / 0.1;
            if (std::abs(e) > 1e-6) b.rows.push_back({c.key, c.label, e, c.unc, std::abs(e) * c.unc});
        }
    } catch (...) { lib.scale = saved; throw; }
    lib.scale = saved;
    std::sort(b.rows.begin(), b.rows.end(), [](const BudgetRow& x, const BudgetRow& y) { return x.effect > y.effect; });
    double s = 0;
    for (const auto& r : b.rows) s += r.effect * r.effect;
    b.combined = std::sqrt(s);
    return b;
}

std::string Budget::text() const {
    char buf[200];
    std::string out = "what the total resistance rests on (ranked by effect):\n";
    std::snprintf(buf, sizeof buf, "  %-54s%11s%12s%8s\n", "input class", "elasticity", "assumed +/-", "effect");
    out += buf;
    for (const auto& r : rows) {
        std::snprintf(buf, sizeof buf, "  %-54s%+11.2f%10.0f %%%6.1f %%\n", r.label.c_str(), r.elasticity,
                      r.assumed_uncertainty * 100, r.effect * 100);
        out += buf;
    }
    std::snprintf(buf, sizeof buf, "  combined (independent classes): +/- %.0f %% on %.4g mm2K/W\n", combined * 100,
                  R_total * 1e6);
    return out + buf;
}

std::vector<ConvergenceRun> convergence(const Stack& stack, const std::vector<double>& factors) {
    std::vector<ConvergenceRun> runs;
    std::vector<std::vector<double>> props;
    for (double f : factors) {
        Stack s = stack;
        s.resolution = stack.resolution.refined(f);
        Result r = simulate(s, false);
        runs.push_back({f, r.model.grid.n(), r.ch.R_total, 0, 0, r.wall_s});
        std::vector<double> p;
        for (const auto& l : r.ch.layers) { p.push_back(l.kx); p.push_back(l.ky); p.push_back(l.kz); }
        props.push_back(p);
    }
    for (size_t i = 0; i < runs.size(); ++i) {
        for (size_t j = 0; j < props[i].size(); ++j)
            runs[i].change = std::max(runs[i].change, std::abs(props[i][j] / props.back()[j] - 1.0));
        runs[i].change_R = std::abs(runs[i].R_total / runs.back().R_total - 1.0);
    }
    return runs;
}

std::map<double, Characterization> temperature_sweep(const Stack& stack, const Vec& temperatures) {
    std::map<double, Characterization> out;
    for (double T : temperatures) {
        Stack s = stack;
        s.temperature = T;
        out[T] = simulate(s, false).ch;
    }
    return out;
}

// ---------------------------------------------------------------- writers
std::string layers_csv(const Characterization& ch) {
    std::ostringstream o;
    o << "layer,thickness_m,kx_W_mK,ky_W_mK,kz_W_mK,R_m2K_W,R_interface_above_m2K_W,kz_folded_W_mK,"
         "kz_rule_of_mixtures_W_mK,rhoc_J_m3K\n";
    for (const auto& p : ch.layers)
        o << p.name << ',' << num(p.thickness) << ',' << num(p.kx) << ',' << num(p.ky) << ',' << num(p.kz) << ','
          << num(p.R) << ',' << num(p.R_above) << ',' << num(p.thickness / (p.R + p.R_above)) << ','
          << num(p.kz_mix) << ',' << num(p.rhoc) << '\n';
    return o.str();
}

// Ansys Mechanical APDL orthotropic material cards, one material per layer.
// KZZ is the folded value (layer plus the interface above it).  Density is
// written as 1 and specific heat as the volumetric heat capacity, so DENS * C is
// correct for transient runs.  With a sweep the conductivities become tables.
std::string layers_apdl(const Characterization& ch, const std::string& name,
                        const std::map<double, Characterization>* sweep) {
    char buf[160];
    std::ostringstream o;
    o << "! stacktherm effective layer properties for '" << name << "'\n"
      << "! units: W/m/K, J/m^3/K (DENS = 1, C = rho*c), K\n"
      << "! KZZ includes the thermal resistance of the interface above each layer\n\n";
    Vec temps;
    if (sweep) for (const auto& kv : *sweep) temps.push_back(kv.first);
    for (size_t n = 0; n < ch.layers.size(); ++n) {
        const LayerProps& p = ch.layers[n];
        const int id = int(n) + 1;
        std::snprintf(buf, sizeof buf, "! layer %zu: %s  thickness %.6e m\n", n, p.name.c_str(), p.thickness);
        o << buf;
        auto folded = [](const LayerProps& l) { return l.thickness / (l.R + l.R_above); };
        if (!temps.empty()) {
            o << "MPTEMP\n";
            for (size_t i = 0; i < temps.size(); i += 6) {
                o << "MPTEMP," << i + 1;
                for (size_t j = i; j < std::min(i + 6, temps.size()); ++j) { std::snprintf(buf, sizeof buf, ",%.2f", temps[j]); o << buf; }
                o << '\n';
            }
            for (int which = 0; which < 3; ++which)
                for (size_t i = 0; i < temps.size(); i += 6) {
                    o << "MPDATA," << (which == 0 ? "KXX" : which == 1 ? "KYY" : "KZZ") << ',' << id << ',' << i + 1;
                    for (size_t j = i; j < std::min(i + 6, temps.size()); ++j) {
                        const LayerProps& l = sweep->at(temps[j]).layers[n];
                        std::snprintf(buf, sizeof buf, ",%.6e", which == 0 ? l.kx : which == 1 ? l.ky : folded(l));
                        o << buf;
                    }
                    o << '\n';
                }
        } else {
            std::snprintf(buf, sizeof buf, "MP,KXX,%d,%.6e\nMP,KYY,%d,%.6e\nMP,KZZ,%d,%.6e\n", id, p.kx, id, p.ky, id,
                          folded(p));
            o << buf;
        }
        if (p.rhoc > 0) {
            std::snprintf(buf, sizeof buf, "MP,DENS,%d,1.0\nMP,C,%d,%.6e\n", id, id, p.rhoc);
            o << buf;
        }
        o << '\n';
    }
    return o.str();
}

std::string result_json(const Result& r, bool with_field) {
    const Grid& g = r.model.grid;
    const Characterization& c = r.ch;
    std::ostringstream o;
    o << "{\"name\":" << q(r.stack.name) << ",\"temperature_K\":" << num(r.stack.temperature) << ",\"cell_m\":["
      << num(g.lx()) << ',' << num(g.ly()) << "],\"thickness_m\":" << num(g.lz()) << ",\"grid\":{\"nx\":" << g.nx()
      << ",\"ny\":" << g.ny() << ",\"nz\":" << g.nz() << ",\"unknowns\":" << g.n() << "},\"wall_s\":"
      << num(r.wall_s) << ",\"effective\":{\"kx_W_mK\":" << num(c.kx) << ",\"ky_W_mK\":" << num(c.ky)
      << ",\"kz_W_mK\":" << num(c.kz) << ",\"R_total_m2K_W\":" << num(c.R_total) << "},\"layers\":"
      << layers_json(c) << ",\"materials\":[";
    bool first = true;
    for (const auto& sm : scoped(r.model))
        for (const Instance& i : sm.second->instances) {
            o << (first ? "" : ",") << "{\"name\":" << q(sm.first + i.name) << ",\"layer\":" << q(sm.first + i.layer)
              << ",\"material\":" << q(i.material) << ",\"kx_W_mK\":" << num(i.kx) << ",\"ky_W_mK\":" << num(i.ky)
              << ",\"kz_W_mK\":" << num(i.kz) << ",\"area_fraction\":" << num(i.fraction) << ",\"model\":"
              << q(i.note) << "}";
            first = false;
        }
    o << "],\"interfaces\":[";
    first = true;
    for (const auto& sm : scoped(r.model))
        for (const InterfaceUse& i : sm.second->interfaces) {
            o << (first ? "" : ",") << "{\"where\":" << q(sm.first + i.where) << ",\"a\":" << q(sm.first + i.a)
              << ",\"b\":" << q(sm.first + i.b) << ",\"orientation\":" << q(i.orientation) << ",\"R_m2K_W\":"
              << num(i.R) << ",\"provenance\":" << q(i.kind) << ",\"source\":" << q(i.source) << "}";
            first = false;
        }
    o << "],\"substacks\":{";
    for (size_t n = 0; n < r.model.sub.size(); ++n) {
        const Characterization& s = r.model.sub[n].second;
        o << (n ? "," : "") << q(r.model.sub[n].first) << ":{\"kx_W_mK\":" << num(s.kx) << ",\"ky_W_mK\":"
          << num(s.ky) << ",\"kz_W_mK\":" << num(s.kz) << ",\"R_total_m2K_W\":" << num(s.R_total) << ",\"layers\":"
          << layers_json(s) << "}";
    }
    o << "},\"solver\":{";
    first = true;
    for (const auto& kv : c.solves) {
        o << (first ? "" : ",") << q(kv.first) << ":{\"method\":\"amg-cg\",\"iterations\":" << kv.second.iterations
          << ",\"residual\":" << num(kv.second.residual) << ",\"setup_s\":" << num(kv.second.setup_s)
          << ",\"solve_s\":" << num(kv.second.solve_s) << "}";
        first = false;
    }
    o << "},\"warnings\":[";
    for (size_t n = 0; n < r.model.warnings.size(); ++n) o << (n ? "," : "") << q(r.model.warnings[n]);
    o << "]";
    if (r.op) {
        const Operating& op = *r.op;
        auto rj = r.junction_resistance();
        o << ",\"operating\":{\"power_W\":" << num(op.power) << ",\"flux_W_m2\":" << num(op.power / (g.lx() * g.ly()))
          << ",\"q_top_W\":" << num(op.q_top) << ",\"q_bottom_W\":" << num(op.q_bottom) << ",\"energy_imbalance\":"
          << num(op.imbalance) << ",\"T_max_K\":" << num(op.T_max) << ",\"hotspot_m\":[" << num(op.hotspot[0]) << ','
          << num(op.hotspot[1]) << ',' << num(op.hotspot[2]) << "],\"junction_resistance_m2K_W\":"
          << (rj ? num(*rj) : "null") << ",\"layers\":[";
        for (size_t n = 0; n < op.layer_T_mean.size(); ++n)
            o << (n ? "," : "") << "{\"layer\":" << q(g.layer_names[n]) << ",\"T_mean_K\":" << num(op.layer_T_mean[n])
              << ",\"T_max_K\":" << num(op.layer_T_max[n]) << ",\"T_min_K\":" << num(op.layer_T_min[n]) << "}";
        o << "]}";
    }
    if (with_field && r.op) {   // temperature field, strided and quantised to 16 bits
        const Vec& T = r.op->T;
        auto stride = [](int n, int m) { return std::max(1, (n + m - 1) / m); };
        const int sx = stride(g.nx(), 48), sy = stride(g.ny(), 48), sz = stride(g.nz(), 240);
        double lo = *std::min_element(T.begin(), T.end()), hi = *std::max_element(T.begin(), T.end());
        double span = hi > lo ? hi - lo : 1.0;
        std::vector<unsigned char> bytes;
        int mx = 0, my = 0, mz = 0;
        for (int k = 0; k < g.nz(); k += sz, ++mz) {
            my = 0;
            for (int j = 0; j < g.ny(); j += sy, ++my) {
                mx = 0;
                for (int i = 0; i < g.nx(); i += sx, ++mx) {
                    unsigned v = unsigned(std::lround((T[g.idx(i, j, k)] - lo) / span * 65535.0));
                    bytes.push_back(static_cast<unsigned char>(v & 255));
                    bytes.push_back(static_cast<unsigned char>(v >> 8));
                }
            }
        }
        o << ",\"field\":{\"shape\":[" << mz << ',' << my << ',' << mx << "],\"tmin\":" << num(lo) << ",\"tmax\":"
          << num(hi) << ",\"data\":\"" << base64(bytes.data(), bytes.size()) << "\",\"x_nm\":[";
        for (int i = 0, n = 0; i < g.nx(); i += sx, ++n) o << (n ? "," : "") << fmt(0.5e9 * (g.xe[i] + g.xe[i + 1]), 8);
        o << "],\"y_nm\":[";
        for (int j = 0, n = 0; j < g.ny(); j += sy, ++n) o << (n ? "," : "") << fmt(0.5e9 * (g.ye[j] + g.ye[j + 1]), 8);
        o << "],\"z_nm\":[";
        for (int k = 0, n = 0; k < g.nz(); k += sz, ++n) o << (n ? "," : "") << fmt(0.5e9 * (g.ze[k] + g.ze[k + 1]), 8);
        o << "],\"layer\":[";
        for (int k = 0, n = 0; k < g.nz(); k += sz, ++n) o << (n ? "," : "") << g.layer_of_k[k];
        o << "],\"layer_names\":[";
        for (size_t n = 0; n < g.layer_names.size(); ++n) o << (n ? "," : "") << q(g.layer_names[n]);
        o << "]}";
    } else if (with_field) {
        o << ",\"field\":null";
    }
    o << "}";
    return o.str();
}

std::string report_html(const Result& r) {
    std::ifstream f(data_dir() + "/report_template.html");
    if (!f) throw std::runtime_error("report template not found in " + data_dir());
    std::stringstream ss;
    ss << f.rdbuf();
    std::string html = ss.str(), payload = result_json(r, true);
    for (size_t p = 0; (p = payload.find("</", p)) != std::string::npos; p += 3) payload.replace(p, 2, "<\\/");
    auto replace = [&](const std::string& key, const std::string& value) {
        size_t p = html.find(key);
        if (p != std::string::npos) html.replace(p, key.size(), value);
    };
    replace("__TITLE__", r.stack.name);
    replace("__DATA__", payload);
    return html;
}

// VTK XML rectilinear grid with cell data (binary).  Coordinates are written in
// micrometres: nanometre values in metres fall below the tolerances of viewers.
void write_vtr(const std::string& path, const Result& r) {
    const Grid& g = r.model.grid;
    std::vector<std::string> blocks;
    std::ostringstream cells, coords;
    uint64_t offset = 0;
    auto add = [&](const void* data, size_t bytes) {
        std::string b(8 + bytes, '\0');
        uint64_t n = bytes;
        std::memcpy(&b[0], &n, 8);
        std::memcpy(&b[8], data, bytes);
        blocks.push_back(std::move(b));
        uint64_t start = offset;
        offset += 8 + bytes;
        return start;
    };
    auto floats = [&](const char* name, const Vec& v) {
        std::vector<float> f(v.begin(), v.end());
        cells << "<DataArray type=\"Float32\" Name=\"" << name << "\" format=\"appended\" offset=\""
              << add(f.data(), f.size() * 4) << "\"/>\n";
    };
    auto ints = [&](const char* name, const std::vector<int32_t>& v) {
        cells << "<DataArray type=\"Int32\" Name=\"" << name << "\" format=\"appended\" offset=\""
              << add(v.data(), v.size() * 4) << "\"/>\n";
    };
    if (r.op) floats("T_K", r.op->T);
    floats("kx", g.kx); floats("ky", g.ky); floats("kz", g.kz); floats("q_W_m3", g.q);
    ints("material", std::vector<int32_t>(g.mat.begin(), g.mat.end()));
    std::vector<int32_t> layer(g.n());
    const size_t nxy = size_t(g.nx()) * g.ny();
    for (int k = 0; k < g.nz(); ++k) std::fill(layer.begin() + k * nxy, layer.begin() + (k + 1) * nxy, g.layer_of_k[k]);
    ints("layer", layer);
    for (auto [name, e] : {std::pair<const char*, const Vec*>{"x_um", &g.xe}, {"y_um", &g.ye}, {"z_um", &g.ze}}) {
        Vec um(e->size());
        for (size_t i = 0; i < um.size(); ++i) um[i] = (*e)[i] * 1e6;
        coords << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"appended\" offset=\""
               << add(um.data(), um.size() * 8) << "\"/>\n";
    }
    std::ofstream f(path, std::ios::binary);
    f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
         "header_type=\"UInt64\">\n<RectilinearGrid WholeExtent=\"0 " << g.nx() << " 0 " << g.ny() << " 0 " << g.nz()
      << "\">\n<Piece Extent=\"0 " << g.nx() << " 0 " << g.ny() << " 0 " << g.nz() << "\">\n<CellData Scalars=\""
      << (r.op ? "T_K" : "kx") << "\">\n" << cells.str() << "</CellData>\n<Coordinates>\n" << coords.str()
      << "</Coordinates>\n</Piece>\n</RectilinearGrid>\n<AppendedData encoding=\"raw\">\n_";
    for (const std::string& b : blocks) f.write(b.data(), std::streamsize(b.size()));
    f << "\n</AppendedData>\n</VTKFile>\n";
}

std::vector<std::string> write_all(const Result& r, const std::string& dir, bool report,
                                   const std::map<double, Characterization>* sweep) {
    std::filesystem::create_directories(dir);
    std::vector<std::string> files;
    auto text = [&](const std::string& name, const std::string& body) {
        std::ofstream(dir + "/" + name) << body;
        files.push_back(dir + "/" + name);
    };
    text("result.json", result_json(r));
    text("layers.csv", layers_csv(r.ch));
    text("layers_apdl.mac", layers_apdl(r.ch, r.stack.name, sweep));
    if (sweep) {
        std::ostringstream o;
        o << "temperature_K,layer,thickness_m,kx_W_mK,ky_W_mK,kz_W_mK,kz_folded_W_mK,R_m2K_W,R_interface_above_m2K_W\n";
        for (const auto& kv : *sweep)
            for (const auto& p : kv.second.layers)
                o << num(kv.first) << ',' << p.name << ',' << num(p.thickness) << ',' << num(p.kx) << ',' << num(p.ky)
                  << ',' << num(p.kz) << ',' << num(p.thickness / (p.R + p.R_above)) << ',' << num(p.R) << ','
                  << num(p.R_above) << '\n';
        text("layers_vs_temperature.csv", o.str());
    }
    write_vtr(dir + "/field.vtr", r);
    files.push_back(dir + "/field.vtr");
    if (report) text("report.html", report_html(r));
    return files;
}

}  // namespace st
