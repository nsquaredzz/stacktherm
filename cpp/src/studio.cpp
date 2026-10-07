// What the interactive studio asks of the core: one solve packaged for 3D
// viewing, the bulk-property variant of a stack, and the material-engine and
// atomistic results as JSON.
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>

#include <yaml-cpp/yaml.h>

namespace st {

namespace {
std::string q(const std::string& s) { return "\"" + json_escape(s) + "\""; }
std::string num(double v) { return fmt(v, 10); }

template <class It, class F>
std::string list(It begin, It end, F f) {
    std::string out = "[";
    for (It it = begin; it != end; ++it) out += (it == begin ? "" : ",") + f(*it);
    return out + "]";
}
std::string numbers(const Vec& v, double scale = 1.0) {
    return list(v.begin(), v.end(), [&](double x) { return num(x * scale); });
}

void strip_layers(std::vector<Layer>& layers) {
    std::vector<Layer> kept;
    for (Layer l : layers) {
        if (l.type == Layer::Interface) continue;
        l.size_effect = false;
        l.barrier = Barrier();
        strip_layers(l.layers);
        kept.push_back(l);
    }
    layers = kept;
}

std::string slurp(const std::string& path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
}  // namespace

Stack bulk_variant(const Stack& stack) {
    Stack s = stack;
    s.lib = std::make_shared<MaterialLibrary>(*stack.lib);
    s.lib->without_interfaces();
    strip_layers(s.layers);
    return s;
}

StudioRun studio_run(const Stack& stack, bool budget, bool with_files) {
    auto t0 = std::chrono::steady_clock::now();
    Result r = simulate(stack);
    const Model& model = r.model;
    const Grid& g = model.grid;
    const Characterization& c = r.ch;
    StudioRun out;

    Vec T;
    std::string note;
    double power = 0;
    std::unique_ptr<System> sys;
    if (r.op) {
        sys.reset(new System(assemble(g, stack.bottom, stack.top)));
        T = r.op->T;
        note = "operating point";
        power = r.op->power;
    } else {   // nothing to drive it: show the unit-drop solve used for the layer table
        sys.reset(new System(assemble(g, BC::dirichlet(1.0), BC::dirichlet(0.0), true, true, 0, 0, false)));
        T = LinearSolver(sys->A).solve(sys->b);
        note = "unit temperature drop (the stack file has no heat source)";
    }
    Vec qx, qy, qz, qb, qt;
    sys->cell_flux(T, qx, qy, qz);
    sys->boundary_flows(T, qb, qt);
    double q_top = 0;
    for (double v : qt) q_top += v;
    out.T.assign(T.begin(), T.end());
    out.qx.assign(qx.begin(), qx.end());
    out.qy.assign(qy.begin(), qy.end());
    out.qz.assign(qz.begin(), qz.end());
    out.kz.assign(g.kz.begin(), g.kz.end());
    out.feature.resize(g.n());
    out.mat.resize(g.n());
    for (size_t i = 0; i < g.n(); ++i) { out.feature[i] = model.instances[g.mat[i]].feature; out.mat[i] = short(g.mat[i]); }

    const size_t nxy = size_t(g.nx()) * g.ny();
    int heat_lo = -1, heat_hi = -1;
    for (int k = 0; k < g.nz(); ++k)
        for (size_t col = 0; col < nxy; ++col)
            if (g.q[size_t(k) * nxy + col] != 0.0) { if (heat_lo < 0) heat_lo = k; heat_hi = k + 1; break; }
    double sink = kInf, tmax = *std::max_element(T.begin(), T.end()), tmin = *std::min_element(T.begin(), T.end());
    for (const BC* bc : {&stack.bottom, &stack.top})
        if (bc->fixes_temperature()) sink = std::min(sink, bc->T);
    const bool has_sink = std::isfinite(sink);
    double thinnest = kInf;
    for (const auto& p : c.layers) thinnest = std::min(thinnest, p.thickness);
    const bool tall = g.lz() > 4 * std::max(g.lx(), g.ly()) || thinnest < 0.01 * g.lz();

    std::ostringstream o;
    o << "{\"scene\":\"custom\",\"title\":" << q(stack.name) << ",\"zmode\":" << q(tall ? "schematic" : "true")
      << ",\"field_note\":" << q(note) << ",\"grid\":{\"nx\":" << g.nx() << ",\"ny\":" << g.ny() << ",\"nz\":" << g.nz()
      << ",\"xe\":" << numbers(g.xe, 1e9) << ",\"ye\":" << numbers(g.ye, 1e9) << ",\"ze\":" << numbers(g.ze, 1e9)
      << ",\"layer\":" << list(g.layer_of_k.begin(), g.layer_of_k.end(), [](int v) { return std::to_string(v); })
      << ",\"layer_names\":" << list(g.layer_names.begin(), g.layer_names.end(), q) << "},\"materials\":"
      << list(model.instances.begin(), model.instances.end(), [](const Instance& i) {
             return "{\"name\":" + q(i.name) + ",\"material\":" + q(i.material) + ",\"kx\":" + num(i.kx) + ",\"ky\":" +
                    num(i.ky) + ",\"kz\":" + num(i.kz) + ",\"model\":" + q(i.note) + ",\"feature\":" +
                    (i.feature ? "true" : "false") + "}";
         });
    // tables include the materials and interfaces resolved inside sub-stacks
    std::string mats = "[", ifaces = "[";
    std::vector<std::pair<std::string, const Model*>> scoped = {{"", &model}};
    for (const auto& sm : model.sub_models) {
        bool seen = false;
        for (const auto& e : scoped) seen = seen || e.second == sm.second.get();
        if (!seen) scoped.emplace_back(sm.first + "/", sm.second.get());
    }
    for (const auto& sm : scoped) {
        for (const Instance& i : sm.second->instances) {
            if (i.material == "homogenised") continue;
            mats += std::string(mats.size() > 1 ? "," : "") + "{\"name\":" + q(sm.first + i.name) + ",\"material\":" +
                    q(i.material) + ",\"kx\":" + num(i.kx) + ",\"ky\":" + num(i.ky) + ",\"kz\":" + num(i.kz) +
                    ",\"model\":" + q(i.note) + "}";
        }
        for (const InterfaceUse& i : sm.second->interfaces)
            ifaces += std::string(ifaces.size() > 1 ? "," : "") + "{\"a\":" + q(sm.first + i.a) + ",\"b\":" +
                      q(sm.first + i.b) + ",\"face\":" + q(i.orientation) + ",\"R\":" + num(i.R) + ",\"provenance\":" +
                      q(i.kind) + ",\"source\":" + q(i.source) + "}";
    }
    o << ",\"materials_table\":" << mats << "],\"interfaces\":" << ifaces << "],\"layers\":"
      << list(c.layers.begin(), c.layers.end(), [](const LayerProps& p) {
             return "{\"layer\":" + q(p.name) + ",\"thickness_m\":" + num(p.thickness) + ",\"kx_W_mK\":" + num(p.kx) +
                    ",\"ky_W_mK\":" + num(p.ky) + ",\"kz_W_mK\":" + num(p.kz) + ",\"R_m2K_W\":" + num(p.R) +
                    ",\"R_interface_above_m2K_W\":" + num(p.R_above) + ",\"kz_rule_of_mixtures_W_mK\":" + num(p.kz_mix) + "}";
         })
      << ",\"effective\":{\"kx\":" << num(c.kx) << ",\"ky\":" << num(c.ky) << ",\"kz\":" << num(c.kz) << ",\"R\":"
      << num(c.R_total) << ",\"thickness\":" << num(g.lz()) << "},\"bulk\":null,\"operating\":{\"T_ref\":"
      << num(has_sink ? sink : 0.0) << ",\"T_max\":" << num(tmax) << ",\"rise\":"
      << num(tmax - (has_sink && power > 0 ? sink : tmin)) << ",\"power\":" << num(power) << ",\"flux\":"
      << num(power / (g.lx() * g.ly())) << ",\"up\":" << (power > 0 ? num(q_top / power) : "null") << ",\"seed\":";
    if (heat_lo >= 0) o << "[" << heat_lo << "," << heat_hi << "]"; else o << "null";
    o << "},\"warnings\":" << list(model.warnings.begin(), model.warnings.end(), q) << ",\"cells\":" << g.n()
      << ",\"budget\":";
    if (budget) {
        Budget b = error_budget(stack);
        o << "{\"R_total\":" << num(b.R_total) << ",\"combined\":" << num(b.combined) << ",\"rows\":"
          << list(b.rows.begin(), b.rows.end(), [](const BudgetRow& x) {
                 return "{\"key\":" + q(x.key) + ",\"label\":" + q(x.label) + ",\"elasticity\":" + num(x.elasticity) +
                        ",\"assumed_uncertainty\":" + num(x.assumed_uncertainty) + ",\"effect\":" + num(x.effect) + "}";
             }) << "}";
    } else {
        o << "null";
    }
    o << ",\"seconds\":" << num(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()) << "}";
    out.json = o.str();
    if (with_files) {
        out.files["result.json"] = result_json(r);
        out.files["layers.csv"] = layers_csv(c);
        out.files["layers_apdl.mac"] = layers_apdl(c, stack.name);
        out.files["report.html"] = report_html(r);
    }
    return out;
}

// conductivity against feature size for the size-dependent materials, with the
// published silicon measurements and the interface comparison
std::string material_curves_json(double T, double p, double R, double aspect) {
    const int N = 70;
    Vec widths(N), thick(N);
    for (int i = 0; i < N; ++i) {
        widths[i] = std::pow(10.0, std::log10(5e-9) + (std::log10(1e-6) - std::log10(5e-9)) * i / (N - 1));
        thick[i] = std::pow(10.0, std::log10(5e-9) + (std::log10(1e-5) - std::log10(5e-9)) * i / (N - 1));
    }
    std::ostringstream o;
    o << "{\"T\":" << num(T) << ",\"widths_nm\":" << numbers(widths, 1e9) << ",\"metals\":[";
    int first = 1;
    for (const char* name : {"Cu", "Ru", "Co", "W"}) {
        Metal m = builtin_metals().at(name);
        m.p = p;
        m.R = R;
        Vec k(N);
        for (int i = 0; i < N; ++i) k[i] = m.k_line(widths[i], aspect * widths[i], T);
        o << (first ? "" : ",") << "{\"name\":" << q(name) << ",\"bulk\":" << num(m.k_bulk_T(T)) << ",\"mfp_nm\":"
          << num(m.mfp(T) * 1e9) << ",\"k\":" << numbers(k) << "}";
        first = 0;
    }
    Vec in(N), cross(N), fin(N);
    for (int i = 0; i < N; ++i) {
        auto f = silicon().k_film(thick[i], T);
        in[i] = f.first;
        cross[i] = f.second;
        fin[i] = silicon().k_wire(thick[i], thick[i], T);
    }
    YAML::Node data = YAML::LoadFile(data_dir() + "/published.yaml");
    o << "],\"thickness_nm\":" << numbers(thick, 1e9) << ",\"silicon\":{\"bulk\":" << num(silicon().k_bulk(T))
      << ",\"in_plane\":" << numbers(in) << ",\"cross_plane\":" << numbers(cross) << ",\"fin\":" << numbers(fin)
      << ",\"points\":[";
    first = 1;
    for (const auto& set : data["silicon_films"]["sets"]) {
        std::string source = set["source"].as<std::string>();
        source = source.substr(0, source.find(','));
        for (const auto& pt : set["points"]) {
            o << (first ? "" : ",") << "{\"source\":" << q(source) << ",\"how\":" << q(set["how"].as<std::string>())
              << ",\"d_nm\":" << num(pt[0].as<double>()) << ",\"k\":" << num(pt[1].as<double>()) << ",\"err\":"
              << num(pt[2].as<double>()) << "}";
            first = 0;
        }
    }
    o << "]},\"interfaces\":[";
    first = 1;
    for (const auto& pt : data["interfaces"]["points"]) {
        std::string a = pt["pair"][0].as<std::string>(), b = pt["pair"][1].as<std::string>();
        double dmm = a == b ? 0.0 : dmm_conductance(a, b, 300.0) / 1e6;
        o << (first ? "" : ",") << "{\"pair\":" << q(a + "/" + b) << ",\"measured\":" << num(pt["tbc"].as<double>())
          << ",\"bound\":" << (pt["bound"] ? q(pt["bound"].as<std::string>()) : "null") << ",\"err\":"
          << (pt["err"] ? num(pt["err"].as<double>()) : "null") << ",\"dmm\":" << (dmm > 0 ? num(dmm) : "null")
          << ",\"source\":" << q(pt["source"].as<std::string>()) << "}";
        first = 0;
    }
    o << "]}";
    return o.str();
}

// the stored atomistic result, the same glass through kALDo when available, and
// the library value against the reference data
std::string wigner_info_json(const std::string& name) {
    const std::string base = data_dir() + "/glass/" + name;
    std::string native = slurp(base + ".wigner.json");
    if (native.empty()) throw std::invalid_argument("no stored result for '" + name + "'");
    std::string kaldo = slurp(base + ".kaldo.json");
    YAML::Node res = YAML::Load(native), ref = YAML::LoadFile(data_dir() + "/published.yaml")["vitreous_silica"];
    Vec temps, ks;
    for (const auto& r : res["temperatures"]) { temps.push_back(r["T"].as<double>()); ks.push_back(r["k_wigner"].as<double>()); }
    auto interp = [&](double t) {
        if (t <= temps.front()) return ks.front();
        if (t >= temps.back()) return ks.back();
        size_t i = size_t(std::upper_bound(temps.begin(), temps.end(), t) - temps.begin());
        return ks[i - 1] + (ks[i] - ks[i - 1]) * (t - temps[i - 1]) / (temps[i] - temps[i - 1]);
    };
    MaterialLibrary lib;
    std::ostringstream o;
    o << "{\"result\":" << native << ",\"kaldo\":" << (kaldo.empty() ? "null" : kaldo) << ",\"reference\":[";
    std::vector<std::tuple<double, double, bool>> rows;
    for (const char* group : {"points", "unverified"})
        for (const auto& pt : ref[group]) rows.emplace_back(pt["T"].as<double>(), pt["k"].as<double>(), std::string(group) == "points");
    std::sort(rows.begin(), rows.end());
    for (size_t i = 0; i < rows.size(); ++i) {
        auto [T, k, verified] = rows[i];
        double model = interp(T), used = lib.props("SiO2_wigner", "bulk", 0, 0, T).kx;
        o << (i ? "," : "") << "{\"T\":" << num(T) << ",\"measured\":" << num(k) << ",\"model\":" << num(model)
          << ",\"deviation\":" << num(model / k - 1) << ",\"library\":" << num(used) << ",\"library_deviation\":"
          << num(used / k - 1) << ",\"verified\":" << (verified ? "true" : "false") << "}";
    }
    o << "],\"samples\":[";
    int first = 1;
    for (const std::string& n : {name, name + "-s2"}) {
        std::string text = slurp(data_dir() + "/glass/" + n + ".wigner.json");
        if (text.empty()) continue;
        for (const auto& r : YAML::Load(text)["temperatures"])
            if (r["T"].as<double>() == 300.0) {
                o << (first ? "" : ",") << "{\"name\":" << q(n) << ",\"k300\":" << num(r["k_wigner"].as<double>()) << "}";
                first = 0;
            }
    }
    o << "],\"library_source\":" << q(lib.table_note("SiO2_wigner")) << "}";
    return o.str();
}

}  // namespace st
