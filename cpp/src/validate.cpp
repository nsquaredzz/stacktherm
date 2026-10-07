// Checks against closed-form solutions and against published measurements
// (data/published.yaml records each number and where it came from).
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <cstdio>
#include <sstream>

#include <yaml-cpp/yaml.h>

namespace st {

// Transverse k_eff / k_matrix of a square array of circular cylinders: Rayleigh's
// multipole series as truncated by Perrins, McKenzie & McPhedran, Proc. R. Soc. A
// 369, 207 (1979).
double rayleigh_cylinders(double f, double sigma) {
    if (sigma == 1.0) return 1.0;
    double t = (1.0 + sigma) / (1.0 - sigma);
    return 1.0 - 2.0 * f / (t + f - 0.305827 * std::pow(f, 4) * t / (t * t - 1.402958 * std::pow(f, 8)) -
                            0.013362 * std::pow(f, 8) / t);
}

// the same cylinder array through the regular stack voxeliser (round posts, staircase walls)
double cylinder_array_fv(double f, double sigma, int n, double k_matrix, double radius, double r_int) {
    Stack s;
    s.lib = std::make_shared<MaterialLibrary>();
    MaterialSpec in, mat;
    in.k = sigma * k_matrix;
    mat.k = k_matrix;
    s.lib->add("inclusion", in);
    s.lib->add("matrix", mat);
    s.lib->set_interface("inclusion", "matrix", r_int);
    const double length = radius * std::sqrt(kPi / f);
    s.cell = {length, length};
    Layer posts;
    posts.type = Layer::Posts; posts.name = "posts"; posts.thickness = length / n;
    posts.material = "inclusion"; posts.matrix = "matrix"; posts.size = {2 * radius, 2 * radius}; posts.circle = true;
    s.layers = {posts};
    s.resolution.n_xy = n; s.resolution.min_cells = 1; s.resolution.n_z = 1; s.resolution.circle_cells = 1;
    return characterize(s.build().grid).kx;
}

// max-norm error for an anisotropic manufactured solution on a non-uniform grid:
// T = cos(2 pi x / Lx) cos(2 pi y / Ly) cos(pi z / 2H)
double manufactured_error(int n) {
    const double lx = 1.0e-6, ly = 1.5e-6, h = 0.4e-6, kx = 3.0, ky = 7.0, kz = 0.5;
    auto stretched = [&](double length) {
        Vec e(n + 1);
        for (int i = 0; i <= n; ++i) { double s = double(i) / n; e[i] = length * (s + 0.12 * std::sin(2 * kPi * s) / (2 * kPi)); }
        return e;
    };
    Grid g;
    g.xe = stretched(lx); g.ye = stretched(ly); g.ze = stretched(h);
    g.allocate();
    const double a = 2 * kPi / lx, b = 2 * kPi / ly, c = kPi / (2 * h);
    Vec exact(g.n());
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                size_t id = g.idx(i, j, k);
                double x = 0.5 * (g.xe[i] + g.xe[i + 1]), y = 0.5 * (g.ye[j] + g.ye[j + 1]), z = 0.5 * (g.ze[k] + g.ze[k + 1]);
                exact[id] = std::cos(c * z) * std::cos(b * y) * std::cos(a * x);
                g.kx[id] = kx; g.ky[id] = ky; g.kz[id] = kz;
                g.q[id] = (kx * a * a + ky * b * b + kz * c * c) * exact[id];
            }
    Operating op = operate(g, BC(), BC::dirichlet(0.0));
    double err = 0;
    for (size_t i = 0; i < exact.size(); ++i) err = std::max(err, std::abs(op.T[i] - exact[i]));
    return err;
}

Check validate_analytic() {
    Check c;
    char buf[240];
    double e16 = manufactured_error(16), e32 = manufactured_error(32), e64 = manufactured_error(64);
    double order = std::log2(e32 / e64);
    std::snprintf(buf, sizeof buf, "  manufactured solution, max error at 16/32/64 cells: %.2e / %.2e / %.2e   "
                                   "observed order %.2f\n", e16, e32, e64, order);
    c.text += buf;
    c.ok = order > 1.9;
    for (auto [sigma, tol] : {std::pair<double, double>{1000.0, 4e-3}, {0.01, 1e-2}}) {
        double err = std::abs(cylinder_array_fv(0.3, sigma, 128) / rayleigh_cylinders(0.3, sigma) - 1.0);
        std::snprintf(buf, sizeof buf, "  cylinder array vs Rayleigh series, k ratio %g: error %.2e\n", sigma, err);
        c.text += buf;
        c.ok = c.ok && err < tol;
    }
    return c;
}

namespace {

// Cross-plane conductivity of one metal level (lines plus the via-level dielectric
// under them) inside a stack of crossed line levels.  first_order = false switches
// off the Cu size effect and the boundary resistance: the bulk-property model.
double beol_layer_kz(double cu_fraction, double k_dielectric, double tbr, double h_metal, double h_total,
                     double linewidth, bool first_order = true) {
    const double pitch = linewidth / (cu_fraction * h_total / h_metal);
    Stack s;
    s.lib = std::make_shared<MaterialLibrary>();
    MaterialSpec ild;
    ild.k = k_dielectric;
    s.lib->add("ild", ild);
    s.lib->set_interface("Cu", "ild", first_order ? tbr : 0.0);
    s.cell = {pitch, pitch};
    for (int n = 0; n < 3; ++n) {   // three identical levels; the middle one sees realistic neighbours
        Layer film, lines;
        film.type = Layer::Film; film.name = "ild" + std::to_string(n); film.thickness = h_total - h_metal; film.material = "ild";
        lines.type = Layer::Lines; lines.name = "metal" + std::to_string(n); lines.thickness = h_metal;
        lines.material = "Cu"; lines.matrix = "ild"; lines.width = linewidth; lines.pitch = pitch;
        lines.direction = n % 2 ? 'y' : 'x'; lines.size_effect = first_order;
        s.layers.push_back(film);
        s.layers.push_back(lines);
    }
    Layer cap;
    cap.type = Layer::Film; cap.name = "cap"; cap.thickness = h_total - h_metal; cap.material = "ild";
    s.layers.push_back(cap);
    s.resolution.n_xy = 24; s.resolution.min_cells = 4; s.resolution.n_z = 6;
    Characterization c = characterize(s.build().grid);
    const LayerProps& film = c.layers[2];
    const LayerProps& metal = c.layers[3];
    return h_total / (film.R + film.R_above + metal.R + metal.R_above);
}

}  // namespace

Check validate_published() {
    Check c;
    char buf[400];
    YAML::Node data = YAML::LoadFile(data_dir() + "/published.yaml");

    c.text += "Silicon films, in-plane k at 300 K (model calibrated on bulk k(T) only)\n\n";
    std::snprintf(buf, sizeof buf, "  %8s%10s%6s%8s%8s  source\n", "d [nm]", "measured", "+/-", "model", "dev");
    c.text += buf;
    double sum = 0, worst = 0;
    int count = 0;
    for (const auto& set : data["silicon_films"]["sets"]) {
        std::string source = set["source"].as<std::string>();
        std::string label = source.substr(0, source.find(','));
        bool figure = set["how"].as<std::string>() == "figure";
        for (const auto& p : set["points"]) {
            double d = p[0].as<double>(), k = p[1].as<double>(), err = p[2].as<double>();
            double model = silicon().k_film(d * 1e-9, 300.0).first, dev = model / k - 1.0;
            std::snprintf(buf, sizeof buf, "  %8g%10g%6g%8.1f%+7.0f%%  %s%s\n", d, k, err, model, dev * 100, label.c_str(),
                          figure ? " (read from figure)" : "");
            c.text += buf;
            if (label.compare(0, 5, "Cuffe") == 0) {
                sum += dev;
                ++count;
                if (std::abs(dev) > std::abs(worst)) worst = dev;
            }
        }
    }
    std::snprintf(buf, sizeof buf, "\n  Cuffe et al. membranes: mean deviation %+.1f %%, worst %+.1f %%\n", sum / count * 100,
                  worst * 100);
    c.text += buf;
    bool ok_si = std::abs(sum / count) < 0.12 && std::abs(worst) < 0.20;

    const YAML::Node b = data["beol_layers"], rep = b["reported"], ass = b["assumed"];
    c.text += "\n\nBEOL metal levels, cross-plane k (W/m/K), 22 nm node, TDTR\n  " + b["source"].as<std::string>() + "\n\n";
    std::snprintf(buf, sizeof buf, "  %-6s%8s%14s%14s%8s%17s%12s\n", "layer", "Cu frac", "TDTR 5x", "TDTR 20x", "model",
                  "linewidth span", "bulk model");
    c.text += buf;
    const double tbr = rep["tbr_cu_dielectric"].as<double>(), hm = ass["h_metal_nm"].as<double>() * 1e-9,
                 ht = ass["h_total_nm"].as<double>() * 1e-9, width = ass["linewidth_nm"].as<double>() * 1e-9;
    int inside = 0, total = 0;
    for (const auto& lay : rep["layers"]) {
        double cu = lay["cu_fraction"].as<double>(), kd = lay["k_dielectric"].as<double>();
        double model = beol_layer_kz(cu, kd, tbr, hm, ht, width);
        double lo = beol_layer_kz(cu, kd, tbr, hm, ht, ass["linewidth_range_nm"][0].as<double>() * 1e-9);
        double hi = beol_layer_kz(cu, kd, tbr, hm, ht, ass["linewidth_range_nm"][1].as<double>() * 1e-9);
        double bulk = beol_layer_kz(cu, kd, tbr, hm, ht, width, false);
        double k5 = lay["k_5x"].as<double>(), e5 = lay["err_5x"].as<double>(), k20 = lay["k_20x"].as<double>(),
               e20 = lay["err_20x"].as<double>();
        std::snprintf(buf, sizeof buf, "  %-6s%8.2f%8.2f +/-%.2f%8.2f +/-%.2f%8.2f%9.2f-%.2f%12.2f\n",
                      lay["name"].as<std::string>().c_str(), cu, k5, e5, k20, e20, model, std::min(lo, hi),
                      std::max(lo, hi), bulk);
        c.text += buf;
        inside += std::abs(model - k5) <= e5 || std::abs(model - k20) <= e20;
        ++total;
    }
    std::snprintf(buf, sizeof buf,
                  "\n  %d of %d layers inside the experimental error bars.\n"
                  "  Inputs taken from the paper: Cu fraction, dielectric k and Cu/dielectric TBR (the\n"
                  "  authors' own fit).  Line height, level height and linewidth are estimated here;\n"
                  "  'linewidth span' shows the model over 40-100 nm.  'bulk model' drops the Cu size\n"
                  "  effect and the boundary resistance.\n", inside, total);
    c.text += buf;
    bool ok_beol = inside >= total - 1;

    c.text += "\n\nInterface conductance, MW/m2K: measurement vs diffuse mismatch model\n\n";
    for (const auto& p : data["interfaces"]["points"]) {
        std::string a = p["pair"][0].as<std::string>(), bb = p["pair"][1].as<std::string>();
        std::ostringstream meas;
        if (p["bound"]) meas << ">= ";
        meas << p["tbc"].as<double>();
        if (p["err"]) meas << " +/- " << p["err"].as<double>();
        double dmm = a == bb ? 0.0 : dmm_conductance(a, bb, 300.0) / 1e6;
        char dm[32];
        if (dmm > 0) std::snprintf(dm, sizeof dm, "%.0f", dmm); else std::snprintf(dm, sizeof dm, "n/a");
        std::snprintf(buf, sizeof buf, "  %-11s%12s   DMM %5s   %s\n", (a + "/" + bb).c_str(), meas.str().c_str(), dm,
                      p["source"].as<std::string>().c_str());
        c.text += buf;
    }
    c.text += "\n  The DMM is consistent with the Si/SiO2 lower limit but overestimates the weakly\n"
              "  bonded Cu/SiO2 interface roughly tenfold, which is why the library prefers\n"
              "  measured values and labels everything else as assumed.\n";

    MaterialLibrary lib;
    c.text += "\n\nAmorphous SiO2 from atoms: material SiO2_wigner as used by the solver\n  (" +
              lib.table_note("SiO2_wigner") + ")\n\n";
    std::snprintf(buf, sizeof buf, "  %7s%10s%9s%8s\n", "T [K]", "measured", "library", "dev");
    c.text += buf;
    bool ok_glass = true;
    const YAML::Node vs = data["vitreous_silica"];
    for (const char* group : {"points", "unverified"})
        for (const auto& p : vs[group]) {
            double T = p["T"].as<double>(), k = p["k"].as<double>(), used = lib.props("SiO2_wigner", "bulk", 0, 0, T).kx;
            bool verified = std::string(group) == "points";
            std::snprintf(buf, sizeof buf, "  %7.0f%10.2f%9.3f%+7.0f%%%s\n", T, k, used, (used / k - 1) * 100,
                          verified ? "" : "  (handbook value recalled, not re-fetched)");
            c.text += buf;
            if (verified) ok_glass = ok_glass && std::abs(used / k - 1) < 0.15;
        }
    c.ok = ok_si && ok_beol && ok_glass;
    return c;
}

}  // namespace st
