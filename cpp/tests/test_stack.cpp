// Material models, stack voxeliser, sub-stacks, error budget and exports.
// Reference numbers marked "reference" were recorded from the original Python
// implementation, which the C++ core reproduced to solver tolerance.
#include "stacktherm/stacktherm.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>

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
static std::string slurp(const std::string& path) {
    std::ifstream f(path);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}
static Layer film(const std::string& name, double t, const std::string& material) {
    Layer l;
    l.type = Layer::Film; l.name = name; l.thickness = t; l.material = material;
    return l;
}
static std::shared_ptr<MaterialLibrary> lib_with(std::initializer_list<std::pair<const char*, double>> solids) {
    auto lib = std::make_shared<MaterialLibrary>();
    for (auto& s : solids) { MaterialSpec m; m.k = s.second; lib->add(s.first, m); }
    return lib;
}

int main() {
    // ---- material models (reference values)
    const Metal& cu = builtin_metals().at("Cu");
    expect_close("Cu line 20 x 40 nm: conductivity", cu.k_line(20e-9, 40e-9, 300.0), 115.7293351, 1e-8);
    expect_close("Cu line 20 x 40 nm: resistivity", cu.rho_line(20e-9, 40e-9, 300.0), 5.970761099e-08, 1e-8);
    expect_close("Cu film 100 nm: in-plane", cu.k_film(100e-9, 300.0).first, 284.8609054, 1e-8);
    expect_close("Cu film 100 nm: cross-plane", cu.k_film(100e-9, 300.0).second, 228.0218561, 1e-8);
    expect_close("Cu recovers the handbook bulk value", cu.k_line(1e-3, 1e-3, 300.0), 401.0, 1e-3);
    expect_close("Fuchs-Sondheimer at one mean free path", fuchs_sondheimer(1.0), 0.6838565946, 1e-9);
    expect_close("Ru line 10 x 20 nm (short mean free path)", builtin_metals().at("Ru").k_line(10e-9, 20e-9, 300.0),
                 64.15531516, 1e-8);
    expect_close("Si bulk at 300 K", silicon().k_bulk(300.0), 147.5531969, 1e-8);
    expect_close("Si bulk at 500 K", silicon().k_bulk(500.0), 76.10491499, 1e-8);
    expect_close("Si film 100 nm: in-plane", silicon().k_film(100e-9, 300.0).first, 52.33288763, 1e-8);
    expect_close("Si film 100 nm: cross-plane", silicon().k_film(100e-9, 300.0).second, 28.65958003, 1e-8);
    expect_close("Si fin 20 x 20 nm", silicon().k_wire(20e-9, 20e-9, 300.0), 12.27050937, 1e-8);
    expect_close("Si bulk within 2 % of the measured 148 W/m/K", silicon().k_bulk(300.0), 148.0, 0.02);
    expect_close("diffuse mismatch Si/SiO2", dmm_conductance("Si", "SiO2"), 977464158.7, 1e-8);
    expect_close("diffuse mismatch is symmetric", dmm_conductance("SiO2", "Si"), dmm_conductance("Si", "SiO2"), 1e-12);
    expect_close("Cahill-Pohl minimum conductivity of SiO2", cahill_pohl(5953.0, 3743.0, 6.62e28), 1.056014156, 1e-8);

    {   // library: provenance, derived materials, atomistic tables
        MaterialLibrary lib;
        expect("Si/SiO2 interface is a measured value", lib.tbr("Si", "SiO2").kind == "measured");
        expect("unlisted pair falls back to the stated default", lib.tbr("Ni", "SnAg").kind == "assumed");
        expect("same material has no interface resistance", lib.tbr("Cu", "Cu").R == 0.0);
        MaterialSpec ulk;
        ulk.base = "SiCOH";
        ulk.porosity = 0.3;
        lib.add("ULK", ulk);
        expect_close("porous dielectric scales as (1 - phi)^1.5", lib.props("ULK").kx, 0.54 * std::pow(0.7, 1.5), 1e-12);
        expect("derived material inherits its base's interfaces", lib.tbr("Cu", "ULK").kind == "measured");
        expect("size effect can be switched off",
               lib.props("Cu", "wire", 16e-9, 32e-9, 300.0, false).kx == lib.props("Cu").kx);
        double k300 = lib.props("SiO2_wigner", "bulk", 0, 0, 300.0).kx;
        expect("glass conductivity from atoms is within 15 % of the measured 1.38", std::abs(k300 / 1.38 - 1) < 0.15);
        expect("it follows the run temperature", lib.props("SiO2_wigner", "bulk", 0, 0, 400.0).kx > k300);
        bool threw = false;
        try { lib.props("Unobtainium"); } catch (const std::invalid_argument&) { threw = true; }
        expect("unknown material is reported", threw);
    }

    {   // film stack: one-dimensional and exact, with a library interface and a bond entry
        Stack s;
        s.lib = lib_with({{"a", 2.0}, {"b", 0.5}});
        s.lib->set_interface("a", "b", 4e-9);
        s.cell = {1e-6, 1e-6};
        Layer bond;
        bond.type = Layer::Interface; bond.name = "bond"; bond.R = 1e-9;
        s.layers = {film("f1", 100e-9, "a"), film("f2", 50e-9, "b"), bond, film("f3", 80e-9, "b")};
        Model m = s.build();
        expect("film stack needs a single column", m.grid.nx() == 1 && m.grid.ny() == 1);
        Characterization c = characterize(m.grid);
        expect_close("film stack: total resistance", c.R_total, 100e-9 / 2 + 130e-9 / 0.5 + 4e-9 + 1e-9, 1e-9);
        expect_close("film stack: library interface", c.layers[0].R_above, 4e-9, 1e-8);
        expect_close("film stack: bond entry between equal materials", c.layers[1].R_above, 1e-9, 1e-8);
    }

    {   // lines with a barrier: parallel rule along, series rule across
        Stack s;
        s.lib = lib_with({{"metal", 200.0}, {"ild", 0.4}});
        s.lib->set_interface("metal", "ild", 5e-9);
        const double w = 20e-9, p = 50e-9, t = 40e-9;
        s.cell = {p, p};
        Layer l;
        l.type = Layer::Lines; l.name = "m"; l.thickness = t; l.material = "metal"; l.matrix = "ild";
        l.width = w; l.pitch = p; l.direction = 'x';
        l.barrier.present = true; l.barrier.thickness = 2e-9; l.barrier.k = 4.0;
        s.layers = {l};
        Characterization c = characterize(s.build().grid);
        const double core = (w - 4e-9) * (t - 2e-9) / (w * t), k_axial = 200.0 * core + 4.0 * (1 - core);
        expect_close("lines: kx is the parallel rule with the barrier", c.kx, (w * k_axial + (p - w) * 0.4) / p, 1e-8);
        expect_close("lines: ky is the series rule with both walls", c.ky,
                     p / (w / 200.0 + (p - w) / 0.4 + 2 * (5e-9 + 2e-9 / 4.0)), 1e-8);
    }

    {   // a sub-stack reproduces the explicit stack; lumping and flipping
        auto lib = lib_with({{"metal", 150.0}, {"ild", 0.5}});
        lib->set_interface("metal", "ild", 1e-8);
        const double p = 60e-9;
        Layer m1, v1, m2;
        m1.type = Layer::Lines; m1.name = "m1"; m1.thickness = 30e-9; m1.material = "metal"; m1.matrix = "ild";
        m1.width = 20e-9; m1.pitch = p; m1.direction = 'x';
        m2 = m1; m2.name = "m2"; m2.direction = 'y';
        v1.type = Layer::Posts; v1.name = "v1"; v1.thickness = 30e-9; v1.material = "metal"; v1.matrix = "ild";
        v1.size = {20e-9, 20e-9};
        Resolution res;
        res.n_xy = 18; res.min_cells = 3; res.n_z = 3;
        auto total = [&](bool grouped, bool lump, bool flip, Model* out = nullptr) {
            Stack s;
            s.lib = lib; s.cell = {p, p}; s.resolution = res;
            if (grouped) {
                Layer sub;
                sub.type = Layer::SubStack; sub.name = "beol"; sub.cell = {p, p}; sub.layers = {m1, v1, m2};
                sub.lump = lump; sub.flip = flip;
                s.layers = {sub};
            } else {
                s.layers = {m1, v1, m2};
            }
            Model m = s.build();
            Characterization c = characterize(m.grid);
            if (out) *out = m;
            return c;
        };
        Characterization a = total(false, false, false), b = total(true, false, false), c = total(true, true, false);
        expect_close("sub-stack: same total resistance as explicit layers", b.R_total, a.R_total, 1e-7);
        expect_close("sub-stack: same in-plane conductivity", b.kx, a.kx, 1e-7);
        expect_close("lumped sub-stack: same total resistance", c.R_total, a.R_total, 1e-7);
        Model flipped;
        total(true, false, true, &flipped);
        expect("flipped sub-stack reverses the layer order", flipped.grid.layer_names.front() == "beol/m2" &&
                                                              flipped.grid.layer_names.back() == "beol/m1");
    }

    {   // a via on a line crossing conducts; a stray one is flagged
        auto lib = lib_with({{"metal", 150.0}, {"ild", 0.5}});
        const double p = 80e-9;
        auto run = [&](double offset, std::vector<std::string>& warnings) {
            Layer m1, v1, m2;
            m1.type = Layer::Lines; m1.name = "m1"; m1.thickness = 30e-9; m1.material = "metal"; m1.matrix = "ild";
            m1.width = 20e-9; m1.pitch = p / 2; m1.direction = 'x';
            m2 = m1; m2.name = "m2"; m2.direction = 'y';
            v1.type = Layer::Posts; v1.name = "v1"; v1.thickness = 30e-9; v1.material = "metal"; v1.matrix = "ild";
            v1.size = {20e-9, 20e-9}; v1.post_offset = {offset, offset};
            Stack s;
            s.lib = lib; s.cell = {p, p}; s.layers = {m1, v1, m2}; s.resolution.n_xy = 16; s.resolution.n_z = 3;
            Model m = s.build();
            warnings = m.warnings;
            return characterize(m.grid).kz;
        };
        std::vector<std::string> w_on, w_off;
        double on = run(0.25, w_on), off = run(0.5, w_off);
        expect("via on a crossing conducts more than twice as well", on > 2.0 * off);
        expect("aligned via raises no warning", w_on.empty());
        expect("stray via is flagged by name", !w_off.empty() && w_off[0].find("v1") != std::string::npos);
    }

    {   // hot spot: same total power, energy conserved, hotter than uniform heating
        Stack s;
        s.lib = lib_with({{"a", 5.0}});
        s.cell = {1e-6, 1e-6};
        s.top = BC::convection(1e5, 300.0);
        s.resolution.n_xy = 16;
        Layer heater = film("heater", 50e-9, "a");
        heater.heat.present = true; heater.heat.flux = 2e6; heater.heat.has_region = true;
        heater.heat.region = {0.25, 0.5, 0.25, 0.75};
        s.layers = {heater, film("spreader", 400e-9, "a")};
        Result r = simulate(s);
        expect_close("hot spot: total power is the cell-average flux", r.op->power, 2e6 * 1e-12, 1e-10);
        expect("hot spot: energy balance closes", r.op->imbalance < 1e-7);
        expect_close("hot spot: sink-side face temperature", r.op->T_face_top, 300.0 + 2e6 / 1e5, 1e-6);
        expect("hot spot is inside the heated region", r.op->hotspot[0] > 0.25e-6 && r.op->hotspot[0] < 0.5e-6);
        s.layers[0].heat.has_region = false;
        expect("uniform heating runs cooler", simulate(s).op->T_max < r.op->T_max);
    }

    {   // error budget recovers analytic elasticities and leaves the library untouched
        Stack s;
        s.lib = std::make_shared<MaterialLibrary>();
        s.cell = {1e-7, 1e-7};
        const double t = 100e-9, r_bond = 2e-8;
        Layer bond;
        bond.type = Layer::Interface; bond.name = "bond"; bond.R = r_bond;
        s.layers = {film("a", t, "Underfill"), bond, film("b", t, "SiO2")};
        const double ra = t / s.lib->props("Underfill").kx, rb = t / s.lib->props("SiO2").kx;
        const double total = ra + rb + r_bond + s.lib->tbr("Underfill", "SiO2").R;
        Budget b = error_budget(s);
        expect_close("budget: base resistance", b.R_total, total, 1e-8);
        for (const auto& row : b.rows) {
            if (row.key == "estimated_k") expect_close("budget: elasticity of the estimated conductivity", row.elasticity, -(1 - 1 / 1.1) / 0.1 * ra / total, 1e-5);
            if (row.key == "tabulated_k") expect_close("budget: elasticity of the tabulated conductivity", row.elasticity, -(1 - 1 / 1.1) / 0.1 * rb / total, 1e-5);
            if (row.key == "user_tbr") expect_close("budget: elasticity of the stack-file interface", row.elasticity, r_bond / total, 1e-5);
            expect("budget: no metal in this stack", row.key != "metal_k");
        }
        expect("budget: rows are ranked by effect", b.rows.size() >= 3 && b.rows[0].effect >= b.rows[1].effect);
        expect("budget: library scale factors are restored", s.lib->scale.empty());
    }

    {   // geometry errors are reported, not simulated
        auto lib = lib_with({{"a", 1.0}, {"b", 2.0}});
        auto fails = [&](Layer l) {
            Stack s;
            s.lib = lib; s.cell = {100e-9, 100e-9}; s.layers = {l};
            try { s.build(); } catch (const std::invalid_argument&) { return true; }
            return false;
        };
        Layer l;
        l.type = Layer::Lines; l.name = "m"; l.thickness = 10e-9; l.material = "a"; l.matrix = "b"; l.width = 10e-9;
        l.pitch = 30e-9;
        expect("pitch that does not divide the cell is rejected", fails(l));
        l.pitch = 50e-9; l.barrier.present = true; l.barrier.thickness = 6e-9;
        expect("barrier thicker than the conductor is rejected", fails(l));
    }

    {   // the shipped examples (reference totals) and the export files
        const std::string root = std::string(STACKTHERM_SOURCE_DIR);
        Stack hybrid = load_stack(root + "/examples/hybrid_bond.yaml");
        Result r = simulate(hybrid);
        expect_close("hybrid-bond example: total resistance", r.ch.R_total, 8.430124065e-07, 1e-6);
        expect_close("hybrid-bond example: hottest temperature", r.op->T_max, 338.421511, 1e-7);
        expect("hybrid-bond example: 28 layers, no warnings", r.ch.layers.size() == 28 && r.model.warnings.empty());
        double sum = 0, pad_kz = 0, pad_mix = 0;
        for (const auto& p : r.ch.layers) { sum += p.R + p.R_above; if (p.name == "pad-b") { pad_kz = p.kz; pad_mix = p.kz_mix; } }
        expect_close("layer and interface terms add up to the total", sum, r.ch.R_total, 1e-8);
        expect("the pad is throttled by spreading in the wiring below", pad_kz < 0.3 * pad_mix);
        const std::string dir = root + "/build/test_out";
        auto files = write_all(r, dir, true);
        std::string json = slurp(dir + "/result.json"), apdl = slurp(dir + "/layers_apdl.mac"), vtr = slurp(dir + "/field.vtr");
        expect("export: five files written", files.size() == 5);
        expect("export: JSON carries measured-interface provenance", json.find("\"provenance\":\"measured") != std::string::npos);
        expect("export: JSON lists the size-dependent metal model", json.find("FS+MS wire") != std::string::npos);
        size_t cards = 0;
        for (size_t pos = 0; (pos = apdl.find("MP,KZZ", pos)) != std::string::npos; ++pos) ++cards;
        expect("export: one APDL card per layer", cards == r.ch.layers.size());
        expect("export: VTK file holds the temperature array", vtr.find("Name=\"T_K\"") != std::string::npos &&
                                                               vtr.size() > 4 * r.model.grid.n());
        expect("export: report template was filled", slurp(dir + "/report.html").find("__DATA__") == std::string::npos);

        Stack bump = load_stack(root + "/examples/microbump.yaml");
        expect_close("microbump example: total resistance", simulate(bump, false).ch.R_total, 3.998675609e-06, 1e-6);
        Stack coarse = hybrid;
        coarse.resolution.n_xy = 16; coarse.resolution.n_z = 3;
        auto runs = convergence(coarse, {1.0, 2.0});
        expect("mesh refinement changes the total by under 2 %", runs[0].change_R < 0.02);
        auto sweep = temperature_sweep(coarse, {300.0, 400.0});
        expect("resistance rises with temperature", sweep[400.0].R_total > sweep[300.0].R_total);
        double bulk_r = through_plane(bulk_variant(coarse).build().grid).r_total;
        expect("the bulk-property model underestimates the resistance", bulk_r < 0.9 * simulate(coarse, false).ch.R_total);
    }

    Check analytic = validate_analytic();
    expect("closed-form checks pass (manufactured solution, Rayleigh cylinders)", analytic.ok);
    Check published = validate_published();
    expect("published-data checks pass", published.ok);
    std::printf("\n%s\n", failures ? "FAILED" : "all passed");
    return failures ? 1 : 0;
}
