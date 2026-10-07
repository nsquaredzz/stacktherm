// Command line: stacktherm run | validate | init
#include "stacktherm/stacktherm.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace st;

static const char* TEMPLATE = R"(# stacktherm stack file.  Layers are listed bottom to top inside one laterally
# periodic unit cell.  Lengths take units (nm, um), power W/cm2, interfaces
# either a conductance (MW/m2K) or a resistance (m2K/GW).
name: my-stack
cell: [200 nm, 200 nm]
temperature: 350 K              # material models are evaluated here
resolution: {n_xy: 32, n_z: 4}

boundary:
  top: {type: convection, h: 5.0e4, T: 318 K}
  bottom: {type: adiabatic}

layers:
  - {name: substrate, type: film, material: Si, thickness: 2 um}
  - {name: devices, type: film, material: Si, thickness: 50 nm, heat: {flux: 100 W/cm2}}
  - {name: M1, type: lines, material: Cu, matrix: SiCOH, width: 50 nm, pitch: 100 nm,
     direction: x, thickness: 100 nm, barrier: {thickness: 3 nm}}
  # offset is the via centre as a fraction of the cell: put it on a line crossing
  - {name: V1, type: posts, material: Cu, matrix: SiCOH, size: 50 nm, offset: 0.25,
     thickness: 100 nm, barrier: {thickness: 3 nm}}
  - {name: M2, type: lines, material: Cu, matrix: SiCOH, width: 50 nm, pitch: 100 nm,
     direction: y, thickness: 100 nm, barrier: {thickness: 3 nm}}
  - {name: passivation, type: film, material: SiO2_wigner, thickness: 200 nm}
)";

static const char* USAGE =
    "usage: stacktherm <command> [options]\n\n"
    "  run STACK.yaml [-o DIR] [--report] [--budget] [--converge] [--temps K...] [--refine F] [--no-operating]\n"
    "      solve a stack file; -o writes result.json, layers.csv, APDL cards and a VTK field\n"
    "  validate          closed-form and published-data checks\n"
    "  wigner [--model NAME] [--recompute] [--sample N]\n"
    "      conductivity of a stored glass model from atoms (native Wigner engine)\n"
    "  glass --name NAME [--seed N] [--quench-ps P] [--melt-ps P] [--cells N]\n"
    "      melt, quench and relax a new amorphous SiO2 sample\n"
    "  init [PATH]       write a commented starter stack file\n"
    "  --version\n";

static int run(const std::vector<std::string>& a) {
    std::string path, out;
    bool report = false, budget = false, converge = false, operating = true;
    double refine = 1.0;
    Vec temps;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == "-o" || a[i] == "--output") { if (++i >= a.size()) throw std::invalid_argument("-o needs a directory"); out = a[i]; }
        else if (a[i] == "--report") report = true;
        else if (a[i] == "--budget") budget = true;
        else if (a[i] == "--converge") converge = true;
        else if (a[i] == "--no-operating") operating = false;
        else if (a[i] == "--refine") { if (++i >= a.size()) throw std::invalid_argument("--refine needs a factor"); refine = std::stod(a[i]); }
        else if (a[i] == "--temps") { while (i + 1 < a.size() && a[i + 1][0] != '-') temps.push_back(std::stod(a[++i])); }
        else if (a[i][0] == '-') throw std::invalid_argument("unknown option " + a[i]);
        else path = a[i];
    }
    if (path.empty()) throw std::invalid_argument("run needs a stack file");
    Stack stack = load_stack(path);
    if (refine != 1.0) stack.resolution = stack.resolution.refined(refine);
    Result result = simulate(stack, operating);
    std::cout << result.summary();
    if (budget) std::cout << "\n" << error_budget(stack).text();
    if (converge) {
        std::cout << "\nmesh convergence (largest relative change of any layer k against the finest mesh):\n";
        for (const auto& r : convergence(stack))
            std::printf("  x%-4g %10zu cells   change %.2e   R_total change %.2e   %.1f s\n", r.factor, r.cells, r.change,
                        r.change_R, r.wall_s);
    }
    std::map<double, Characterization> sweep;
    if (!temps.empty()) {
        sweep = temperature_sweep(stack, temps);
        std::cout << "\neffective through-plane resistance vs temperature:\n";
        for (const auto& kv : sweep)
            std::printf("  %7.1f K   R_total %.4g mm2K/W   kx %.3g  kz %.3g W/m/K\n", kv.first, kv.second.R_total * 1e6,
                        kv.second.kx, kv.second.kz);
    }
    if (!out.empty()) {
        std::cout << "\nwrote";
        for (const auto& f : write_all(result, out, report, temps.empty() ? nullptr : &sweep)) std::cout << " " << f;
        std::cout << "\n";
    }
    return 0;
}

static void say(const std::string& m) { std::cout << "  " << m << std::endl; }

static int wigner(const std::vector<std::string>& a) {
    std::string model = "a-SiO2-648";
    bool recompute = false;
    int sample = 160;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == "--model" && i + 1 < a.size()) model = a[++i];
        else if (a[i] == "--recompute") recompute = true;
        else if (a[i] == "--sample" && i + 1 < a.size()) sample = std::stoi(a[++i]);
        else throw std::invalid_argument("unknown option " + a[i]);
    }
    if (recompute) wigner_compute(model, {50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800}, sample, say);
    std::cout << wigner_report(model);
    return 0;
}

static int glass(const std::vector<std::string>& a) {
    std::string name;
    unsigned seed = 1;
    double quench = 60.0, melt = 8.0;
    int cells = 3;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == "--name" && i + 1 < a.size()) name = a[++i];
        else if (a[i] == "--seed" && i + 1 < a.size()) seed = unsigned(std::stoul(a[++i]));
        else if (a[i] == "--quench-ps" && i + 1 < a.size()) quench = std::stod(a[++i]);
        else if (a[i] == "--melt-ps" && i + 1 < a.size()) melt = std::stod(a[++i]);
        else if (a[i] == "--cells" && i + 1 < a.size()) cells = std::stoi(a[++i]);
        else throw std::invalid_argument("unknown option " + a[i]);
    }
    if (name.empty()) throw std::invalid_argument("glass needs --name");
    Glass g = relax(melt_quench(cristobalite(cells), 6000.0, 300.0, melt, quench, 2.0, 1.0, seed, say), say);
    std::vector<int> pi, pj;
    g.pairs(10.17, pi, pj);
    Vec f;
    double e = g.energy_forces(pi, pj, g.pos, f), fmax = 0;
    for (double v : f) fmax = std::max(fmax, std::abs(v));
    Coordination c = g.coordination();
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "{\"atoms\": %d, \"box_A\": %.12g, \"density_g_cm3\": %.6g, \"energy_eV_per_atom\": %.8g, "
                  "\"max_force_eV_A\": %.3g, \"melt_K\": 6000.0, \"melt_ps\": %g, \"quench_ps\": %g, "
                  "\"quench_rate_K_per_s\": %.4g, \"seed\": %u, \"si_four_fold\": %.6g, \"o_two_fold\": %.6g, "
                  "\"si_o_bond\": %.6g, \"si_o_bond_std\": %.6g}",
                  g.n(), g.box, g.density(), e / g.n(), fmax, melt, quench, 5700.0 / (quench * 1e-12), seed,
                  c.si_four_fold, c.o_two_fold, c.si_o_bond, c.si_o_bond_std);
    g.info = buf;
    save_glass(g, name);
    std::cout << "stored glass '" << name << "': " << buf << "\n";
    return 0;
}

static int validate() {
    std::cout << "Solver verification against closed-form solutions\n\n";
    Check a = validate_analytic();
    std::cout << a.text << "\n\nComparison with published measurements\n\n";
    Check p = validate_published();
    std::cout << p.text << "\nRESULT: " << (a.ok && p.ok ? "all checks passed" : "CHECKS FAILED") << "\n";
    return a.ok && p.ok ? 0 : 1;
}

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    try {
        if (args.empty() || args[0] == "-h" || args[0] == "--help") { std::cout << USAGE; return args.empty() ? 2 : 0; }
        if (args[0] == "--version") { std::cout << "stacktherm 0.2.0\n"; return 0; }
        std::vector<std::string> rest(args.begin() + 1, args.end());
        if (args[0] == "run") return run(rest);
        if (args[0] == "validate") return validate();
        if (args[0] == "wigner") return wigner(rest);
        if (args[0] == "glass") return glass(rest);
        if (args[0] == "init") {
            std::string path = rest.empty() ? "stack.yaml" : rest[0];
            if (std::ifstream(path).good()) throw std::invalid_argument(path + " already exists");
            std::ofstream(path) << TEMPLATE;
            std::cout << "wrote " << path << "; run it with: stacktherm run " << path << " -o out --report\n";
            return 0;
        }
        throw std::invalid_argument("unknown command '" + args[0] + "'\n\n" + USAGE);
    } catch (const std::exception& e) {
        // problems with the input, not with the program: say what and stop
        std::cerr << "error: " << e.what() << "\n";
        return 2;
    }
}
