// stacktherm core: 3D thermal simulation of chip-stack unit cells with
// size-dependent material inputs and interface resistance.  All quantities SI
// unless a name says otherwise.
#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace st {

using Vec = std::vector<double>;
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kPi = 3.14159265358979323846;

// ---------------------------------------------------------------- numerics
// "42 nm", "100 W/cm2", "150 MW/m2K", "85 C" -> SI.  A bare number is scaled
// by default_unit when one is given.
double parse_quantity(const std::string& text, const std::string& default_unit = "");
// "200 MW/m2K" (a conductance) or "5e-9" / "5 m2K/GW" (a resistance) -> m2K/W
double parse_interface(const std::string& text);
void gauss_legendre(int n, double a, double b, Vec& x, Vec& w);
double expint(int n, double x);                       // exponential integral E_n(x)
double voigt(double x, double sigma, double gamma);   // unit-area Voigt profile
std::string data_dir();                               // library, published data, glass models
std::string json_escape(const std::string& s);
std::string fmt(double v, int significant = 6);

// ---------------------------------------------------------------- grid
// Rectilinear voxel grid.  Cell (i, j, k) lives at index (k * ny + j) * nx + i.
struct Grid {
    Vec xe, ye, ze;                    // cell edges
    Vec kx, ky, kz;                    // conductivity per cell, W/m/K
    Vec rx, ry, rz;                    // interface resistance on the +x/+y/+z face, m2K/W
    Vec q, rhoc;                       // heat source W/m3, heat capacity J/m3/K
    std::vector<int> mat;              // material-instance id per cell
    std::vector<int> layer_of_k;       // layer index of each z cell
    std::vector<std::string> layer_names;

    int nx() const { return int(xe.size()) - 1; }
    int ny() const { return int(ye.size()) - 1; }
    int nz() const { return int(ze.size()) - 1; }
    size_t n() const { return size_t(nx()) * ny() * nz(); }
    size_t idx(int i, int j, int k) const { return (size_t(k) * ny() + j) * nx() + i; }
    double dx(int i) const { return xe[i + 1] - xe[i]; }
    double dy(int j) const { return ye[j + 1] - ye[j]; }
    double dz(int k) const { return ze[k + 1] - ze[k]; }
    double lx() const { return xe.back() - xe.front(); }
    double ly() const { return ye.back() - ye.front(); }
    double lz() const { return ze.back() - ze.front(); }
    void allocate();                   // size the per-cell arrays (zero-filled) from the edges
    std::vector<std::pair<int, int>> layer_slices() const;   // [k0, k1) of every layer
};

Vec graded_edges(double t, double h_bot, double h_top, double ratio = 1.4, double h_max = -1.0);

// ---------------------------------------------------------------- finite volumes
struct BC {
    enum Kind { Adiabatic, Dirichlet, Convection, Flux } kind = Adiabatic;
    double T = 0.0, h = kInf, q = 0.0, r = 0.0;
    bool fixes_temperature() const { return kind == Dirichlet || kind == Convection; }
    double resistance() const { return r + (kind == Convection ? 1.0 / h : 0.0); }
    static BC dirichlet(double T) { BC b; b.kind = Dirichlet; b.T = T; return b; }
    static BC convection(double h, double T) { BC b; b.kind = Convection; b.h = h; b.T = T; return b; }
    static BC flux(double q) { BC b; b.kind = Flux; b.q = q; return b; }
};

struct Csr {
    int n = 0;
    std::vector<ptrdiff_t> ptr, col;
    Vec val;
    void mul(const Vec& x, Vec& y) const;
};

struct System {
    const Grid* grid = nullptr;
    Csr A;
    Vec b;
    Vec gx, gy, gz;        // face conductances; gx[.., nx-1] and gy[.., ny-1, ..] are the wrap faces
    Vec gbot, gtop;        // per column
    BC bottom, top;
    double jump[2] = {0.0, 0.0};
    bool pinned = false;

    void face_flows(const Vec& T, Vec& fx, Vec& fy, Vec& fz) const;     // W, positive towards +axis
    void boundary_flows(const Vec& T, Vec& qbot, Vec& qtop) const;      // W leaving, per column
    void cell_flux(const Vec& T, Vec& qx, Vec& qy, Vec& qz) const;      // W/m2 at cell centres
};

System assemble(const Grid& g, const BC& bottom = BC(), const BC& top = BC(), bool periodic_x = true,
                bool periodic_y = true, double jump_x = 0.0, double jump_y = 0.0, bool sources = true);

struct SolveInfo { int iterations = 0; double residual = 0.0, setup_s = 0.0, solve_s = 0.0; };

// Conjugate gradients preconditioned by algebraic multigrid; the hierarchy is
// built once and reused for every right-hand side.
class LinearSolver {
public:
    explicit LinearSolver(const Csr& A, double tol = 1e-10, int maxiter = 500);
    ~LinearSolver();
    Vec solve(const Vec& b, SolveInfo* info = nullptr) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ---------------------------------------------------------------- cell problems
struct LayerProps {
    std::string name;
    double thickness = 0, kx = 0, ky = 0, kz = 0, R = 0, R_above = 0, rhoc = 0, kz_mix = 0;
};

struct Characterization {
    std::vector<LayerProps> layers;
    double kx = 0, ky = 0, kz = 0, R_total = 0, thickness = 0;
    std::map<std::string, SolveInfo> solves;
    std::string table() const;
};

struct ThroughPlane { Vec r_layer, r_int; double r_total = 0; SolveInfo info; Vec T; };
ThroughPlane through_plane(const Grid& g, double tol = 1e-10);
Characterization characterize(const Grid& g, double tol = 1e-10);

struct Operating {
    Vec T;
    double power = 0, q_bottom = 0, q_top = 0, imbalance = 0, T_max = 0;
    std::array<double, 3> hotspot{};
    Vec layer_T_mean, layer_T_max, layer_T_min;
    double T_face_bottom = 0, T_face_top = 0;
    SolveInfo info;
};
Operating operate(const Grid& g, const BC& bottom, const BC& top, double tol = 1e-10);

// ---------------------------------------------------------------- material models
double fuchs_sondheimer(double kappa, double p = 0.0);   // film conductivity ratio, in-plane
double mayadas_shatzkes(double alpha);                   // grain-boundary resistivity ratio
double fs_phonon(double delta);                          // diffuse film suppression for phonons

struct Metal {
    std::string name;
    double rho0 = 0, tcr = 0, rho_lambda = 0, k_bulk = 0, p = 0.0, R = 0.3, grain_factor = 1.0,
           density = 0, cp = 0;
    double rho_bulk(double T) const { return rho0 * (1.0 + tcr * (T - 293.15)); }
    double mfp(double T) const { return rho_lambda / rho_bulk(T); }
    double lorenz() const { return k_bulk * rho_bulk(300.0) / 300.0; }
    double rho_line(double w, double h, double T) const;
    double rho_film(double d, double T) const;
    double rho_film_cross(double d, double T) const;
    double k_bulk_T(double T) const { return lorenz() * T / rho_bulk(T); }
    double k_line(double w, double h, double T) const { return lorenz() * T / rho_line(w, h, T); }
    std::pair<double, double> k_film(double d, double T) const;   // in-plane, cross-plane
};
const std::map<std::string, Metal>& builtin_metals();

struct PhononCrystal {
    std::string name;
    double v_s = 0, n_prim = 0, A = 0, P = 0, C_U = 0, density = 0, cp = 0;
    double k_bulk(double T) const;
    std::pair<double, double> k_film(double d, double T) const;   // in-plane, cross-plane
    double k_wire(double w, double h, double T) const;
};
const PhononCrystal& silicon();

double dmm_conductance(const std::string& a, const std::string& b, double T = 300.0);   // 0 if unknown
double cahill_pohl(double v_l, double v_t, double n_atoms, double T = 300.0);

struct Props { double kx = 0, ky = 0, kz = 0, rhoc = 0; std::string base, note; };
struct Tbr { double R = 0; std::string kind, source; };

struct MaterialSpec {                // a custom material from a stack file
    std::optional<double> k, density, cp, porosity, p, R, grain_factor;
    std::string base, source;
    std::vector<std::string> tables;
};

class MaterialLibrary {
public:
    MaterialLibrary();               // loads <data_dir>/library.yaml
    void add(const std::string& name, const MaterialSpec& spec);
    bool add_table(const std::string& name, const std::vector<std::string>& models,
                   const std::string& base = "", double density = 0, double cp = 0);
    void set_interface(const std::string& a, const std::string& b, double resistance);
    void without_interfaces();
    std::string canonical(const std::string& name) const;
    // shape: "bulk", "film" (dims = thickness) or "wire" (dims = width, height)
    Props props(const std::string& name, const std::string& shape = "bulk", double d0 = 0, double d1 = 0,
                double T = 300.0, bool size_effect = true) const;
    Tbr tbr(const std::string& a, const std::string& b, double T = 300.0) const;
    double scale_of(const std::string& key) const;
    std::string table_note(const std::string& name) const;

    std::map<std::string, double> scale;   // per input class, for the error budget
    std::map<std::string, Metal> metals;
    double default_tbr = 1e-8;
    std::string fallback = "assumed";
private:
    struct Solid { double k = 0, density = 0, cp = 0; std::string source; };
    struct Table { Vec T, k; double rhoc = 0; std::string note; };
    std::map<std::string, Solid> solids_;
    std::map<std::string, Table> tables_;
    std::map<std::string, std::string> base_;
    std::map<std::pair<std::string, std::string>, std::pair<double, std::string>> measured_;
    std::map<std::pair<std::string, std::string>, double> user_tbr_;
};

// ---------------------------------------------------------------- stack description
struct Resolution {
    int n_xy = 48, min_cells = 2, n_z = 4;
    double grading = 1.3;
    int circle_cells = 10;
    Resolution refined(double f) const;
};

struct Barrier {
    bool present = false;
    double thickness = 0, k = 3.0, R_extra = 0;
    bool bottom = true;
    double R() const { return thickness / k + R_extra; }
};

struct Heat {
    bool present = false, has_region = false;
    double flux = 0;                     // W/m2 averaged over the unit cell
    std::array<double, 4> region{};      // x0, x1, y0, y1 as fractions of the cell
};

struct Layer {
    enum Type { Film, Aniso, Lines, Posts, Interface, SubStack } type = Film;
    std::string name;
    double thickness = 0;
    std::string material, matrix;
    double kx = 0, ky = 0, kz = 0, rhoc = 0;                 // Aniso
    double width = 0, pitch = 0, offset = 0.5;               // Lines
    char direction = 'x';
    std::array<double, 2> size{}, post_pitch{}, post_offset{0.5, 0.5};   // Posts
    bool has_pitch = false, circle = false;
    Barrier barrier;
    bool inverted = false, size_effect = true;
    Heat heat;
    int n_z = 0;                                             // 0: use the resolution default
    double R = 0;                                            // Interface: added to every pair
    std::vector<std::pair<std::string, double>> pairs;       // Interface: "A|B" -> resistance
    std::array<double, 2> cell{};                            // SubStack
    std::vector<Layer> layers;
    bool flip = false, lump = false;
    std::optional<Resolution> resolution;
};

struct Instance {
    std::string name, layer, material, base, note;
    double kx = 0, ky = 0, kz = 0, fraction = 1.0;
    double shell_side = 0, shell_bottom = 0, shell_top = 0;
    bool feature = false;
};

struct InterfaceUse { std::string where, a, b, orientation, kind, source; double R = 0; };

struct Model {
    Grid grid;
    std::vector<Layer> layers;
    std::vector<Instance> instances;
    std::vector<InterfaceUse> interfaces;
    std::vector<std::pair<std::string, Characterization>> sub;
    std::vector<std::pair<std::string, std::shared_ptr<Model>>> sub_models;
    std::vector<std::string> warnings;
};

struct Stack {
    std::string name = "stack";
    std::array<double, 2> cell{};
    std::vector<Layer> layers;
    double temperature = 300.0;
    BC bottom, top = BC::dirichlet(300.0);
    std::shared_ptr<MaterialLibrary> lib;
    Resolution resolution;
    Model build(const Resolution* res = nullptr) const;
};

Stack stack_from_yaml(const std::string& text);
Stack load_stack(const std::string& path);

// ---------------------------------------------------------------- runs and exports
struct Result {
    Stack stack;
    Model model;
    Characterization ch;
    std::optional<Operating> op;
    double wall_s = 0;
    std::optional<double> junction_resistance() const;
    std::string summary() const;
};
Result simulate(const Stack& stack, bool operating = true, double tol = 1e-10);

struct BudgetRow { std::string key, label; double elasticity = 0, assumed_uncertainty = 0, effect = 0; };
struct Budget { double R_total = 0, combined = 0; std::vector<BudgetRow> rows; std::string text() const; };
Budget error_budget(const Stack& stack);

struct ConvergenceRun { double factor = 1; size_t cells = 0; double R_total = 0, change = 0, change_R = 0, wall_s = 0; };
std::vector<ConvergenceRun> convergence(const Stack& stack, const std::vector<double>& factors = {1.0, 1.5, 2.0});
std::map<double, Characterization> temperature_sweep(const Stack& stack, const Vec& temperatures);

std::string layers_csv(const Characterization& ch);
std::string layers_apdl(const Characterization& ch, const std::string& name,
                        const std::map<double, Characterization>* sweep = nullptr);
std::string result_json(const Result& r, bool with_field = false);   // with_field embeds the T field
std::string report_html(const Result& r);
void write_vtr(const std::string& path, const Result& r);
std::vector<std::string> write_all(const Result& r, const std::string& dir, bool report,
                                   const std::map<double, Characterization>* sweep = nullptr);

// ---------------------------------------------------------------- atomistic glass model
// Amorphous silica with the BKS pair potential (truncated, screened Coulomb and a
// short-range repulsion).  Units here: eV, angstrom, atomic mass units.
void pair_potential(double r, int kind, double& v, double& v1, double& v2, double* v3 = nullptr);

struct Coordination { double si_four_fold = 0, o_two_fold = 0, si_o_bond = 0, si_o_bond_std = 0; };

struct Glass {
    Vec pos;                        // 3N coordinates, angstrom
    std::vector<int> species;       // 0 = Si, 1 = O
    double box = 0;                 // cubic cell edge, angstrom
    std::string info;               // provenance, a JSON object

    int n() const { return int(species.size()); }
    double mass(int i) const;
    double density() const;         // g/cm3
    void pairs(double cutoff, std::vector<int>& i, std::vector<int>& j) const;
    double energy_forces(const std::vector<int>& i, const std::vector<int>& j, const Vec& x, Vec& forces) const;
    // dense (3N x 3N) second derivatives, eV/A^2; optionally the pair vectors and blocks
    void hessian(const std::vector<int>& i, const std::vector<int>& j, const Vec& x, Vec& h,
                 Vec* d = nullptr, Vec* blocks = nullptr) const;
    // third derivatives as (row, col, slice) index triples and values, eV/A^3
    void third_order(double cutoff, std::vector<long long>& coords, Vec& data) const;
    Coordination coordination(double cutoff = 2.0) const;
};

using Log = void (*)(const std::string&);
Glass cristobalite(int cells = 3, double density = 2.20);
Glass melt_quench(const Glass& start, double t_melt = 6000.0, double t_end = 300.0, double melt_ps = 8.0,
                  double quench_ps = 60.0, double hold_ps = 2.0, double dt_fs = 1.0, unsigned seed = 1,
                  Log log = nullptr);
Glass relax(const Glass& g, Log log = nullptr);
Glass load_glass(const std::string& name);
void save_glass(const Glass& g, const std::string& name);

// ---------------------------------------------------------------- Wigner conductivity
struct Modes {
    const Glass* glass = nullptr;
    int m = 0;                      // number of modes, 3N; the first three are translations
    Vec omega;                      // rad/s
    Vec vec;                        // eigenvectors, column-major (3N x M)
    Vec vel2;                       // direction-averaged |v_ss'|^2, m^2/s^2 (M x M)
    double volume = 0;              // m^3
    double level_spacing() const { return omega.back() / (m - 3); }
};
Modes glass_modes(const Glass& g);
Vec cubic_row(const Modes& m, int s, const std::vector<int>& pi, const std::vector<int>& pj,
              double amplitude = 0.01);
// three-phonon linewidths (rad/s) of `sample` modes at each temperature: gamma[t * chosen.size() + c]
void linewidths(const Modes& m, const Vec& temperatures, int sample, double sigma_cm, std::vector<int>& chosen,
                Vec& gamma, Log log = nullptr);
Vec smooth_linewidths(const Modes& m, const std::vector<int>& chosen, const Vec& values, double width_cm = 40.0);
// Wigner conductivity (W/m/K); gamma = nullptr gives the harmonic Allen-Feldman limit
double wigner_conductivity(const Modes& m, double T, const Vec* gamma = nullptr, double eta_cm = -1.0,
                           Vec* diffusivity = nullptr);
// full pipeline for a stored glass; writes <data>/glass/<name>.wigner.json and returns a text table
std::string wigner_compute(const std::string& name, const Vec& temperatures, int sample = 160, Log log = nullptr);
std::string wigner_report(const std::string& name);

// ---------------------------------------------------------------- studio support
// The same stack with bulk conductivities, no barriers and no interface
// resistance: the conventional model every result is compared against.
Stack bulk_variant(const Stack& stack);

struct StudioRun {                  // everything the browser view needs for one solve
    std::string json;               // grid, tables, numbers
    std::vector<float> T, qx, qy, qz, kz;
    std::vector<unsigned char> feature;
    std::vector<short> mat;
    std::map<std::string, std::string> files;   // export files as text
};
StudioRun studio_run(const Stack& stack, bool budget = false, bool with_files = false);
std::string material_curves_json(double T = 300.0, double p = 0.0, double R = 0.3, double aspect = 2.0);
std::string wigner_info_json(const std::string& name = "a-SiO2-648");

// ---------------------------------------------------------------- validation
struct Check { std::string text; bool ok = true; };
Check validate_analytic();
Check validate_published();
double rayleigh_cylinders(double f, double sigma);
double cylinder_array_fv(double f, double sigma, int n, double k_matrix = 1.0, double radius = 1.0,
                         double r_int = 0.0);
double manufactured_error(int n);

}  // namespace st
