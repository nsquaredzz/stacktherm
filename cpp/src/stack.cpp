// Stack description -> voxel grid.
//
// A stack is a bottom-to-top list of layers in a laterally periodic unit cell.
// Each layer is extruded through its thickness: a uniform film, parallel lines
// in a matrix, a 2D array of posts in a matrix, a pre-homogenised block, or a
// sub-stack that is characterised in its own small cell first.  Material
// conductivities are resolved per feature from the layer geometry, so a 20 nm
// line and a 2 um pad of the same metal get different values.
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>

namespace st {

Resolution Resolution::refined(double f) const {
    Resolution r = *this;
    r.n_xy = int(std::ceil(n_xy * f));
    r.min_cells = int(std::ceil(min_cells * f));
    r.n_z = int(std::ceil(n_z * f));
    r.circle_cells = int(std::ceil(circle_cells * f));
    r.grading = 1.0 + (grading - 1.0) / f;
    return r;
}

namespace {

double periodic_offset(double c, double pitch, double offset) {   // signed distance to the nearest centre
    double u = c / pitch - offset + 0.5;
    return (u - std::floor(u) - 0.5) * pitch;
}

int count(double length, double pitch, const std::string& what) {
    double n = length / pitch;
    if (std::abs(n - std::round(n)) > 1e-6 || std::round(n) < 1) {
        std::ostringstream m;
        m << what << ": pitch " << pitch << " m does not divide the unit cell " << length << " m";
        throw std::invalid_argument(m.str());
    }
    return int(std::round(n));
}

Vec edges(double length, const Vec& breaks, double h, int min_cells) {
    std::set<double> pts = {0.0, length};
    for (double b : breaks) {
        double m = std::fmod(b, length);
        if (m < 0) m += length;
        pts.insert(m);
    }
    Vec merged;
    for (double p : pts)
        if (merged.empty() || p - merged.back() > 1e-7 * length) merged.push_back(p);
    merged.back() = length;
    if (merged.size() == 2) return {0.0, length};   // nothing varies along this axis
    Vec out = {0.0};
    for (size_t s = 0; s + 1 < merged.size(); ++s) {
        double a = merged[s], b = merged[s + 1];
        int n = std::max(min_cells, int(std::ceil((b - a) / h - 1e-9)));
        for (int i = 1; i <= n; ++i) out.push_back(a + (b - a) * i / n);
    }
    return out;
}

std::string describe(const Layer& l) {   // identity of a layer definition, for the sub-stack cache
    std::ostringstream s;
    s.precision(17);
    s << l.type << '|' << l.name << '|' << l.thickness << '|' << l.material << '|' << l.matrix << '|' << l.kx
      << '|' << l.ky << '|' << l.kz << '|' << l.width << '|' << l.pitch << '|' << l.offset << '|' << l.direction
      << '|' << l.size[0] << ',' << l.size[1] << '|' << l.post_pitch[0] << ',' << l.post_pitch[1] << '|'
      << l.post_offset[0] << ',' << l.post_offset[1] << '|' << l.has_pitch << l.circle << l.barrier.present << '|'
      << l.barrier.thickness << '|' << l.barrier.k << '|' << l.barrier.R_extra << '|' << l.barrier.bottom
      << l.inverted << l.size_effect << '|' << l.n_z << '|' << l.R << '|' << l.cell[0] << ',' << l.cell[1] << '|';
    for (const auto& p : l.pairs) s << p.first << '=' << p.second << ';';
    for (const auto& c : l.layers) s << '{' << describe(c) << '}';
    return s.str();
}

struct Voxels {                      // one layer on the lateral mesh
    Vec kx, ky, kz, rhoc;
    std::vector<int> inst;
    std::vector<Instance> instances;
    double side_scale = 1.0;
};

struct Builder {
    const Stack& stack;
    const MaterialLibrary& lib;
    std::vector<Instance>& instances;
    std::vector<InterfaceUse>& used;
    std::set<std::string> seen;

    Voxels voxelise(const Layer& lay, const Vec& xe, const Vec& ye, int first_id) const {
        const int nx = int(xe.size()) - 1, ny = int(ye.size()) - 1;
        const size_t nxy = size_t(nx) * ny;
        const double lx = stack.cell[0], ly = stack.cell[1], t = lay.thickness, T = stack.temperature;
        Voxels v;
        v.kx.resize(nxy); v.ky.resize(nxy); v.kz.resize(nxy); v.rhoc.resize(nxy); v.inst.assign(nxy, first_id);
        auto fill = [&](double a, double b, double c, double rc) {
            std::fill(v.kx.begin(), v.kx.end(), a); std::fill(v.ky.begin(), v.ky.end(), b);
            std::fill(v.kz.begin(), v.kz.end(), c); std::fill(v.rhoc.begin(), v.rhoc.end(), rc);
        };
        if (lay.type == Layer::Aniso) {
            Instance in;
            in.name = in.layer = lay.name; in.material = "homogenised"; in.base = "homog:" + lay.name;
            in.kx = lay.kx; in.ky = lay.ky; in.kz = lay.kz; in.note = "sub-stack characterisation";
            v.instances.push_back(in);
            fill(lay.kx, lay.ky, lay.kz, lay.rhoc);
            return v;
        }
        if (lay.type == Layer::Film) {
            Props p = lib.props(lay.material, "film", t, 0, T, lay.size_effect);
            Instance in;
            in.name = lay.name + ":" + lay.material; in.layer = lay.name; in.material = lay.material;
            in.base = p.base; in.kx = p.kx; in.ky = p.ky; in.kz = p.kz; in.note = p.note;
            v.instances.push_back(in);
            fill(p.kx, p.ky, p.kz, p.rhoc);
            return v;
        }
        if (lay.type != Layer::Lines && lay.type != Layer::Posts)
            throw std::invalid_argument("unsupported layer type for '" + lay.name + "'");

        const double tb = lay.barrier.present ? lay.barrier.thickness : 0.0;
        const double tb_bottom = (lay.barrier.present && lay.barrier.bottom) ? tb : 0.0;
        Vec frac(nxy, 0.0);
        double w_c, h_c, drawn;
        char axial;
        auto xc = [&](int i) { return 0.5 * (xe[i] + xe[i + 1]); };
        auto yc = [&](int j) { return 0.5 * (ye[j] + ye[j + 1]); };
        if (lay.type == Layer::Lines) {
            if (lay.direction != 'x' && lay.direction != 'y')
                throw std::invalid_argument(lay.name + ": direction must be 'x' or 'y'");
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    double across = lay.direction == 'x' ? yc(j) : xc(i);
                    frac[size_t(j) * nx + i] =
                        std::abs(periodic_offset(across, lay.pitch, lay.offset)) < lay.width / 2 ? 1.0 : 0.0;
                }
            w_c = lay.width - 2 * tb; h_c = t - tb_bottom; drawn = lay.width * t; axial = lay.direction;
        } else {
            const double px = lay.has_pitch ? lay.post_pitch[0] : lx, py = lay.has_pitch ? lay.post_pitch[1] : ly;
            if (lay.circle) {
                const int sub = 8;
                const double r2 = (lay.size[0] / 2) * (lay.size[0] / 2);
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i) {
                        int inside = 0;
                        for (int b = 0; b < sub; ++b) {
                            double uy = periodic_offset(ye[j] + (b + 0.5) / sub * (ye[j + 1] - ye[j]), py, lay.post_offset[1]);
                            for (int a = 0; a < sub; ++a) {
                                double ux = periodic_offset(xe[i] + (a + 0.5) / sub * (xe[i + 1] - xe[i]), px, lay.post_offset[0]);
                                if (ux * ux + uy * uy <= r2) ++inside;
                            }
                        }
                        frac[size_t(j) * nx + i] = double(inside) / (sub * sub);
                    }
            } else {
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                        frac[size_t(j) * nx + i] =
                            (std::abs(periodic_offset(yc(j), py, lay.post_offset[1])) < lay.size[1] / 2 &&
                             std::abs(periodic_offset(xc(i), px, lay.post_offset[0])) < lay.size[0] / 2) ? 1.0 : 0.0;
            }
            w_c = lay.size[0] - 2 * tb; h_c = lay.size[1] - 2 * tb; drawn = lay.size[0] * lay.size[1]; axial = 'z';
        }
        if (w_c <= 0 || h_c <= 0) throw std::invalid_argument(lay.name + ": the barrier leaves no conductor");
        std::vector<char> feat(nxy);
        size_t n_feat = 0;
        for (size_t c = 0; c < nxy; ++c) { feat[c] = frac[c] >= 0.5; n_feat += feat[c]; }
        if (n_feat == 0 || n_feat == nxy)
            throw std::invalid_argument(lay.name + ": feature is not resolved by the lateral mesh");

        Props f = lib.props(lay.material, "wire", w_c, h_c, T, lay.size_effect);
        Props m = lib.props(lay.matrix, "film", t, 0, T, lay.size_effect);
        // along its axis the drawn feature is the conductor core in parallel with the barrier
        const double core = w_c * h_c / drawn;
        const double k_axial = f.kx * core + (lay.barrier.present ? lay.barrier.k : 0.0) * (1.0 - core);
        double fk[3] = {f.kx, f.ky, f.kz};
        fk[axial == 'x' ? 0 : axial == 'y' ? 1 : 2] = k_axial;

        if (lay.type == Layer::Posts && lay.circle) {
            // the staircase overstates the wall area; rescale the wall resistance
            double stair = 0.0;
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    if (feat[size_t(j) * nx + i] != feat[size_t(j) * nx + (i + 1) % nx]) stair += ye[j + 1] - ye[j];
                    if (feat[size_t(j) * nx + i] != feat[size_t((j + 1) % ny) * nx + i]) stair += xe[i + 1] - xe[i];
                }
            const double px = lay.has_pitch ? lay.post_pitch[0] : lx, py = lay.has_pitch ? lay.post_pitch[1] : ly;
            v.side_scale = stair / (count(lx, px, lay.name) * count(ly, py, lay.name) * kPi * lay.size[0]);
        }

        const double shell = lay.barrier.present ? lay.barrier.R() : 0.0;
        const double shell_end = (lay.barrier.present && lay.barrier.bottom) ? shell : 0.0;
        double covered = 0, total = 0;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                double a = (xe[i + 1] - xe[i]) * (ye[j + 1] - ye[j]);
                covered += frac[size_t(j) * nx + i] * a;
                total += a;
            }
        Instance fi, mi;
        fi.name = lay.name + ":" + lay.material; fi.layer = lay.name; fi.material = lay.material; fi.base = f.base;
        fi.kx = fk[0]; fi.ky = fk[1]; fi.kz = fk[2]; fi.note = f.note; fi.fraction = covered / total;
        if (tb > 0) {
            char buf[64];
            std::snprintf(buf, sizeof buf, ", core %.0fx%.0f nm", w_c * 1e9, h_c * 1e9);
            fi.note += buf;
        }
        fi.shell_side = shell; fi.shell_bottom = lay.inverted ? 0.0 : shell_end;
        fi.shell_top = lay.inverted ? shell_end : 0.0; fi.feature = true;
        mi.name = lay.name + ":" + lay.matrix; mi.layer = lay.name; mi.material = lay.matrix; mi.base = m.base;
        mi.kx = m.kx; mi.ky = m.ky; mi.kz = m.kz; mi.note = m.note; mi.fraction = 1.0 - covered / total;
        v.instances = {fi, mi};
        for (size_t c = 0; c < nxy; ++c) {
            v.kx[c] = feat[c] ? fk[0] : m.kx;
            v.ky[c] = feat[c] ? fk[1] : m.ky;
            // axial conduction scales with the exact covered area, also in cut cells
            v.kz[c] = frac[c] * fk[2] + (1.0 - frac[c]) * m.kz;
            v.rhoc[c] = frac[c] * f.rhoc + (1.0 - frac[c]) * m.rhoc;
            v.inst[c] = feat[c] ? first_id : first_id + 1;
        }
        return v;
    }

    // interface resistance between two instances; records what was applied
    double pair(int ia, int ib, const std::string& where, bool side, const std::vector<const Layer*>& entries) {
        if (ia == ib) return 0.0;
        const Instance& A = instances[ia];
        const Instance& B = instances[ib];
        const double given = lib.scale_of("user_tbr");   // bond interfaces and barriers from the stack file
        double extra = 0.0, r = 0.0;
        bool overridden = false;
        double override_value = 0.0;
        for (const Layer* e : entries) {
            extra += e->R * given;
            for (const auto& p : e->pairs) {
                auto bar = p.first.find('|');
                auto trim = [](std::string s) {
                    s.erase(0, s.find_first_not_of(' '));
                    s.erase(s.find_last_not_of(' ') + 1);
                    return s;
                };
                std::string x = lib.canonical(trim(p.first.substr(0, bar))), y = lib.canonical(trim(p.first.substr(bar + 1)));
                if ((x == A.base && y == B.base) || (x == B.base && y == A.base)) {
                    overridden = true;
                    override_value = p.second * given;
                }
            }
        }
        std::string kind, source;
        if (side) {
            Tbr t = lib.tbr(A.base, B.base, stack.temperature);
            r = t.R + (A.shell_side + B.shell_side) * given;
            kind = t.kind; source = t.source;
            if (A.shell_side + B.shell_side > 0) kind += "+barrier";
        } else if (overridden) {
            r = override_value + (A.shell_top + B.shell_bottom) * given;
            kind = "stack file";
        } else {
            Tbr t = lib.tbr(A.base, B.base, stack.temperature);
            r = t.R + extra + (A.shell_top + B.shell_bottom) * given;
            kind = t.kind; source = t.source;
            if (extra > 0) kind += "+interface";
            if (A.shell_top + B.shell_bottom > 0) kind += "+barrier";
        }
        if (r == 0.0) return 0.0;
        std::string key = where + "|" + (side && B.name < A.name ? B.name + "|" + A.name : A.name + "|" + B.name) +
                          (side ? "|s" : "|v");
        if (seen.insert(key).second)
            used.push_back({where, A.name, B.name, side ? "lateral" : "vertical", kind, source, r});
        return r;
    }
};

}  // namespace

Model Stack::build(const Resolution* res_in) const {
    if (!lib) throw std::invalid_argument("the stack has no material library");
    const Resolution res = res_in ? *res_in : resolution;
    const double lx = cell[0], ly = cell[1];
    Model model;

    // ---- sub-stacks: characterise each distinct definition once, insert its layers
    std::vector<Layer> flat;
    std::map<std::string, std::pair<std::shared_ptr<Model>, Characterization>> cache;
    for (const Layer& layer : layers) {
        if (layer.type != Layer::SubStack) { flat.push_back(layer); continue; }
        Resolution sub_res = layer.resolution ? *layer.resolution : res;
        std::ostringstream key;
        key << layer.cell[0] << ',' << layer.cell[1] << '|' << sub_res.n_xy << ',' << sub_res.min_cells << ','
            << sub_res.n_z << ',' << sub_res.grading << ',' << sub_res.circle_cells << '|';
        for (const Layer& c : layer.layers) key << '{' << describe(c) << '}';
        auto it = cache.find(key.str());
        if (it == cache.end()) {
            Stack inner;
            inner.name = layer.name; inner.cell = layer.cell; inner.layers = layer.layers;
            inner.temperature = temperature; inner.lib = lib; inner.resolution = sub_res;
            auto m = std::make_shared<Model>(inner.build());
            it = cache.emplace(key.str(), std::make_pair(m, characterize(m->grid))).first;
        }
        const Model& sm = *it->second.first;
        const Characterization& c = it->second.second;
        model.sub.emplace_back(layer.name, c);
        model.sub_models.emplace_back(layer.name, it->second.first);
        std::vector<Layer> prims;
        for (size_t n = 0; n < c.layers.size(); ++n) {
            const LayerProps& p = c.layers[n];
            Layer a;
            a.type = Layer::Aniso; a.name = layer.name + "/" + p.name; a.thickness = p.thickness;
            a.kx = p.kx; a.ky = p.ky; a.kz = p.kz; a.rhoc = p.rhoc; a.heat = sm.layers[n].heat; a.n_z = sm.layers[n].n_z;
            prims.push_back(a);
            if (p.R_above > 0) {
                Layer e;
                e.type = Layer::Interface; e.name = a.name + ":top"; e.R = p.R_above;
                prims.push_back(e);
            }
        }
        if (layer.lump) {
            double flux = 0, rc = 0;
            bool heated = false;
            for (const Layer& p : prims) if (p.type == Layer::Aniso && p.heat.present) { flux += p.heat.flux; heated = true; }
            for (const LayerProps& p : c.layers) rc += p.rhoc * p.thickness;
            Layer a;
            a.type = Layer::Aniso; a.name = layer.name; a.thickness = c.thickness;
            a.kx = c.kx; a.ky = c.ky; a.kz = c.kz; a.rhoc = rc / c.thickness;
            if (heated) { a.heat.present = true; a.heat.flux = flux; }
            prims = {a};
        }
        if (layer.flip) std::reverse(prims.begin(), prims.end());
        flat.insert(flat.end(), prims.begin(), prims.end());
    }

    // ---- split primitives from the interface entries between them
    std::vector<Layer>& lays = model.layers;
    std::vector<std::vector<Layer>> planes;
    std::vector<Layer> pending;
    for (const Layer& item : flat) {
        if (item.type == Layer::Interface) {
            if (lays.empty()) throw std::invalid_argument("an interface entry needs a layer below it");
            pending.push_back(item);
        } else {
            if (!lays.empty()) planes.push_back(pending);
            lays.push_back(item);
            pending.clear();
        }
    }
    if (!pending.empty()) throw std::invalid_argument("an interface entry needs a layer above it");
    if (lays.empty()) throw std::invalid_argument("the stack has no layers");
    const size_t nl = lays.size();

    // ---- lateral mesh, conforming to every feature edge of every layer
    Vec bx, by;
    double h = std::min(lx, ly) / res.n_xy;
    for (const Layer& lay : lays) {
        if (lay.heat.present && lay.heat.has_region) {
            bx.push_back(lay.heat.region[0] * lx); bx.push_back(lay.heat.region[1] * lx);
            by.push_back(lay.heat.region[2] * ly); by.push_back(lay.heat.region[3] * ly);
        }
        if (lay.type == Layer::Lines) {
            bool along_x = lay.direction == 'x';
            Vec& target = along_x ? by : bx;
            int n = count(along_x ? ly : lx, lay.pitch, lay.name);
            for (int m = 0; m < n; ++m) {
                double c = (m + lay.offset) * lay.pitch;
                target.push_back(c - lay.width / 2); target.push_back(c + lay.width / 2);
            }
        } else if (lay.type == Layer::Posts) {
            const double p[2] = {lay.has_pitch ? lay.post_pitch[0] : lx, lay.has_pitch ? lay.post_pitch[1] : ly};
            const double len[2] = {lx, ly};
            for (int a = 0; a < 2; ++a) {
                Vec& target = a == 0 ? bx : by;
                for (int m = 0; m < count(len[a], p[a], lay.name); ++m) {
                    double c = (m + lay.post_offset[a]) * p[a];
                    target.push_back(c - lay.size[a] / 2); target.push_back(c + lay.size[a] / 2);
                }
            }
            if (lay.circle) h = std::min(h, std::min(lay.size[0], lay.size[1]) / res.circle_cells);
        }
    }
    Grid& g = model.grid;
    g.xe = edges(lx, bx, h, res.min_cells);
    g.ye = edges(ly, by, h, res.min_cells);
    const int nx = g.nx(), ny = g.ny();
    const bool structured = nx > 1 || ny > 1;
    double h_xy = kInf;
    if (nx > 1) for (int i = 0; i < nx; ++i) h_xy = std::min(h_xy, g.dx(i));
    if (ny > 1) for (int j = 0; j < ny; ++j) h_xy = std::min(h_xy, g.dy(j));

    // ---- vertical mesh: uniform in thin layers, graded in thick ones
    std::vector<int> nz_min(nl);
    Vec base(nl);
    for (size_t n = 0; n < nl; ++n) { nz_min[n] = lays[n].n_z ? lays[n].n_z : res.n_z; base[n] = lays[n].thickness / nz_min[n]; }
    g.ze = {0.0};
    std::vector<int> first_k(nl), last_k(nl);
    for (size_t n = 0; n < nl; ++n) {
        const double cap = structured ? 2.0 * h_xy : kInf;
        double lo = std::min({base[n], n > 0 ? 1.5 * base[n - 1] : kInf, cap});
        double hi = std::min({base[n], n + 1 < nl ? 1.5 * base[n + 1] : kInf, cap});
        Vec local;
        if (lo < 0.999 * base[n] || hi < 0.999 * base[n]) local = graded_edges(lays[n].thickness, lo, hi, res.grading);
        else { local.resize(nz_min[n] + 1); for (int i = 0; i <= nz_min[n]; ++i) local[i] = lays[n].thickness * i / nz_min[n]; }
        first_k[n] = int(g.ze.size()) - 1;
        const double z0 = g.ze.back();
        for (size_t i = 1; i < local.size(); ++i) { g.ze.push_back(z0 + local[i]); g.layer_of_k.push_back(int(n)); }
        last_k[n] = int(g.ze.size()) - 2;
        g.layer_names.push_back(lays[n].name);
    }
    g.allocate();
    const size_t nxy = size_t(nx) * ny;

    Builder b{*this, *lib, model.instances, model.interfaces, {}};
    std::vector<std::vector<int>> inst_maps(nl);
    for (size_t n = 0; n < nl; ++n) {
        const Layer& lay = lays[n];
        Voxels v = b.voxelise(lay, g.xe, g.ye, int(model.instances.size()));
        model.instances.insert(model.instances.end(), v.instances.begin(), v.instances.end());
        inst_maps[n] = v.inst;
        // lateral interface resistance inside the layer
        Vec rx2(nxy, 0.0), ry2(nxy, 0.0);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                size_t c = size_t(j) * nx + i;
                rx2[c] = b.pair(v.inst[c], v.inst[size_t(j) * nx + (i + 1) % nx], lay.name, true, {}) * v.side_scale;
                ry2[c] = b.pair(v.inst[c], v.inst[size_t((j + 1) % ny) * nx + i], lay.name, true, {}) * v.side_scale;
            }
        Vec source(nxy, 0.0);
        if (lay.heat.present) {
            std::vector<char> mask(nxy, 1);
            double area = 0.0;
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    size_t c = size_t(j) * nx + i;
                    if (lay.heat.has_region) {
                        double x = 0.5 * (g.xe[i] + g.xe[i + 1]), y = 0.5 * (g.ye[j] + g.ye[j + 1]);
                        mask[c] = x >= lay.heat.region[0] * lx && x <= lay.heat.region[1] * lx &&
                                  y >= lay.heat.region[2] * ly && y <= lay.heat.region[3] * ly;
                    }
                    if (mask[c]) area += g.dx(i) * g.dy(j);
                }
            if (area <= 0) throw std::invalid_argument(lay.name + ": heat region is empty");
            for (size_t c = 0; c < nxy; ++c) if (mask[c]) source[c] = lay.heat.flux * lx * ly / (area * lay.thickness);
        }
        for (int k = first_k[n]; k <= last_k[n]; ++k)
            for (size_t c = 0; c < nxy; ++c) {
                size_t id = size_t(k) * nxy + c;
                g.kx[id] = v.kx[c]; g.ky[id] = v.ky[c]; g.kz[id] = v.kz[c]; g.rhoc[id] = v.rhoc[c];
                g.mat[id] = v.inst[c]; g.rx[id] = rx2[c]; g.ry[id] = ry2[c]; g.q[id] += source[c];
            }
    }

    // ---- resistance on the planes between layers
    for (size_t n = 0; n + 1 < nl; ++n) {
        std::vector<const Layer*> entries;
        for (const Layer& e : planes[n]) entries.push_back(&e);
        const std::string where = lays[n].name + " | " + lays[n + 1].name;
        for (size_t c = 0; c < nxy; ++c)
            g.rz[size_t(last_k[n]) * nxy + c] = b.pair(inst_maps[n][c], inst_maps[n + 1][c], where, false, entries);
    }

    // A post that lands on nothing but matrix above and below is almost always a
    // misplaced via: say so instead of silently simulating it.
    auto patterned = [&](size_t n) { return lays[n].type == Layer::Lines || lays[n].type == Layer::Posts; };
    auto is_feature = [&](size_t n, size_t c) { return model.instances[inst_maps[n][c]].feature; };
    for (size_t n = 0; n < nl; ++n) {
        if (lays[n].type != Layer::Posts) continue;
        std::vector<size_t> near;
        if (n > 0 && patterned(n - 1)) near.push_back(n - 1);
        if (n + 1 < nl && patterned(n + 1)) near.push_back(n + 1);
        if (near.empty()) continue;
        bool touches = false;
        for (size_t m : near)
            for (size_t c = 0; c < nxy && !touches; ++c) touches = is_feature(n, c) && is_feature(m, c);
        if (!touches) {
            std::string names;
            for (size_t m : near) names += (names.empty() ? "" : " or ") + lays[m].name;
            model.warnings.push_back("posts in '" + lays[n].name + "' do not touch any feature of " + names +
                                     "; check the offsets");
        }
    }
    for (const auto& sm : model.sub_models)
        for (const std::string& w : sm.second->warnings) {
            std::string msg = sm.first + ": " + w;
            if (std::find(model.warnings.begin(), model.warnings.end(), msg) == model.warnings.end())
                model.warnings.push_back(msg);
        }
    return model;
}

}  // namespace st
