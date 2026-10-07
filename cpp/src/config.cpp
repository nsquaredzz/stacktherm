// Stack description from YAML.
#include "stacktherm/stacktherm.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

namespace st {

namespace {

// Lookup that also follows YAML merge keys ("<<: *anchor"), which yaml-cpp
// leaves unresolved.
YAML::Node find(const YAML::Node& n, const std::string& key) {
    if (!n.IsMap()) return YAML::Node(YAML::NodeType::Undefined);
    if (n[key].IsDefined()) return n[key];
    const YAML::Node merged = n["<<"];
    if (!merged.IsDefined()) return YAML::Node(YAML::NodeType::Undefined);
    if (merged.IsMap()) return find(merged, key);
    if (merged.IsSequence())
        for (const auto& m : merged) {
            YAML::Node hit = find(m, key);
            if (hit.IsDefined()) return hit;
        }
    return YAML::Node(YAML::NodeType::Undefined);
}

double quantity(const YAML::Node& v, const std::string& unit = "") { return parse_quantity(v.as<std::string>(), unit); }

std::array<double, 2> pair_of(const YAML::Node& v, bool plain = false) {
    auto one = [&](const YAML::Node& x) { return plain ? x.as<double>() : quantity(x); };
    if (v.IsSequence()) {
        if (v.size() != 2) throw std::invalid_argument("expected two values in a pair");
        return {one(v[0]), one(v[1])};
    }
    return {one(v), one(v)};
}

Resolution resolution_of(const YAML::Node& d, Resolution r) {
    if (!d.IsDefined()) return r;
    if (find(d, "n_xy").IsDefined()) r.n_xy = find(d, "n_xy").as<int>();
    if (find(d, "min_cells").IsDefined()) r.min_cells = find(d, "min_cells").as<int>();
    if (find(d, "n_z").IsDefined()) r.n_z = find(d, "n_z").as<int>();
    if (find(d, "grading").IsDefined()) r.grading = find(d, "grading").as<double>();
    if (find(d, "circle_cells").IsDefined()) r.circle_cells = find(d, "circle_cells").as<int>();
    return r;
}

BC bc_of(const YAML::Node& d, BC fallback) {
    if (!d.IsDefined()) return fallback;
    BC b;
    std::string kind = find(d, "type").IsDefined() ? find(d, "type").as<std::string>() : "adiabatic";
    if (kind == "adiabatic") b.kind = BC::Adiabatic;
    else if (kind == "dirichlet") b.kind = BC::Dirichlet;
    else if (kind == "convection") b.kind = BC::Convection;
    else if (kind == "flux") b.kind = BC::Flux;
    else throw std::invalid_argument("unknown boundary kind '" + kind + "'");
    if (find(d, "T").IsDefined()) b.T = quantity(find(d, "T"));
    if (find(d, "h").IsDefined()) b.h = find(d, "h").as<double>();
    if (find(d, "q").IsDefined()) b.q = quantity(find(d, "q"), "W/m2");
    if (find(d, "r").IsDefined()) b.r = parse_interface(find(d, "r").as<std::string>());
    return b;
}

Layer layer_of(const YAML::Node& d, int index) {
    Layer l;
    std::string kind = find(d, "type").IsDefined() ? find(d, "type").as<std::string>() : "film";
    l.name = find(d, "name").IsDefined() ? find(d, "name").as<std::string>() : kind + std::to_string(index);
    if (kind == "interface") {
        l.type = Layer::Interface;
        if (find(d, "pairs").IsDefined())
            for (const auto& kv : find(d, "pairs"))
                l.pairs.emplace_back(kv.first.as<std::string>(), parse_interface(kv.second.as<std::string>()));
        if (find(d, "tbc").IsDefined()) l.R = parse_interface(find(d, "tbc").as<std::string>());
        else if (find(d, "R").IsDefined()) l.R = parse_interface(find(d, "R").as<std::string>());
        return l;
    }
    if (kind == "substack") {
        l.type = Layer::SubStack;
        l.cell = pair_of(find(d, "cell"));
        int i = 0;
        for (const auto& c : find(d, "layers")) l.layers.push_back(layer_of(c, i++));
        if (find(d, "flip").IsDefined()) l.flip = find(d, "flip").as<bool>();
        if (find(d, "lump").IsDefined()) l.lump = find(d, "lump").as<bool>();
        if (find(d, "resolution").IsDefined()) l.resolution = resolution_of(find(d, "resolution"), Resolution());
        return l;
    }
    if (!find(d, "thickness").IsDefined()) throw std::invalid_argument("layer '" + l.name + "' needs a thickness");
    l.thickness = quantity(find(d, "thickness"));
    if (find(d, "n_z").IsDefined()) l.n_z = find(d, "n_z").as<int>();
    const YAML::Node heat = find(d, "heat");
    if (heat.IsDefined()) {
        l.heat.present = true;
        l.heat.flux = quantity(find(heat, "flux"), "W/m2");
        const YAML::Node region = find(heat, "region");
        if (region.IsDefined()) {
            if (region.size() != 4)
                throw std::invalid_argument("heat region must be [x0, x1, y0, y1] as fractions of the cell");
            l.heat.has_region = true;
            for (int i = 0; i < 4; ++i) l.heat.region[i] = region[i].as<double>();
        }
    }
    if (kind == "aniso") {
        l.type = Layer::Aniso;
        const YAML::Node k = find(d, "k");
        l.kx = find(d, "kx").IsDefined() ? find(d, "kx").as<double>() : k.as<double>();
        l.ky = find(d, "ky").IsDefined() ? find(d, "ky").as<double>() : k.as<double>();
        l.kz = find(d, "kz").IsDefined() ? find(d, "kz").as<double>() : k.as<double>();
        if (find(d, "rhoc").IsDefined()) l.rhoc = find(d, "rhoc").as<double>();
        return l;
    }
    if (find(d, "size_effect").IsDefined()) l.size_effect = find(d, "size_effect").as<bool>();
    if (!find(d, "material").IsDefined()) throw std::invalid_argument("layer '" + l.name + "' needs a material");
    l.material = find(d, "material").as<std::string>();
    if (kind == "film") { l.type = Layer::Film; return l; }
    if (kind != "lines" && kind != "posts")
        throw std::invalid_argument("layer '" + l.name + "': unknown type '" + kind + "'");
    if (!find(d, "matrix").IsDefined()) throw std::invalid_argument("layer '" + l.name + "' needs a matrix material");
    l.matrix = find(d, "matrix").as<std::string>();
    if (find(d, "inverted").IsDefined()) l.inverted = find(d, "inverted").as<bool>();
    const YAML::Node barrier = find(d, "barrier");
    if (barrier.IsDefined()) {
        l.barrier.present = true;
        l.barrier.thickness = quantity(find(barrier, "thickness"));
        if (find(barrier, "k").IsDefined()) l.barrier.k = find(barrier, "k").as<double>();
        if (find(barrier, "R_extra").IsDefined()) l.barrier.R_extra = parse_interface(find(barrier, "R_extra").as<std::string>());
        if (find(barrier, "bottom").IsDefined()) l.barrier.bottom = find(barrier, "bottom").as<bool>();
    }
    if (kind == "lines") {
        l.type = Layer::Lines;
        l.width = quantity(find(d, "width"));
        l.pitch = quantity(find(d, "pitch"));
        if (find(d, "direction").IsDefined()) l.direction = find(d, "direction").as<std::string>()[0];
        if (find(d, "offset").IsDefined()) l.offset = find(d, "offset").as<double>();
        return l;
    }
    l.type = Layer::Posts;
    if (find(d, "diameter").IsDefined()) { l.size = pair_of(find(d, "diameter")); l.circle = true; }
    else {
        l.size = pair_of(find(d, "size"));
        l.circle = find(d, "shape").IsDefined() && find(d, "shape").as<std::string>() == "circle";
    }
    if (find(d, "pitch").IsDefined()) { l.post_pitch = pair_of(find(d, "pitch")); l.has_pitch = true; }
    if (find(d, "offset").IsDefined()) l.post_offset = pair_of(find(d, "offset"), true);
    return l;
}

}  // namespace

Stack stack_from_yaml(const std::string& text) {
    YAML::Node d = YAML::Load(text);
    if (!d.IsMap()) throw std::invalid_argument("the stack file is not a YAML mapping");
    Stack s;
    s.lib = std::make_shared<MaterialLibrary>();
    if (d["materials"])
        for (const auto& kv : d["materials"]) {
            const YAML::Node& m = kv.second;
            MaterialSpec spec;
            auto opt = [&](const char* key) -> std::optional<double> {
                return find(m, key).IsDefined() ? std::optional<double>(find(m, key).as<double>()) : std::nullopt;
            };
            spec.k = opt("k"); spec.density = opt("density"); spec.cp = opt("cp"); spec.porosity = opt("porosity");
            spec.p = opt("p"); spec.R = opt("R"); spec.grain_factor = opt("grain_factor");
            if (find(m, "base").IsDefined()) spec.base = find(m, "base").as<std::string>();
            if (find(m, "source").IsDefined()) spec.source = find(m, "source").as<std::string>();
            if (find(m, "table").IsDefined()) spec.tables.push_back(find(m, "table").as<std::string>());
            if (find(m, "tables").IsDefined())
                for (const auto& t : find(m, "tables")) spec.tables.push_back(t.as<std::string>());
            s.lib->add(kv.first.as<std::string>(), spec);
        }
    if (d["interfaces"])
        for (const auto& kv : d["interfaces"]) {
            std::string key = kv.first.as<std::string>();
            if (key == "default") s.lib->default_tbr = parse_interface(kv.second.as<std::string>());
            else if (key == "fallback") s.lib->fallback = kv.second.as<std::string>();
            else {
                auto bar = key.find('|');
                if (bar == std::string::npos) throw std::invalid_argument("interface key '" + key + "' must be 'A|B'");
                auto trim = [](std::string x) { x.erase(0, x.find_first_not_of(' ')); x.erase(x.find_last_not_of(' ') + 1); return x; };
                s.lib->set_interface(trim(key.substr(0, bar)), trim(key.substr(bar + 1)),
                                     parse_interface(kv.second.as<std::string>()));
            }
        }
    if (d["metal_scattering"])
        for (const auto& kv : d["metal_scattering"]) {
            auto it = s.lib->metals.find(kv.first.as<std::string>());
            if (it == s.lib->metals.end())
                throw std::invalid_argument("metal_scattering: '" + kv.first.as<std::string>() +
                                            "' is not a size-dependent metal");
            if (kv.second["p"]) it->second.p = kv.second["p"].as<double>();
            if (kv.second["R"]) it->second.R = kv.second["R"].as<double>();
        }
    if (d["name"]) s.name = d["name"].as<std::string>();
    if (!d["cell"]) throw std::invalid_argument("the stack file needs 'cell'");
    s.cell = pair_of(d["cell"]);
    if (d["temperature"]) s.temperature = quantity(d["temperature"]);
    s.resolution = resolution_of(d["resolution"], Resolution());
    const YAML::Node boundary = d["boundary"];
    s.bottom = bc_of(boundary.IsDefined() ? boundary["bottom"] : YAML::Node(YAML::NodeType::Undefined), BC());
    s.top = bc_of(boundary.IsDefined() ? boundary["top"] : YAML::Node(YAML::NodeType::Undefined),
                  BC::dirichlet(s.temperature));
    if (!d["layers"] || !d["layers"].IsSequence()) throw std::invalid_argument("the stack file needs a list of 'layers'");
    int i = 0;
    for (const auto& l : d["layers"]) s.layers.push_back(layer_of(l, i++));
    return s;
}

Stack load_stack(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::invalid_argument("cannot open '" + path + "'");
    std::stringstream ss;
    ss << f.rdbuf();
    return stack_from_yaml(ss.str());
}

}  // namespace st
