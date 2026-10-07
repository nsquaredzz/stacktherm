// Python bindings: the thin surface the web studio and the kALDo bridge need.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "stacktherm/stacktherm.hpp"

namespace py = pybind11;
using namespace st;

template <class T>
static py::bytes raw(const std::vector<T>& v) {
    return py::bytes(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "stacktherm C++ core";
    m.attr("__version__") = "0.2.0";
    m.def("data_dir", &data_dir, "Directory holding the library, published data and glass models");

    m.def("run_yaml", [](const std::string& text, bool budget, bool with_files) {
        StudioRun r;
        {
            py::gil_scoped_release release;
            r = studio_run(stack_from_yaml(text), budget, with_files);
        }
        py::dict d;
        d["json"] = r.json;
        d["T"] = raw(r.T);
        d["qx"] = raw(r.qx);
        d["qy"] = raw(r.qy);
        d["qz"] = raw(r.qz);
        d["kz"] = raw(r.kz);
        d["feature"] = raw(r.feature);
        d["mat"] = raw(r.mat);
        py::dict files;
        for (const auto& kv : r.files) files[py::str(kv.first)] = kv.second;
        d["files"] = files;
        return d;
    }, py::arg("text"), py::arg("budget") = false, py::arg("with_files") = false,
          "Solve a stack file; returns tables as JSON and the fields as raw float32 / uint8 / int16 buffers");

    m.def("bulk_resistance", [](const std::string& text) {
        py::gil_scoped_release release;
        Stack s = bulk_variant(stack_from_yaml(text));
        Model model = s.build();
        return through_plane(model.grid).r_total;
    }, "Through-plane resistance of the same stack with bulk conductivities and no interface resistance");

    m.def("summary", [](const std::string& text) {
        py::gil_scoped_release release;
        return simulate(stack_from_yaml(text)).summary();
    });
    m.def("material_curves", &material_curves_json, py::arg("T") = 300.0, py::arg("p") = 0.0, py::arg("R") = 0.3,
          py::arg("aspect") = 2.0);
    m.def("wigner_info", &wigner_info_json, py::arg("name") = "a-SiO2-648");
    m.def("wigner_report", &wigner_report, py::arg("name") = "a-SiO2-648");

    // exact force constants of a stored glass, for the kALDo bridge
    m.def("glass_arrays", [](const std::string& name, double third_cutoff) {
        Glass g = load_glass(name);
        const int n = g.n();
        std::vector<int> pi, pj;
        g.pairs(10.17, pi, pj);
        Vec h;
        g.hessian(pi, pj, g.pos, h);
        std::vector<long long> coords;
        Vec data;
        g.third_order(third_cutoff, coords, data);
        py::array_t<double> pos({n, 3}), hess({3 * n, 3 * n}), masses(n), values(data.size());
        py::array_t<int> species(n);
        py::array_t<long long> idx({py::ssize_t(data.size()), py::ssize_t(3)});
        std::copy(g.pos.begin(), g.pos.end(), pos.mutable_data());
        std::copy(h.begin(), h.end(), hess.mutable_data());       // symmetric: layout does not matter
        std::copy(data.begin(), data.end(), values.mutable_data());
        std::copy(coords.begin(), coords.end(), idx.mutable_data());
        for (int i = 0; i < n; ++i) { masses.mutable_at(i) = g.mass(i); species.mutable_at(i) = g.species[i]; }
        py::dict d;
        d["pos"] = pos;
        d["species"] = species;
        d["masses"] = masses;
        d["box"] = g.box;
        d["hessian"] = hess;
        d["third_coords"] = idx;
        d["third_data"] = values;
        d["info"] = g.info;
        return d;
    }, py::arg("name") = "a-SiO2-648", py::arg("third_cutoff") = 6.0);
}
