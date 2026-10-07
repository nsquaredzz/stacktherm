"""The Python shell over the C++ core: studio endpoints and the kALDo bridge.
The numerics themselves are tested in cpp/tests."""
import base64
import json
from pathlib import Path

import numpy as np
import pytest

from stacktherm import _core, kaldo_backend, studio

EXAMPLES = Path(__file__).parent.parent / "examples"


def field(payload, name, dtype="<f4"):
    g = payload["grid"]
    return np.frombuffer(base64.b64decode(payload["fields"][name]), dtype=dtype).reshape(g["nz"], g["ny"], g["nx"])


@pytest.mark.parametrize("scene", sorted(studio.SCENES))
def test_every_scene_solves_and_beats_the_bulk_model(scene):
    out = studio.solve_scene(scene, {})
    g = out["grid"]
    T, feature = field(out, "T"), field(out, "feature", "u1")
    assert T.shape == (g["nz"], g["ny"], g["nx"]) and np.isfinite(T).all()
    assert feature.any() and not feature.all()  # conductors and matrix are both present
    assert out["effective"]["R"] > out["bulk"]["R"] > 0  # size effects and interfaces add resistance
    assert not out["warnings"] and out["operating"]["power"] > 0
    assert len(out["materials"]) > int(field(out, "mat", "<i2").max())


def test_scene_parameters_are_clamped_and_change_the_answer():
    base = studio.solve_scene("wiring", {})
    assert studio.solve_scene("wiring", {"width": 9999})["values"]["width"] == 60
    assert studio.solve_scene("wiring", {"metal": "Unobtainium"})["values"]["metal"] == "Cu"
    none = studio.solve_scene("wiring", {"vias": "none"})
    every = studio.solve_scene("wiring", {"vias": "all"})
    assert none["effective"]["kz"] < base["effective"]["kz"] < every["effective"]["kz"]
    cold = studio.solve_scene("wiring", {"dielectric": "SiO2_wigner", "T": 250})
    hot = studio.solve_scene("wiring", {"dielectric": "SiO2_wigner", "T": 450})
    assert hot["effective"]["kz"] > cold["effective"]["kz"]  # the atom-computed glass follows temperature
    shifted = studio.solve_scene("bond", {"shift": 0.5})
    assert shifted["effective"]["R"] > studio.solve_scene("bond", {})["effective"]["R"]


def test_user_stack_run_returns_files_and_budget():
    text = (EXAMPLES / "microbump.yaml").read_text().replace("n_xy: 40", "n_xy: 16")
    out = studio.run_yaml(text, budget=True)
    assert set(out["files"]) == {"layers.csv", "layers_apdl.mac", "report.html", "result.json"}
    assert json.loads(out["files"]["result.json"])["effective"]["R_total_m2K_W"] == pytest.approx(
        out["effective"]["R"], rel=1e-9)
    rows = out["budget"]["rows"]
    assert rows[0]["key"] == "estimated_k"  # underfill and solder dominate this stack
    assert rows == sorted(rows, key=lambda r: -r["effect"])
    assert any("FS+MS wire" in m["model"] for m in out["materials_table"])
    assert out["zmode"] == "schematic"
    with pytest.raises(Exception, match="unknown material"):
        studio.run_yaml("name: x\ncell: [1 um, 1 um]\nlayers:\n  - {type: film, material: Nope, thickness: 1 um}\n")
    assert set(studio.examples()) == {"hybrid_bond.yaml", "microbump.yaml"}


def test_material_and_wigner_endpoints():
    m = json.loads(_core.material_curves(300.0, 0.0, 0.3, 2.0))
    cu, ru = (next(x for x in m["metals"] if x["name"] == n) for n in ("Cu", "Ru"))
    assert cu["k"][0] < 0.2 * cu["bulk"] and ru["k"][0] / ru["bulk"] > cu["k"][0] / cu["bulk"]
    assert len(m["silicon"]["points"]) == 15 and m["silicon"]["in_plane"][-1] > m["silicon"]["cross_plane"][-1]
    w = json.loads(_core.wigner_info())
    room = next(r for r in w["reference"] if r["verified"] and r["T"] == 300)
    assert abs(room["library_deviation"]) < 0.15 and abs(room["deviation"]) < 0.15
    assert len(w["samples"]) == 2 and "computed from atoms" in w["library_source"]
    assert "native" in _core.wigner_report()


def test_glass_arrays_are_consistent_force_constants():
    g = _core.glass_arrays("a-SiO2-648", 3.0)
    n = len(g["species"])
    h = np.asarray(g["hessian"])
    assert n == 648 and h.shape == (3 * n, 3 * n) and np.allclose(h, h.T)
    assert np.abs(h.reshape(n, 3, n, 3).sum(axis=2)).max() < 1e-8  # acoustic sum rule
    coords = np.asarray(g["third_coords"])
    assert coords.shape[1] == 3 and coords.min() >= 0 and coords.max() < 3 * n
    assert json.loads(g["info"])["atoms"] == 648


def test_kaldo_bridge():
    with pytest.raises(ValueError, match="unknown force-constant format"):
        kaldo_backend.from_folder(".", "not-a-format")
    stored = kaldo_backend.load("a-SiO2-648")["temperatures"]
    widths = [r["mean_linewidth_cm"] for r in stored]
    # linewidths must follow the temperature: a stale kALDo cache once froze them
    assert widths == sorted(widths) and widths[-1] > 2 * widths[0]
    native = {r["T"]: r["k_wigner"] for r in json.loads(_core.wigner_info())["result"]["temperatures"]}
    k300 = next(r["k_wigner"] for r in stored if r["T"] == 300)
    assert k300 == pytest.approx(native[300.0], rel=0.05)  # two independent codes
    pytest.importorskip("kaldo")
    from kaldo.phonons import Phonons

    fc, _ = kaldo_backend.from_glass("a-SiO2-648", third_cutoff=3.0)
    ph = Phonons(forceconstants=fc, temperature=300, storage="memory", folder=fc.folder)
    freq_cm = np.array(ph.frequency).ravel() / kaldo_backend.THZ_PER_CM
    assert freq_cm[3] == pytest.approx(33.37947683, rel=1e-5)  # same harmonic problem as the core
