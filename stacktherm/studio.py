"""Interactive studio: parametric unit cells solved on demand for a browser.

``python -m stacktherm serve`` starts a local HTTP server.  The page asks for a
scene and a set of parameter values; this module turns them into a stack file
and hands it to the C++ core, which returns the 3D fields for rendering.  Every
scene solve is paired with a bulk-property solve of the same geometry, so the
effect of size-dependent conductivity and interface resistance is always shown.
"""
from __future__ import annotations

import base64
import json
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import yaml

from . import _core

_PAGE = Path(__file__).with_name("studio.html")
_FIELDS = ("T", "qx", "qy", "qz", "kz", "feature", "mat")
MAX_VIEW_CELLS = 700_000


def _rng(key, label, lo, hi, step, value, unit=""):
    return dict(key=key, label=label, kind="range", min=lo, max=hi, step=step, value=value, unit=unit)


def _sel(key, label, options, value):
    return dict(key=key, label=label, kind="select", options=options, value=value)


_COMMON = [_rng("T", "Temperature", 250, 450, 10, 350, "K"),
           _rng("flux", "Heat flux", 10, 1000, 10, 100, "W/cm²")]


def _through(p: dict) -> dict:
    """Heat enters at the bottom face and leaves through an isothermal top."""
    return dict(bottom=dict(type="flux", q=f"{p['flux']} W/cm2"), top=dict(type="dirichlet", T=p["T"]))


# -- scenes: each returns a stack description (the same schema as a stack file) -----
def _wiring(p: dict) -> dict:
    w = p["width"]
    pitch, t, metal, diel = 2 * w, 2 * w, p["metal"], p["dielectric"]
    feature = dict(material=metal, matrix=diel, thickness=f"{t} nm")
    if p["barrier"] > 0 and w - 2 * p["barrier"] > 2:
        feature["barrier"] = dict(thickness=f"{p['barrier']} nm")

    def lines(name, direction):
        return dict(name=name, type="lines", width=f"{w} nm", pitch=f"{pitch} nm", direction=direction, **feature)

    def via(name):
        if p["vias"] == "none":
            return dict(name=name, type="film", material=diel, thickness=f"{t} nm")
        post = dict(name=name, type="posts", size=f"{w} nm", **feature)
        if p["vias"] == "one":
            post["offset"] = 0.25
        else:
            post["pitch"] = f"{pitch} nm"
        return post

    cap = dict(type="film", material=diel, thickness=f"{t / 2} nm")
    return dict(
        name="wiring", cell=[f"{2 * pitch} nm"] * 2, temperature=p["T"],
        resolution=dict(n_xy=32, min_cells=3, n_z=5), boundary=_through(p),
        materials=dict(ULK=dict(base="SiCOH", porosity=0.3)),
        interfaces={f"{metal}|{diel}": f"{p['tbc']} MW/m2K"},
        layers=[dict(name="ILD", **cap), lines("M1", "x"), via("V1"), lines("M2", "y"), via("V2"),
                lines("M3", "x"), dict(name="cap", **cap)])


def _bond(p: dict) -> dict:
    pitch = p["pitch"] * 1000.0  # nm
    s, t_si = p["pad"] * pitch, max(2000.0, pitch)

    def wiring(name):
        return dict(name=name, type="aniso", thickness="400 nm", kx=19.0, ky=19.0, kz=p["kz_wiring"], rhoc=1.8e6)

    def pad(name, **extra):
        return dict(name=name, type="posts", material="Cu", matrix="SiO2", size=f"{s} nm",
                    thickness=f"{p['pad_t']} nm", barrier=dict(thickness="5 nm"), **extra)

    return dict(
        name="hybrid bond", cell=[f"{pitch} nm"] * 2, temperature=p["T"],
        resolution=dict(n_xy=32, min_cells=2, n_z=4), boundary=_through(p),
        layers=[dict(name="si-bottom", type="film", material="Si", thickness=f"{t_si} nm"), wiring("wiring-b"),
                pad("pad-b"),
                dict(type="interface", name="bond", pairs={"SiO2|SiO2": f"{p['tbc_ox']} MW/m2K",
                                                           "Cu|Cu": f"{p['tbc_cu']} MW/m2K"}),
                pad("pad-t", offset=[0.5 + p["shift"] * s / pitch, 0.5], inverted=True),
                wiring("wiring-t"), dict(name="si-top", type="film", material="Si", thickness=f"{t_si} nm")])


def _bump(p: dict) -> dict:
    pitch, h = p["pitch"], p["height"]  # um
    t_si = max(20.0, 0.75 * pitch)

    def wiring(name):
        return dict(name=name, type="aniso", thickness="500 nm", kx=16.0, ky=16.0, kz=p["kz_wiring"], rhoc=1.8e6)

    def post(name, frac, material):
        return dict(name=name, type="posts", material=material, matrix="underfill",
                    diameter=f"{p['dia'] * pitch} um", thickness=f"{frac * h} um")

    return dict(
        name="microbump", cell=[f"{pitch} um"] * 2, temperature=p["T"],
        resolution=dict(n_xy=32, min_cells=2, n_z=4, circle_cells=12), boundary=_through(p),
        materials=dict(underfill=dict(k=p["k_uf"], density=1700, cp=1000)),
        interfaces={"Cu|SnAg": "1000 MW/m2K"},
        layers=[dict(name="si-bottom", type="film", material="Si", thickness=f"{t_si} um"), wiring("wiring-b"),
                post("pillar-b", 0.25, "Cu"), post("solder", 0.4, "SnAg"), post("pillar-t", 0.35, "Cu"),
                wiring("wiring-t"), dict(name="si-top", type="film", material="Si", thickness=f"{t_si} um")])


def _device(p: dict) -> dict:
    cpp = 48.0  # contacted gate pitch, nm
    contacts = dict(name="contacts", type="posts", material="W", matrix="SiO2", size=f"{p['contact']} nm",
                    pitch=[f"{cpp} nm", f"{p['fin_pitch']} nm"], thickness="40 nm")
    if p["contact"] > 8:
        contacts["barrier"] = dict(thickness="2 nm")
    sink = dict(type="dirichlet", T=p["T"])
    return dict(
        name="transistor tier", cell=[f"{2 * cpp} nm", f"{2 * p['fin_pitch']} nm"], temperature=p["T"],
        resolution=dict(n_xy=32, min_cells=3, n_z=5), boundary=dict(bottom=sink, top=sink),
        interfaces={"W|SiO2": "100 MW/m2K"},
        layers=[dict(name="substrate", type="film", material="Si", thickness="60 nm", size_effect=False),
                dict(name="fins", type="lines", material="Si", matrix="SiO2", width=f"{p['fin_w']} nm",
                     pitch=f"{p['fin_pitch']} nm", direction="x", thickness=f"{p['fin_h']} nm",
                     heat=dict(flux=f"{p['flux']} W/cm2")),
                contacts, dict(name="M0", type="film", material="Cu", thickness="20 nm")])


SCENES = {
    "wiring": dict(
        title="Wiring levels", build=_wiring, zmode="true",
        blurb="Three crossed line levels joined by vias. Heat enters at the bottom and must cross "
              "the dielectric unless a via carries it.",
        params=[_sel("metal", "Metal", ["Cu", "Ru", "Co", "W"], "Cu"),
                _rng("width", "Line width", 8, 60, 1, 16, "nm"),
                _sel("vias", "Vias", [["all", "every crossing"], ["one", "one per cell"],
                                      ["none", "none"]], "one"),
                _sel("dielectric", "Dielectric", [["SiCOH", "SiCOH low-k"], ["SiO2_PECVD", "SiO2"],
                                                  ["SiO2_wigner", "SiO2, computed from atoms"],
                                                  ["ULK", "porous ULK"]], "SiCOH"),
                _rng("barrier", "Barrier thickness", 0, 4, 0.5, 2, "nm"),
                _rng("tbc", "Metal/dielectric conductance", 10, 500, 10, 50, "MW/m²K")] + _COMMON),
    "bond": dict(
        title="Hybrid bond", build=_bond, zmode="schematic",
        blurb="Cu pads in oxide, bonded face to face between two dies. The wiring on either side "
              "is a homogenised slab; heat has to spread through it to reach the pad.",
        params=[_rng("pitch", "Bond pitch", 0.4, 10, 0.1, 1, "µm"),
                _rng("pad", "Pad size / pitch", 0.3, 0.7, 0.05, 0.5),
                _rng("pad_t", "Pad thickness", 100, 1000, 50, 500, "nm"),
                _rng("shift", "Misalignment / pad size", 0, 0.5, 0.05, 0),
                _rng("tbc_ox", "Oxide bond conductance", 50, 1000, 10, 150, "MW/m²K"),
                _rng("tbc_cu", "Cu–Cu bond conductance", 100, 5000, 100, 1000, "MW/m²K"),
                _rng("kz_wiring", "Wiring kz", 0.5, 10, 0.1, 2.5, "W/m/K")] + _COMMON),
    "bump": dict(
        title="Microbump", build=_bump, zmode="schematic",
        blurb="Cu pillar, SnAg solder and Cu pillar in underfill. The joint conducts only through "
              "the bump footprint.",
        params=[_rng("pitch", "Bump pitch", 20, 150, 5, 40, "µm"),
                _rng("dia", "Diameter / pitch", 0.3, 0.7, 0.05, 0.5),
                _rng("height", "Joint height", 10, 60, 2, 20, "µm"),
                _rng("k_uf", "Underfill conductivity", 0.2, 5, 0.1, 0.5, "W/m/K"),
                _rng("kz_wiring", "Wiring kz", 0.5, 10, 0.1, 2.5, "W/m/K")] + _COMMON),
    "device": dict(
        title="Transistor tier", build=_device, zmode="true",
        blurb="Heat is generated in the Si fins and leaves both ways: down into the substrate and "
              "up through the W contacts.",
        params=[_rng("fin_w", "Fin width", 4, 20, 1, 7, "nm"),
                _rng("fin_h", "Fin height", 30, 80, 5, 50, "nm"),
                _rng("fin_pitch", "Fin pitch", 24, 60, 2, 28, "nm"),
                _rng("contact", "Contact size", 8, 22, 1, 14, "nm")] + _COMMON),
}


def scene_list() -> list[dict]:
    return [dict(id=k, title=v["title"], blurb=v["blurb"], params=v["params"], zmode=v["zmode"])
            for k, v in SCENES.items()]


def _values(scene: dict, given: dict) -> dict:
    """Parameter values, defaulted and clamped to the declared ranges."""
    out = {}
    for spec in scene["params"]:
        v = given.get(spec["key"], spec["value"])
        if spec["kind"] == "range":
            v = min(max(float(v), spec["min"]), spec["max"])
        else:
            allowed = [o[0] if isinstance(o, list) else o for o in spec["options"]]
            v = v if v in allowed else spec["value"]
        out[spec["key"]] = v
    return out


def _payload(raw: dict) -> dict:
    """Core output -> JSON-ready dict with the field buffers base64-encoded."""
    out = json.loads(raw["json"])
    if out["cells"] > MAX_VIEW_CELLS:
        raise ValueError(f"{out['cells']:,} cells is too many to view in the browser; lower "
                         "'resolution: {n_xy: ...}' or use the command line")
    out["fields"] = {k: base64.b64encode(raw[k]).decode("ascii") for k in _FIELDS}
    return out


def scene_yaml(scene_id: str, given: dict) -> tuple[str, dict]:
    scene = SCENES[scene_id]
    values = _values(scene, given)
    return yaml.safe_dump(scene["build"](values), sort_keys=False, allow_unicode=True), values


def solve_scene(scene_id: str, given: dict) -> dict:
    t0 = time.perf_counter()
    text, values = scene_yaml(scene_id, given)
    out = _payload(_core.run_yaml(text))
    r_bulk = _core.bulk_resistance(text)
    out.update(scene=scene_id, title=SCENES[scene_id]["title"], values=values, zmode=SCENES[scene_id]["zmode"],
               bulk=dict(R=r_bulk, kz=out["effective"]["thickness"] / r_bulk),
               seconds=time.perf_counter() - t0)
    return out


def run_yaml(text: str, budget: bool = False) -> dict:
    """Solve a user-supplied stack file: fields, tables and the export files."""
    raw = _core.run_yaml(text, budget, True)
    out = _payload(raw)
    out["files"] = dict(raw["files"])
    return out


def examples() -> dict[str, str]:
    folder = Path(__file__).parent.parent / "examples"
    return {p.name: p.read_text() for p in sorted(folder.glob("*.yaml"))} if folder.is_dir() else {}


# -- server ----------------------------------------------------------------------
class _Handler(BaseHTTPRequestHandler):
    lock = threading.Lock()
    cache: dict[str, bytes] = {}

    def log_message(self, *args):  # keep the terminal quiet
        pass

    def _send(self, body: bytes, ctype: str, status: int = 200):
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, status: int = 200):
        self._send(json.dumps(obj).encode(), "application/json", status)

    def do_GET(self):
        url = urlparse(self.path)
        try:
            if url.path in ("/", "/index.html"):
                self._send(_PAGE.read_bytes(), "text/html; charset=utf-8")
            elif url.path == "/api/scenes":
                self._json(scene_list())
            elif url.path == "/api/examples":
                self._json(examples())
            elif url.path == "/api/wigner":
                self._send(_core.wigner_info().encode(), "application/json")
            elif url.path == "/api/materials":
                q = {k: float(v[0]) for k, v in parse_qs(url.query).items()}
                self._send(_core.material_curves(min(max(q.get("T", 300.0), 100.0), 600.0),
                                                 min(max(q.get("p", 0.0), 0.0), 0.95),
                                                 min(max(q.get("R", 0.3), 0.0), 0.9),
                                                 min(max(q.get("ar", 2.0), 0.5), 5.0)).encode(),
                           "application/json")
            else:
                self._json(dict(error="not found"), 404)
        except Exception as exc:  # report the failure to the page instead of dropping the socket
            self._json(dict(error=str(exc)), 500)

    def do_POST(self):
        path = urlparse(self.path).path
        if path not in ("/api/solve", "/api/run"):
            return self._json(dict(error="not found"), 404)
        try:
            req = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
            if path == "/api/run":
                with self.lock:
                    body = json.dumps(run_yaml(str(req.get("yaml", "")), bool(req.get("budget")))).encode()
                return self._send(body, "application/json")
            if req.get("scene") not in SCENES:
                return self._json(dict(error="unknown scene"), 400)
            key = json.dumps(req, sort_keys=True)
            with self.lock:  # one solve at a time; repeated requests come from the cache
                if key not in self.cache:
                    if len(self.cache) > 24:
                        self.cache.pop(next(iter(self.cache)))
                    self.cache[key] = json.dumps(solve_scene(req["scene"], req.get("values") or {})).encode()
                body = self.cache[key]
            self._send(body, "application/json")
        except Exception as exc:
            self._json(dict(error=str(exc)), 500)


def serve(port: int = 8770, host: str = "127.0.0.1") -> None:
    server = ThreadingHTTPServer((host, port), _Handler)
    print(f"stacktherm studio on http://{host}:{port}  (Ctrl+C to stop)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
