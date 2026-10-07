"""kALDo as the reference transport backend for disordered solids.

kALDo (Barbalinardo, Chen, Lundgren and Donadio, J. Appl. Phys. 128, 135104
(2020)) computes thermal conductivity from harmonic and cubic force constants.
For a single disordered cell it uses the quasi-harmonic Green-Kubo method, the
linear-response counterpart of the Wigner formulation implemented in the C++
core.

Two routes lead here:

* ``from_glass`` hands kALDo the exact force constants of a stored glass model,
  computed analytically by the C++ core;
* ``from_folder`` loads force constants written by LAMMPS, a DFT code or a
  machine-learned potential, in any format kALDo reads.  This is how a material
  the core cannot model itself becomes a solver input.

Tables are written next to the native ones (``<data>/glass/<name>.kaldo.json``),
where the C++ material library picks them up in preference to the native table.
kALDo and TensorFlow are optional: ``pip install kaldo``.
"""
from __future__ import annotations

import json
import os
import tempfile
import time
from pathlib import Path

import numpy as np

from . import _core

THZ_PER_CM = 0.0299792458
FORMATS = ("numpy", "eskm", "lammps", "shengbte", "shengbte-qe", "shengbte-d3q", "hiphive",
           "tdep", "vasp-sheng", "qe-sheng", "qe-d3q", "vasp-d3q")


def data() -> Path:
    return Path(_core.data_dir()) / "glass"


def available() -> bool:
    try:
        os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
        import kaldo  # noqa: F401
        return True
    except ImportError:
        return False


def _require():
    if not available():
        raise RuntimeError("the kALDo backend is not installed: pip install kaldo")
    import logging

    logging.getLogger("kaldo").setLevel(logging.WARNING)


def from_glass(name: str = "a-SiO2-648", third_cutoff: float | None = 6.0, folder: str | None = None):
    """kALDo ForceConstants from the core's analytic derivatives of a stored glass.

    ``third_cutoff`` (angstrom) drops cubic terms of distant pairs.  None keeps
    them all, which is exact but needs about 3 GB for 648 atoms."""
    _require()
    import ase
    import sparse
    from kaldo.forceconstants import ForceConstants
    from kaldo.observables.secondorder import SecondOrder
    from kaldo.observables.thirdorder import ThirdOrder

    g = _core.glass_arrays(name, third_cutoff if third_cutoff else -1.0)
    folder = folder or tempfile.mkdtemp(prefix="stacktherm-kaldo-")
    n = len(g["species"])
    atoms = ase.Atoms(symbols=["Si" if s == 0 else "O" for s in g["species"]], positions=g["pos"],
                      cell=[g["box"]] * 3, pbc=True)
    atoms.set_masses(g["masses"])
    second = SecondOrder.from_supercell(atoms, grid_type="C", supercell=(1, 1, 1),
                                        value=np.asarray(g["hessian"]).reshape(1, n, 3, 1, n, 3), folder=folder)
    third = ThirdOrder.from_supercell(
        atoms, supercell=(1, 1, 1), grid_type="C", folder=folder,
        value=sparse.COO(np.asarray(g["third_coords"]).T, np.asarray(g["third_data"]), shape=(3 * n,) * 3))
    return ForceConstants(atoms=atoms, supercell=(1, 1, 1), folder=folder, second_order=second,
                          third_order=third), g


def from_folder(folder: str, fmt: str = "lammps", supercell=(1, 1, 1)):
    """Force constants produced elsewhere (LAMMPS, DFT, ML potentials)."""
    if fmt not in FORMATS:
        raise ValueError(f"unknown force-constant format {fmt!r}; kALDo reads: {', '.join(FORMATS)}")
    _require()
    from kaldo.forceconstants import ForceConstants

    return ForceConstants.from_folder(folder=folder, supercell=tuple(supercell), format=fmt)


def conductivity_table(forceconstants, temperatures, sigma_cm: float = 3.0, log=None) -> dict:
    """Quasi-harmonic Green-Kubo conductivity at each temperature.

    The cubic projection onto the modes is done once and reused for every
    temperature; only the occupations change."""
    _require()
    from kaldo.conductivity import Conductivity
    from kaldo.phonons import Phonons
    from kaldo.storable import LAZY_PREFIX

    atoms = forceconstants.atoms
    if not np.array_equal(forceconstants.supercell, (1, 1, 1)):
        raise ValueError("the disordered-solid route needs force constants of a single cell")
    temps = [float(t) for t in temperatures]
    ph = Phonons(forceconstants=forceconstants, temperature=temps[0], is_classic=False,
                 third_bandwidth=sigma_cm * THZ_PER_CM, broadening_shape="gauss", storage="memory",
                 folder=tempfile.mkdtemp(prefix="stacktherm-kaldo-"))
    freq = np.array(ph.frequency).ravel()  # THz
    physical = np.array(ph.physical_mode).ravel().astype(bool)
    # With in-memory storage kALDo caches every derived quantity on the object
    # and does not drop the temperature-dependent ones when the temperature
    # changes.  Clear everything except the (temperature-independent, and by far
    # most expensive) cubic projection before each temperature.
    keep = LAZY_PREFIX + "_sparse_phase_and_potential"
    rows, t0 = [], time.perf_counter()
    for T in temps:
        for attr in [a for a in vars(ph) if a.startswith(LAZY_PREFIX) and a != keep]:
            delattr(ph, attr)
        ph.temperature = T
        gamma = np.array(ph.bandwidth).ravel()  # rad/ps
        k = Conductivity(phonons=ph, method="qhgk", storage="memory").conductivity.sum(axis=0)
        rows.append(dict(T=T, k_wigner=float(np.diagonal(k).mean()), k_allen_feldman=None,
                         mean_linewidth_cm=float(gamma[physical].mean() / (2 * np.pi) / THZ_PER_CM)))
        if log:
            log(f"T {T:5.0f} K  kALDo QHGK {rows[-1]['k_wigner']:.3f} W/m/K  "
                f"({time.perf_counter() - t0:.0f} s elapsed)")
    widths = [r["mean_linewidth_cm"] for r in rows]
    if len(rows) > 1 and max(widths) - min(widths) < 1e-9 * max(widths):
        raise RuntimeError("kALDo returned identical linewidths at every temperature: its cache "
                           "was not refreshed, so the table would be wrong")
    volume = float(atoms.get_volume())
    return dict(temperatures=rows, modes=int(physical.sum()), atoms=len(atoms),
                density_g_cm3=float(atoms.get_masses().sum() * 1.66053907 / volume),
                lowest_mode_cm=float(freq[physical].min() / THZ_PER_CM),
                highest_mode_cm=float(freq.max() / THZ_PER_CM))


def compute(name: str = "a-SiO2-648", temperatures=(50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800),
            third_cutoff: float | None = 6.0, log=None) -> dict:
    """Run a stored glass model through kALDo; writes ``<name>.kaldo.json``."""
    fc, g = from_glass(name, third_cutoff)
    table = conductivity_table(fc, temperatures, log=log)
    out = dict(name=name, backend="kaldo", structure=json.loads(g["info"] or "{}"),
               third_cutoff_A=third_cutoff, **table)
    (data() / f"{name}.kaldo.json").write_text(json.dumps(out, indent=1))
    return out


def import_material(folder: str, name: str, fmt: str = "lammps", temperatures=(200, 300, 400, 600),
                    log=None) -> dict:
    """Conductivity table for externally supplied force constants, stored under
    ``name`` so that stack files can use it: ``materials: {name: {table: name}}``."""
    table = conductivity_table(from_folder(folder, fmt), temperatures, log=log)
    out = dict(name=name, backend="kaldo", source=dict(folder=str(folder), format=fmt), **table)
    (data() / f"{name}.kaldo.json").write_text(json.dumps(out, indent=1))
    return out


def load(name: str) -> dict:
    return json.loads((data() / f"{name}.kaldo.json").read_text())
