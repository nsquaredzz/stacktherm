"""Python-side commands: ``python -m stacktherm serve`` and ``... kaldo``.

Everything else (run, validate, wigner, glass, init) is the C++ ``stacktherm``
executable.
"""
from __future__ import annotations

import argparse
import sys


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m stacktherm", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("serve", help="interactive 3D studio in the browser")
    s.add_argument("--port", type=int, default=8770)
    k = sub.add_parser("kaldo", help="glass conductivity through the kALDo reference backend")
    k.add_argument("--model", default="a-SiO2-648", help="stored glass model to run")
    k.add_argument("--import", dest="import_folder", metavar="FOLDER",
                   help="compute a table from force constants written by another code instead")
    k.add_argument("--format", default="lammps", help="force-constant format for --import")
    k.add_argument("--name", help="material name for --import")
    args = p.parse_args(argv)
    try:
        if args.command == "serve":
            from .studio import serve

            serve(args.port)
            return 0
        from . import kaldo_backend

        if args.import_folder:
            if not args.name:
                raise ValueError("--import needs --name, the material name to store the table under")
            out = kaldo_backend.import_material(args.import_folder, args.name, args.format, log=print)
            print(f"use it in a stack file with\n  materials:\n    {args.name}: {{table: {args.name}}}")
        else:
            out = kaldo_backend.compute(args.model, log=print)
        print(f"stored {out['name']}: " + ", ".join(f"{r['T']:.0f} K {r['k_wigner']:.3f}"
                                                    for r in out["temperatures"]) + " W/m/K")
        return 0
    except (ValueError, KeyError, FileNotFoundError, RuntimeError) as exc:
        message = str(exc) if isinstance(exc, OSError) or not exc.args else exc.args[0]
        print(f"error: {message}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
