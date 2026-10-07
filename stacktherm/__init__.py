"""stacktherm: Python shell over the C++ core.

The numerics (solver, material models, glass model, Wigner engine) live in the
C++ library; this package adds the browser studio and the kALDo bridge.
"""
try:
    from . import _core
except ImportError as exc:  # the extension is produced by the CMake build
    raise ImportError(
        "the stacktherm C++ core is not built: run 'cmake -S . -B build && cmake --build build' "
        "in the repository root") from exc

__version__ = _core.__version__
__all__ = ["_core", "__version__"]
