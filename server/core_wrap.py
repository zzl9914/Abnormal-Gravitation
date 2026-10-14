"""ctypes wrapper around gravity_core.dll.

Provides SimCore: create/destroy a simulation, add bodies, step, snapshot.
All numeric conversions happen here so server.py only deals with Python types.
"""

import ctypes
import math
import os
import threading

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_DLL = os.path.normpath(os.path.join(_HERE, "..", "bin", "gravity_core.dll"))


def _load(path: str) -> ctypes.CDLL:
    lib = ctypes.CDLL(path)
    d = ctypes.c_double
    p = ctypes.POINTER
    lib.gs_version.restype = ctypes.c_char_p
    lib.gs_last_error.restype = ctypes.c_char_p
    lib.gs_create.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p]
    lib.gs_create.restype = ctypes.c_void_p
    lib.gs_destroy.argtypes = [ctypes.c_void_p]
    lib.gs_configure.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
    lib.gs_set_tree.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
    lib.gs_set_G.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib.gs_force_at.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib.gs_force_at.restype = ctypes.c_double
    lib.gs_add_body.argtypes = [ctypes.c_void_p, ctypes.c_double, p(d), p(d)]
    lib.gs_add_body_ex.argtypes = [ctypes.c_void_p, ctypes.c_double, p(d), p(d), ctypes.c_int]
    lib.gs_clear_bodies.argtypes = [ctypes.c_void_p]
    lib.gs_step.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_int]
    lib.gs_count.argtypes = [ctypes.c_void_p]
    lib.gs_count.restype = ctypes.c_int
    lib.gs_time.argtypes = [ctypes.c_void_p]
    lib.gs_time.restype = ctypes.c_double
    lib.gs_energy.argtypes = [ctypes.c_void_p]
    lib.gs_energy.restype = ctypes.c_double
    lib.gs_dims.argtypes = [ctypes.c_void_p]
    lib.gs_dims.restype = ctypes.c_int
    lib.gs_get_state.argtypes = [ctypes.c_void_p, p(d), ctypes.c_int]
    lib.gs_sample_potential.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double,
                                        ctypes.c_double, p(d), p(d)]
    lib.gs_unit_ball.argtypes = [ctypes.c_void_p, ctypes.c_int, p(d)]
    return lib


_LIB = None
_LIB_LOCK = threading.Lock()


def get_lib(path: str = _DEFAULT_DLL) -> ctypes.CDLL:
    global _LIB
    with _LIB_LOCK:
        if _LIB is None:
            _LIB = _load(path)
        return _LIB


class SimCore:
    """One gravity simulation instance. Not thread-safe by itself:
    the server serializes calls behind a lock."""

    def __init__(self, dims: int, f_expr: str, force_expr: str,
                 collision: int = 0, density: float = 100.0):
        self.lib = get_lib()
        self._buf = None
        self._g = 1.0
        h = self.lib.gs_create(int(dims), f_expr.encode("utf-8"), force_expr.encode("utf-8"))
        if not h:
            raise ValueError(self.lib.gs_last_error().decode("utf-8"))
        self.handle = ctypes.c_void_p(h)
        self.configure(collision, density)

    def __del__(self):
        h = getattr(self, "handle", None)
        if h is not None and getattr(self, "lib", None) is not None:
            try:
                self.lib.gs_destroy(h)
            except Exception:
                pass
            self.handle = None

    def configure(self, collision: int, density: float) -> None:
        self.lib.gs_configure(self.handle, int(collision), float(density))

    def set_tree(self, tree_mode: int, theta: float) -> None:
        """tree_mode: 0 = auto (Barnes-Hut when N >= 1500), 1 = always direct."""
        self.lib.gs_set_tree(self.handle, int(tree_mode), float(theta))

    def set_G(self, g: float) -> None:
        """Force-law coupling constant (scales the force and derived V)."""
        if g > 0 and math.isfinite(g):
            self.lib.gs_set_G(self.handle, float(g))
            self._g = float(g)

    @property
    def g(self) -> float:
        return self._g

    def force_at(self, r: float) -> float:
        """Raw force-law value F(r) at metric distance r (no G factor)."""
        return self.lib.gs_force_at(self.handle, float(r))

    def add_body(self, mass: float, pos, vel, fixed: bool = False) -> None:
        """`fixed` pins the body (dark-matter halo): no motion, no merging."""
        d = ctypes.c_double
        n = max(len(pos), len(vel))
        p = (d * n)(*[float(x) for x in pos])
        v = (d * n)(*[float(x) for x in vel])
        if self.lib.gs_add_body_ex(self.handle, float(mass), p, v, 1 if fixed else 0) < 0:
            raise ValueError(self.lib.gs_last_error().decode("utf-8"))

    def clear_bodies(self) -> None:
        self.lib.gs_clear_bodies(self.handle)

    def step(self, dt: float, substeps: int) -> None:
        if self.lib.gs_step(self.handle, float(dt), int(substeps)) != 0:
            raise RuntimeError(self.lib.gs_last_error().decode("utf-8"))

    def count(self) -> int:
        return self.lib.gs_count(self.handle)

    def time(self) -> float:
        return self.lib.gs_time(self.handle)

    def energy(self) -> float:
        return self.lib.gs_energy(self.handle)

    def dims(self) -> int:
        return self.lib.gs_dims(self.handle)

    def state(self):
        """Returns list of dicts: {m, p:[x,y,z], v:[vx,vy,vz], flags} for alive bodies."""
        if self._buf is None or len(self._buf) < 8 * 4096:
            self._buf = (ctypes.c_double * (8 * 4096))()
        n = self.lib.gs_get_state(self.handle, self._buf, 4096)
        out = []
        flat = self._buf
        for i in range(n):
            base = 8 * i
            out.append({
                "m": flat[base],
                "p": [flat[base + 1], flat[base + 2], flat[base + 3]],
                "v": [flat[base + 4], flat[base + 5], flat[base + 6]],
                "flags": int(flat[base + 7]),
            })
        return out

    def flat_state(self):
        """Compact list [m,x,y,z,vx,vy,vz,flags, ...] for JSON streaming."""
        if self._buf is None or len(self._buf) < 8 * 4096:
            self._buf = (ctypes.c_double * (8 * 4096))()
        n = self.lib.gs_get_state(self.handle, self._buf, 4096)
        return list(self._buf[: 8 * n])

    def sample_potential(self, n=120, rmin=1e-2, rmax=1e3):
        d = ctypes.c_double
        rr = (d * n)()
        vv = (d * n)()
        self.lib.gs_sample_potential(self.handle, n, rmin, rmax, rr, vv)
        return list(rr), list(vv)

    def unit_ball(self, n=96):
        d = ctypes.c_double
        xy = (d * (2 * n))()
        got = self.lib.gs_unit_ball(self.handle, n, xy)
        return [xy[i] for i in range(2 * got)]

    def __repr__(self):
        return f"<SimCore t={self.time():.4g} n={self.count()}>"
