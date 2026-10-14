"""Gravity simulation server.

stdlib-only (http.server + threading): serves the web UI, streams simulation
frames over SSE, accepts setup/control commands over POST.

Run:  python server.py [--port 8123] [--dll path/to/gravity_core.dll]
"""

import argparse
import json
import math
import os
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from core_wrap import SimCore, get_lib

WEB_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "web"))
TICK_HZ = 30.0


class Manager:
    def __init__(self):
        self.lock = threading.Lock()
        self.core = None          # SimCore or None
        self.setup = None         # last applied setup dict
        self.running = False
        self.error = None
        self.energy0 = None
        self._stop = threading.Event()
        self._snapCount = 0       # frames since last energy evaluation
        self._energyCache = 0.0
        self.thread = threading.Thread(target=self._sim_loop, daemon=True)
        self.thread.start()

    # ------------------------------------------------------------------ loop

    def _sim_loop(self):
        interval = 1.0 / TICK_HZ
        budget = 0.7 * interval  # wall-clock budget per tick (keep API responsive)
        while not self._stop.is_set():
            t0 = time.perf_counter()
            with self.lock:
                if self.core is not None and self.running:
                    try:
                        # Steps per tick = substeps * speed (clamped).  The
                        # batch is split into small chunks and abandoned once
                        # the tick's wall-clock budget is spent, so a high
                        # multiplier on a heavy simulation degrades to "as
                        # fast as the CPU manages" instead of starving the
                        # HTTP API and the SSE stream behind this lock.
                        sub = int(self.setup.get("substeps", 4))
                        speed = float(self.setup.get("speed", 1.0))
                        total = int(round(sub * speed))
                        if total < 1:
                            total = 1
                        if total > 512:
                            total = 512
                        # adaptive chunk: start with a single step, grow it
                        # while chunks stay cheap, stop once the budget is
                        # spent -- so cheap systems take big batches and
                        # expensive ones never hog the lock for long.
                        chunk = 1
                        done = 0
                        while done < total:
                            n = chunk if chunk <= total - done else total - done
                            tc = time.perf_counter()
                            self.core.step(self.setup["dt"], n)
                            done += n
                            if time.perf_counter() - t0 > budget:
                                break
                            if time.perf_counter() - tc < budget * 0.25:
                                chunk = min(chunk * 4, 32, total - done) if total - done > 0 else chunk
                    except RuntimeError as e:
                        self.error = str(e)
                        self.running = False
            dt_used = time.perf_counter() - t0
            time.sleep(max(0.0, interval - dt_used))

    # ------------------------------------------------------------------ api

    def apply_setup(self, cfg):
        """Rebuild the simulation from a full config dict.
        Returns (ok, error)."""
        try:
            dims = int(cfg.get("dims", 3))
            if dims not in (2, 3):
                return False, "维度只能是 2 或 3"
            f_expr = str(cfg.get("f", "t^2"))
            # Force law F(r): F > 0 attractive (points at the source body).
            # Accept "F" (preferred) or legacy "V" key name for the value.
            f_force = str(cfg.get("F", cfg.get("V", "1/r^2")))
            collision = int(cfg.get("collision", 0))
            if collision not in (0, 1, 2):
                return False, "碰撞效果只能是 0(穿过) 1(并合) 2(反弹)"
            density = float(cfg.get("density", 100.0))
            if not (density > 0) or not (density < 1e12):
                return False, "密度必须为正的有限值"
            dt = float(cfg.get("dt", 0.002))
            substeps = int(cfg.get("substeps", 4))
            if not (dt > 0) or substeps < 1 or substeps > 256:
                return False, "dt 必须为正，substeps 在 1..256 之间"
            speed = float(cfg.get("speed", 1.0))
            if not (0.25 <= speed <= 64.0):
                return False, "速度倍率必须在 0.25 .. 64 之间"
            tree = int(cfg.get("tree", 0))
            if tree not in (0, 1):
                return False, "树模式只能是 0(自动) 1(关闭)"
            theta = float(cfg.get("theta", 0.6))
            if not (0.05 <= theta <= 1.5):
                return False, "theta 必须在 0.05 .. 1.5 之间"
            bodies = cfg.get("bodies", [])
            if not bodies:
                return False, "初始条件为空：至少需要一个质点"
            core = SimCore(dims, f_expr, f_force, collision, density)
            core.set_tree(tree, theta)
            # Coupling constant G: explicit "G", or auto-calibrated from
            # "autoG" = the half-radius R/2 of the random cloud so that
            # G * F(R/2) == kGRef (the default R=10, F=1/r^2 configuration)
            # -- the force at the half-radius is then the same in every
            # setup, keeping evolution timescales comparable.
            kGRef = 0.04
            if "G" in cfg:
                g = float(cfg["G"])
                if not (g > 0 and g < 1e12):
                    return False, "G 必须为正的有限值"
                core.set_G(g)
            elif "autoG" in cfg:
                r_half = float(cfg["autoG"])
                if r_half > 0:
                    fval = core.force_at(r_half)
                    if fval != 0 and math.isfinite(fval):
                        g = kGRef / fval
                        if math.isfinite(g) and 1e-6 <= g <= 1e6:
                            core.set_G(g)
            for b in bodies:
                mass = float(b["mass"])
                pos = [float(x) for x in b["pos"]]
                vel = [float(x) for x in b["vel"]]
                if len(pos) < dims or len(vel) < dims:
                    core.__del__()
                    return False, "质点的 pos/vel 维度不足"
                # fixed = dark-matter halo: pinned at its position, never
                # moves, never merges (flags ride along in the state stream)
                core.add_body(mass, pos[:dims], vel[:dims],
                              fixed=bool(b.get("fixed", False)))
        except ValueError as e:
            return False, str(e)
        except (KeyError, TypeError) as e:
            return False, f"配置字段无效: {e}"

        with self.lock:
            old = self.core
            self.core = core
            if old is not None:
                old.__del__()
            self.setup = cfg
            self.running = bool(cfg.get("running", True))
            self.error = None
            self.energy0 = core.energy()
            self._energyCache = self.energy0
            self._snapCount = 0
        return True, ""

    def control(self, cmd):
        with self.lock:
            if self.core is None:
                return False, "尚无模拟（先应用一次配置）"
            if "running" in cmd:
                self.running = bool(cmd["running"])
            if self.setup is not None:
                if "dt" in cmd:
                    dt = float(cmd["dt"])
                    if dt > 0:
                        self.setup["dt"] = dt
                if "substeps" in cmd:
                    ss = int(cmd["substeps"])
                    if 1 <= ss <= 256:
                        self.setup["substeps"] = ss
                if "speed" in cmd:
                    sp = float(cmd["speed"])
                    if 0.25 <= sp <= 64.0:
                        self.setup["speed"] = sp
                if "tree" in cmd or "theta" in cmd:
                    tree = int(cmd.get("tree", self.setup.get("tree", 0)))
                    theta = float(cmd.get("theta", self.setup.get("theta", 0.6)))
                    if tree in (0, 1) and 0.05 <= theta <= 1.5:
                        self.setup["tree"] = tree
                        self.setup["theta"] = theta
                        try:
                            self.core.set_tree(tree, theta)
                        except Exception as e:
                            return False, str(e)
                if "G" in cmd:
                    g = float(cmd["G"])
                    if g > 0 and g < 1e12:
                        self.setup["G"] = g
                        self.core.set_G(g)
            if cmd.get("stepOnce"):
                try:
                    self.core.step(self.setup["dt"], self.setup["substeps"])
                except RuntimeError as e:
                    return False, str(e)
            return True, ""

    def reset(self):
        if self.setup is None:
            return False, "尚无可重置的配置"
        cfg = dict(self.setup)
        return self.apply_setup(cfg)

    def snapshot(self):
        """dict for one SSE frame (call while NOT holding the lock is fine:
        we take it inside).  Energy evaluation rebuilds the tree and runs a
        full PE traversal when the tree is active (~30 ms at N=2500), so it
        is refreshed at a lower cadence and cached between frames."""
        with self.lock:
            if self.core is None:
                return {"t": 0.0, "n": 0, "bodies": [], "error": self.error}
            self._snapCount += 1
            n = self.core.count()
            # small systems: every 2nd frame; large tree systems: every 15th
            cadence = 2 if n < 1000 else 15
            if self._snapCount >= cadence:
                self._snapCount = 0
                self._energyCache = self.core.energy()
            return {
                "t": self.core.time(),
                "n": n,
                "e": self._energyCache,
                "e0": self.energy0,
                "G": self.core.g,
                "bodies": self.core.flat_state(),
                "running": self.running,
                "error": self.error,
            }

    def analysis(self):
        with self.lock:
            core = self.core
            if core is None:
                return None
            rr, vv = core.sample_potential(120, 1e-2, 1e3)
            ball = core.unit_ball(96)
        # round to keep the payload small
        rr = [round(x, 6) for x in rr]
        vv = [None if x != x else round(x, 6) for x in vv]  # NaN -> null
        ball = [None if x != x else round(x, 6) for x in ball]
        return {"potR": rr, "potV": vv, "ball": ball}


MANAGER = Manager()


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    # silence default stderr logging per request (keep errors only)
    def log_message(self, fmt, *args):
        pass

    # ----------------------------------------------------------------- GET

    def do_GET(self):
        path = self.path.split("?")[0]
        if path == "/" or path == "/index.html":
            self._serve_file(os.path.join(WEB_DIR, "index.html"), "text/html; charset=utf-8")
        elif path == "/three.html":
            self._serve_file(os.path.join(WEB_DIR, "three.html"), "text/html; charset=utf-8")
        elif path == "/api/stream":
            self._stream()
        elif path == "/api/state":
            self._json(MANAGER.snapshot())
        elif path == "/api/analysis":
            data = MANAGER.analysis()
            if data is None:
                self._json({"error": "尚无模拟"}, status=400)
            else:
                self._json(data)
        else:
            self._json({"error": "not found"}, status=404)

    # ----------------------------------------------------------------- POST

    def do_POST(self):
        path = self.path.split("?")[0]
        length = int(self.headers.get("Content-Length", 0) or 0)
        if length <= 0 or length > 16 * 1024 * 1024:
            self._json({"ok": False, "error": "请求体为空或过大"}, status=400)
            return
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8"))
        except (ValueError, UnicodeDecodeError):
            self._json({"ok": False, "error": "JSON 解析失败"}, status=400)
            return

        if path == "/api/setup":
            ok, err = MANAGER.apply_setup(body)
            self._json({"ok": ok, "error": err})
        elif path == "/api/control":
            ok, err = MANAGER.control(body)
            self._json({"ok": ok, "error": err})
        elif path == "/api/calibrate":
            # body: {"F": force-law expr, "r": reference (half) radius,
            #        "rs": optional list of sample radii}
            # -> {"g": G with G*F(r) == 0.04 (default R=10, F=1/r^2),
            #     "fs": [F(r) for each r in "rs"]} -- the sampled table lets
            #     the front end build near-equilibrium clouds for ANY force law
            #     (circular speed v = sqrt(G*mEnc*F(r)*r), not the 1/r^2 shortcut).
            try:
                f_force = str(body.get("F", "1/r^2"))
                r = float(body.get("r", 5.0))
                rs = body.get("rs")
                rs = [float(x) for x in rs] if isinstance(rs, list) and rs else None
                if not (0 < r < 1e9):
                    self._json({"g": 1.0})
                else:
                    probe = SimCore(2, "t^2", f_force)
                    fval = probe.force_at(r)
                    out = {}
                    if rs is not None:
                        def _f(x):
                            try:
                                v = probe.force_at(x)
                                return v if math.isfinite(v) else None
                            except Exception:
                                return None
                        out["fs"] = [_f(x) if 0 < x < 1e9 else None for x in rs]
                    probe.__del__()
                    g = 0.04 / fval if (fval != 0 and math.isfinite(fval)) else 1.0
                    if not (math.isfinite(g) and 1e-6 <= g <= 1e6):
                        g = 1.0
                    out["g"] = g
                    self._json(out)
            except ValueError as e:
                self._json({"g": 1.0, "error": str(e)})
        elif path == "/api/reset":
            ok, err = MANAGER.reset()
            self._json({"ok": ok, "error": err})
        else:
            self._json({"error": "not found"}, status=404)

    # ------------------------------------------------------------------ misc

    def _json(self, obj, status=200):
        data = json.dumps(obj, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def _serve_file(self, path, mime):
        try:
            with open(path, "rb") as f:
                data = f.read()
        except OSError:
            self._json({"error": "not found"}, status=404)
            return
        self.send_response(200)
        self.send_header("Content-Type", mime)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _stream(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream; charset=utf-8")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "keep-alive")
        self.end_headers()
        try:
            while True:
                snap = MANAGER.snapshot()
                data = json.dumps(snap, separators=(",", ":")).encode("utf-8")
                self.wfile.write(b"data: " + data + b"\n\n")
                self.wfile.flush()
                # Large systems stream at a lower cadence: the JSON payload
                # (7 doubles x N) and its lock/GIL cost scale with N, and
                # 30 fps is visually wasted on thousands of particles.
                n = snap.get("n", 0)
                hz = TICK_HZ if n < 1500 else (12.0 if n < 3000 else 8.0)
                time.sleep(1.0 / hz)
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError, OSError):
            return


def main():
    ap = argparse.ArgumentParser(description="引力模拟服务器")
    ap.add_argument("--port", type=int, default=8123)
    ap.add_argument("--dll", default=None, help="gravity_core.dll 路径（默认 ../bin/gravity_core.dll）")
    args = ap.parse_args()
    if args.dll:
        get_lib(args.dll)  # preload with custom path
    else:
        get_lib()
    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    srv.daemon_threads = True
    print(f"引力模拟服务器已启动: http://127.0.0.1:{args.port}")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
