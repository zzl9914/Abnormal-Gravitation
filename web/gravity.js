/* ============================================================================
 * gravity.js —— C++ 内核 ↔ JS 的直接桥接层（无 Python、无 HTTP、无 SSE）
 *
 * 内核是 freestanding wasm32（clang -nostdlib 编译，base64 内联在
 * gravity_wasm.js 里），JS 用 WebAssembly 直接调用 gs_* 函数；
 * 状态帧是 wasm 内存上的 Float64Array 视图，**零拷贝、零 JSON**。
 *
 * 对外接口（前端只认这一层）：
 *   await Gravity.ready()            初始化（可在页面加载时调用一次）
 *   Gravity.setup(cfg)               全量重建          → {ok, error}
 *   Gravity.control(cmd)             running/dt/speed… → {ok, error}
 *   Gravity.reset()                  按最近配置重建     → {ok, error}
 *   Gravity.analysis()               派生 V(r) + 单位球 → {potR, potV, ball}
 *   Gravity.calibrate({F, r, rs})    标定 G + F 采样表  → {g, fs}
 *   Gravity.sim                      实时快照 {t,n,e,e0,bodies,error,running,G}
 *   Gravity.onFrame = fn             每渲染帧回调（替代 SSE 的 onmessage）
 *   Gravity.start() / stop()         启停 tick 循环
 * ==========================================================================*/
(function (global) {
  'use strict';

  var W = global && global.GRAVITY_WASM_B64;
  var I = null;          // wasm 导出
  var mem = null;        // WebAssembly.Memory
  var handle = 0;        // 主模拟句柄（slot 0）
  var probe = 0;         // 一次性探针句柄（slot 1）
  var lastCfg = null;
  var dec = null, enc = null;
  var kGRef = 0.04;      // 标定基准：默认 R=10、F=1/r^2、G=1 的配置

  var sim = {
    t: 0, n: 0, e: 0, e0: 0, bodies: new Float64Array(0),
    error: null, running: true, G: 1
  };

  // ---------------------------------------------------------------- 基础工具

  function b64ToBytes(b64) {
    if (typeof atob === 'function') {
      var bin = atob(b64);
      var u = new Uint8Array(bin.length);
      for (var i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
      return u;
    }
    return new Uint8Array(global.Buffer.from(b64, 'base64'));
  }

  function strPtr(s, off) {
    // 把字符串写成 UTF-8（含 '\0'）到内核的字符串暂存区，返回指针
    var cap = I.gs_str_buf_size();
    var base = I.gs_str_buf();
    var bytes = enc.encode(String(s));
    var pos = base + (off || 0);
    if (bytes.length + 1 > cap - (off || 0)) bytes = bytes.subarray(0, cap - (off || 0) - 1);
    var view = new Uint8Array(mem.buffer, pos, bytes.length + 1);
    view.set(bytes);
    view[bytes.length] = 0;
    return pos;
  }

  function doublePtr(arr, n) {
    // 把 double 数组写进内核暂存区，返回指针
    var p = I.gs_dscratch();
    var view = new Float64Array(mem.buffer, p, n);
    for (var i = 0; i < n; i++) view[i] = arr[i] || 0;
    return p;
  }

  function f64(ptr, n) {
    return new Float64Array(mem.buffer, ptr, n);
  }

  function lastError() {
    var p = I.gs_last_error();
    if (!p) return '';
    var u = new Uint8Array(mem.buffer, p, 256);
    var end = u.indexOf(0);
    return dec.decode(u.subarray(0, end < 0 ? 64 : end));
  }

  // ------------------------------------------------------------------ 初始化

  function ready() {
    if (I) return Promise.resolve(true);
    if (!W) return Promise.reject(new Error('缺少 GRAVITY_WASM_B64：先运行 build_wasm.bat'));
    dec = new TextDecoder('utf-8');
    enc = new TextEncoder();
    var imports = {
      env: {
        gs_pow: Math.pow, gs_log: Math.log, gs_exp: Math.exp,
        gs_sin: Math.sin, gs_cos: Math.cos, gs_tan: Math.tan, gs_atan: Math.atan
      }
    };
    return WebAssembly.instantiate(b64ToBytes(W), imports).then(function (r) {
      I = r.instance.exports;
      mem = I.memory;
      return true;
    });
  }

  // -------------------------------------------------------------- setup 族

  function setup(cfg) {
    if (!I) return { ok: false, error: '内核尚未初始化' };
    var dims = +(cfg.dims || 3);
    if (dims !== 2 && dims !== 3) return { ok: false, error: '维度只能是 2 或 3' };
    var fExpr = String(cfg.f || 't^2');
    var fForce = String(cfg.F || cfg.V || '1/r^2');
    var collision = +(cfg.collision || 0);
    var density = +(cfg.density !== undefined ? cfg.density : 100);
    if (!(density > 0) || !(density < 1e12)) return { ok: false, error: '密度必须为正的有限值' };
    var dt = +(cfg.dt || 0.002);
    var substeps = +(cfg.substeps || 4);
    if (!(dt > 0) || substeps < 1 || substeps > 256)
      return { ok: false, error: 'dt 必须为正，substeps 在 1..256 之间' };
    var speed = +(cfg.speed || 1.0);
    if (!(0.25 <= speed && speed <= 64)) return { ok: false, error: '速度倍率必须在 0.25 .. 64 之间' };
    var tree = +(cfg.tree || 0), theta = +(cfg.theta || 0.6);
    if (tree !== 0 && tree !== 1) return { ok: false, error: '树模式只能是 0(自动) 1(关闭)' };
    if (!(0.05 <= theta && theta <= 1.5)) return { ok: false, error: 'theta 必须在 0.05 .. 1.5 之间' };
    var bodies = cfg.bodies || [];
    if (!bodies.length) return { ok: false, error: '初始条件为空：至少需要一个质点' };

    var pF = strPtr(fExpr, 0);
    var pFrc = strPtr(fForce, Math.min(512, I.gs_str_buf_size() / 2));
    var h = I.gs_create(dims, pF, pFrc);
    if (!h) return { ok: false, error: lastError() };
    handle = h;

    I.gs_configure(h, collision, density);
    I.gs_set_tree(tree, theta);

    // 耦合常数 G：显式给定，或按半半径 autoG 标定成 G·F(r) = kGRef
    var G = 1.0;
    if (cfg.G !== undefined && cfg.G !== null) {
      G = +cfg.G;
      if (!(G > 0 && G < 1e12)) return { ok: false, error: 'G 必须为正的有限值' };
      I.gs_set_G(h, G);
    } else if (cfg.autoG !== undefined && cfg.autoG !== null && +cfg.autoG > 0) {
      var fv = I.gs_force_at(h, +cfg.autoG);
      var g = (fv !== 0 && isFinite(fv)) ? kGRef / fv : 1.0;
      if (isFinite(g) && 1e-6 <= g && g <= 1e6) { G = g; I.gs_set_G(h, g); }
    }

    for (var i = 0; i < bodies.length; i++) {
      var b = bodies[i];
      // pos 与 vel 都写进内核的 double 暂存区：先 pos 后 vel（各 dims 个）
      var base = I.gs_dscratch();
      var dv = new Float64Array(mem.buffer, base, dims * 2);
      for (var k = 0; k < dims; k++) {
        dv[k] = +(b.pos[k] || 0);
        dv[dims + k] = +(b.vel[k] || 0);
      }
      var idx = I.gs_add_body_ex(h, +b.mass, base, base + dims * 8, b.fixed ? 1 : 0);
      if (idx < 0) return { ok: false, error: lastError() || '质点参数无效' };
    }

    lastCfg = cfg;
    sim.running = cfg.running !== undefined ? !!cfg.running : true;
    sim.error = null;
    sim.G = I.gs_get_G(h);
    sim.e0 = I.gs_energy(h);
    sim.e = sim.e0;
    snapCount = 0;
    refresh(true);
    return { ok: true, error: '' };
  }

  function control(cmd) {
    if (!I || !handle) return { ok: false, error: '尚无模拟（先应用一次配置）' };
    if (cmd.running !== undefined) sim.running = !!cmd.running;
    if (lastCfg) {
      if (cmd.dt !== undefined && +cmd.dt > 0) lastCfg.dt = +cmd.dt;
      if (cmd.substeps !== undefined && cmd.substeps >= 1 && cmd.substeps <= 256)
        lastCfg.substeps = +cmd.substeps;
      if (cmd.speed !== undefined && +cmd.speed >= 0.25 && +cmd.speed <= 64)
        lastCfg.speed = +cmd.speed;
      if (cmd.tree !== undefined || cmd.theta !== undefined) {
        var tree = +(cmd.tree !== undefined ? cmd.tree : lastCfg.tree || 0);
        var theta = +(cmd.theta !== undefined ? cmd.theta : lastCfg.theta || 0.6);
        if ((tree === 0 || tree === 1) && 0.05 <= theta && theta <= 1.5) {
          lastCfg.tree = tree;
          lastCfg.theta = theta;
          I.gs_set_tree(handle, tree, theta);
        }
      }
      if (cmd.G !== undefined && +cmd.G > 0 && +cmd.G < 1e12) {
        lastCfg.G = +cmd.G;
        I.gs_set_G(handle, +cmd.G);
        sim.G = +cmd.G;
      }
    }
    if (cmd.stepOnce) I.gs_step(handle, lastCfg.dt, lastCfg.substeps);
    return { ok: true, error: '' };
  }

  function reset() {
    if (!lastCfg) return { ok: false, error: '尚无可重置的配置' };
    var cfg = {};
    for (var k in lastCfg) cfg[k] = lastCfg[k];
    return setup(cfg);
  }

  // ------------------------------------------------------------- 分析 / 标定

  function analysis() {
    if (!I || !handle) return null;
    var n = I.gs_analysis(120, 1e-2, 1e3, handle);
    var rr = Array.prototype.slice.call(f64(I.gs_pot_r_ptr(), n));
    var vv = Array.prototype.slice.call(f64(I.gs_pot_v_ptr(), n));
    var bn = I.gs_ball(96, handle);
    var ball = Array.prototype.slice.call(f64(I.gs_ball_ptr(), bn * 2));
    for (var i = 0; i < rr.length; i++) { rr[i] = round6(rr[i]); }
    for (var j = 0; j < vv.length; j++) vv[j] = isFinite(vv[j]) ? round6(vv[j]) : null;
    for (var k = 0; k < ball.length; k++) ball[k] = isFinite(ball[k]) ? round6(ball[k]) : null;
    return { potR: rr, potV: vv, ball: ball };
  }

  function round6(x) { return Math.round(x * 1e6) / 1e6; }

  // G·F(r) = 0.04 的标定 + F 采样表（任意力律的近平衡建云用）
  function calibrate(opt) {
    var F = String((opt && opt.F) || '1/r^2');
    var r = +((opt && opt.r) || 5.0);
    var rs = (opt && opt.rs) || null;
    var out = { g: 1, fs: null };
    if (!I) return out;
    var pf = strPtr(F, 0);
    var pt = strPtr('t^2', 256);
    var h = I.gs_create_probe(2, pt, pf);
    if (!h) return out;
    var fval = I.gs_force_at(h, r);
    if (rs) {
      out.fs = rs.map(function (x) {
        if (!(x > 0) || !(x < 1e9)) return null;
        var v = I.gs_force_at(h, x);
        return isFinite(v) ? v : null;
      });
    }
    var g = (fval !== 0 && isFinite(fval)) ? kGRef / fval : 1.0;
    if (!(isFinite(g) && 1e-6 <= g && g <= 1e6)) g = 1.0;
    out.g = g;
    I.gs_destroy(h);
    return out;
  }

  // ------------------------------------------------------------- tick 循环
  // 与旧服务端一致：目标 30 Hz 物理 tick（不是 60），每 tick 步数 = substeps×speed，
  // 并按墙钟预算 0.7×interval 分批下发；渲染由调用方的 rAF 驱动。

  var TICK_HZ = 30.0;
  var snapCount = 0;
  var acc = 0, lastT = 0, running = false, rafId = 0;

  function refresh(force) {
    if (!handle) return;
    var n = I.gs_count(handle);
    snapCount++;
    var cadence = n < 1000 ? 2 : 15;   // energy() 在大 N 很贵：降频
    if (force || snapCount >= cadence) {
      snapCount = 0;
      sim.e = I.gs_energy(handle);
    }
    sim.n = n;
    sim.t = I.gs_time(handle);
    sim.G = I.gs_get_G(handle);
    var cnt = I.gs_state_fill(handle);
    var buf = I.gs_state_buf();
    // 零拷贝视图：直接指向 wasm 内存（每帧新建视图，不复制数据）
    sim.bodies = new Float64Array(mem.buffer, buf, cnt * 8);
  }

  function tick(now) {
    rafId = requestAnimationFrame(tick);
    if (!lastT) lastT = now;
    var dtms = now - lastT;
    lastT = now;
    if (dtms > 250) dtms = 250;  // 标签页切回来：不要一次性补几百毫秒
    acc += dtms;
    var interval = 1000 / TICK_HZ;
    var budget = 0.7 * interval;
    var guard = 0;
    while (acc >= interval && guard++ < 4) {
      acc -= interval;
      if (sim.running && handle) {
        var sub = +(lastCfg && lastCfg.substeps) || 4;
        var speed = +(lastCfg && lastCfg.speed) || 1;
        var total = Math.max(1, Math.min(512, Math.round(sub * speed)));
        var t0 = performance.now();
        var chunk = 1, done = 0;
        while (done < total) {
          var n = Math.min(chunk, total - done);
          I.gs_step(handle, lastCfg.dt, n);
          done += n;
          if (performance.now() - t0 > budget) break;      // 预算耗尽：优雅降级
          if (performance.now() - t0 < budget * 0.25) chunk = Math.min(chunk * 4, 32, total - done) || chunk;
        }
      }
      refresh(false);
    }
    if (global.Gravity.onFrame) global.Gravity.onFrame(sim);
  }

  function start() {
    if (running) return;
    running = true;
    lastT = 0;
    rafId = requestAnimationFrame(tick);
  }
  function stop() {
    running = false;
    if (rafId) cancelAnimationFrame(rafId);
    rafId = 0;
  }

  global.Gravity = {
    ready: ready,
    setup: setup,
    control: control,
    reset: reset,
    analysis: analysis,
    calibrate: calibrate,
    start: start,
    stop: stop,
    refresh: refresh,
    sim: sim,
    onFrame: null
  };
})(typeof window !== 'undefined' ? window : globalThis);
