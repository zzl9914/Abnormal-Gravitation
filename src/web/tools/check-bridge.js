'use strict';
// 桥接层自检（Node，无需浏览器）：验证 wasm 内核 + gravity.js 的端到端行为。
//   node src/web/tools/check-bridge.js
//
// 覆盖：初始化、setup（含 fixed 晕）、标定 G、control、步进、状态帧、
// analysis、reset、错误路径（非法 f / 非法 F）。

var path = require('path');

// 模拟浏览器最小环境（gravity.js 在 Node 下挂到 globalThis）
globalThis.GRAVITY_WASM_B64 = require('../../web/gravity_wasm.js');
globalThis.requestAnimationFrame = function () { return 0; };
globalThis.cancelAnimationFrame = function () {};
globalThis.performance = require('perf_hooks').performance;
globalThis.TextEncoder = require('util').TextEncoder;
globalThis.TextDecoder = require('util').TextDecoder;

require(path.join(__dirname, '..', 'gravity.js'));

var G = globalThis;
var Gravity = globalThis.Gravity;
var pass = 0, fail = 0;
function check(cond, name, extra) {
  if (cond) { pass++; console.log('[PASS] ' + name + (extra ? '  ' + extra : '')); }
  else { fail++; console.log('[FAIL] ' + name + (extra ? '  ' + extra : '')); }
}

Gravity.ready().then(function () {
  // ---- setup：中央晕（fixed）+ 一颗圆轨道质点 ----
  var cfg = {
    dims: 2, f: 't^2', F: '1/r^2', collision: 0, density: 1e9,
    dt: 1e-4, substeps: 10, speed: 1, tree: 1, theta: 0.6,
    bodies: [
      { mass: 1000, pos: [0, 0], vel: [0, 0], fixed: true },
      { mass: 1, pos: [1, 0], vel: [0, Math.sqrt(1000)], vel0: 0 }
    ],
    running: true
  };
  var r = Gravity.setup(cfg);
  check(r.ok, 'setup ok', r.error || '');
  check(Gravity.sim.n === 2, 'alive = 2', 'n=' + Gravity.sim.n);

  var e0 = Gravity.sim.e0;
  check(isFinite(e0) && e0 < 0, 'e0 有限且为负（束缚态）', 'e0=' + e0);

  // ---- 步进 2000 步 ----
  for (var i = 0; i < 200; i++) Gravity.control({ stepOnce: true });
  var e1 = Gravity.sim.e;
  var drift = Math.abs((e1 - e0) / e0);
  check(drift < 1e-6, '圆轨道能量守恒（20×10 步）', 'dE/E=' + drift.toExponential(3));

  // ---- 状态帧：stride 8、flags ----
  Gravity.refresh(true);
  var b = Gravity.sim.bodies;
  check(b.length === 16, '状态帧长度 = n×8', 'len=' + b.length);
  check(b[0] === 1000 && b[7] === 1, 'body0 = 晕（fixed flags=1）',
        'm=' + b[0] + ' flags=' + b[7]);
  check(b[8] === 1 && b[15] === 0, 'body1 = 卫星（flags=0）', 'm=' + b[8]);

  // 轨道半径仍在 1 附近（圆轨道没漂走）
  var rr = Math.hypot(b[9], b[10]);
  check(Math.abs(rr - 1) < 0.05, '轨道半径 ≈ 1', 'r=' + rr.toFixed(4));

  // ---- 标定：G·F(5) = 0.04 → 1/r^2 时 G = 0.04×25 = 1 ----
  var cal = Gravity.calibrate({ F: '1/r^2', r: 5, rs: [1, 2, 5, 10] });
  check(Math.abs(cal.g - 1.0) < 1e-9, '标定 G（F=1/r^2, r=5）= 1', 'g=' + cal.g);
  check(cal.fs && cal.fs.length === 4 && Math.abs(cal.fs[2] - 0.04) < 1e-12,
        'F 采样表正确', 'fs=' + JSON.stringify(cal.fs));
  var cal2 = Gravity.calibrate({ F: 'r', r: 5 });
  check(Math.abs(cal2.g - 0.008) < 1e-12, '标定 G（F=r, r=5）= 0.008', 'g=' + cal2.g);

  // ---- analysis ----
  var an = Gravity.analysis();
  check(an && an.potR.length === 120 && an.potV.length === 120,
        'analysis 采样 120 点', an ? 'n=' + an.potR.length : 'null');
  check(an && an.ball.length === 192, '单位球 96 点', an ? 'n=' + (an.ball.length / 2) : '');
  // f=t^2 → 单位球是圆：|ball[i]| ≈ 1
  var ok1 = true;
  for (var k = 0; k < 96; k++) {
    var x = an.ball[2 * k], y = an.ball[2 * k + 1];
    if (Math.abs(Math.hypot(x, y) - 1) > 1e-4) ok1 = false;
  }
  check(ok1, 'f=t^2 单位球是欧氏圆');

  // ---- reset ----
  var rs = Gravity.reset();
  check(rs.ok && Gravity.sim.n === 2, 'reset 后重建', rs.error || '');

  // ---- 错误路径 ----
  var bad = Gravity.setup({ dims: 2, f: 't-1', F: '1/r^2', bodies: [{ mass: 1, pos: [0, 0], vel: [0, 0] }] });
  check(!bad.ok && /f\(0\)|严格递增/.test(bad.error || ''), '拒绝 f(0)≠0', bad.error);
  var bad2 = Gravity.setup({ dims: 2, f: 't^2', F: '1/', bodies: [{ mass: 1, pos: [0, 0], vel: [0, 0] }] });
  check(!bad2.ok, '拒绝非法 F 表达式', bad2.error);

  // ---- 多体 + 树（N=2000 走 auto 树）----
  var bodies = [];
  for (var j = 0; j < 2000; j++) {
    var a = 2 * Math.PI * (j / 2000), rad = 0.2 + 8 * Math.sqrt((j % 97) / 97);
    bodies.push({ mass: 1, pos: [rad * Math.cos(a), rad * Math.sin(a)], vel: [0, 0] });
  }
  var big = Gravity.setup({
    dims: 2, f: 't^2', F: '1/r^2', collision: 0, density: 100,
    dt: 1e-4, substeps: 1, tree: 0, theta: 0.6, bodies: bodies, running: false
  });
  check(big.ok, '2000 体 setup（树模式）', big.error || '');
  if (big.ok) {
    var t0 = G.performance.now();
    for (var s = 0; s < 5; s++) Gravity.control({ stepOnce: true });
    var ms = (G.performance.now() - t0) / 5;
    check(Gravity.sim.n === 2000, '2000 体存活', 'n=' + Gravity.sim.n);
    console.log('      wasm 2000 体：' + ms.toFixed(1) + ' ms/步');
  }

  // ---- 容量上限：N = kMaxBodies = 16384（前端 #genN 滑块上限）----
  // ⚠ 这里必须建满：超限的 addBody 会返回 -1 而 setup 仍报 ok，只有数 n 才看得出来
  var many = [];
  for (var k = 0; k < 16384; k++) {
    var a2 = 2 * Math.PI * (k / 16384), rad2 = 0.2 + 8 * Math.sqrt((k % 97) / 97);
    many.push({ mass: 1, pos: [rad2 * Math.cos(a2), rad2 * Math.sin(a2)], vel: [0, 0] });
  }
  var cap = Gravity.setup({
    dims: 2, f: 't^2', F: '1/r^2', collision: 0, density: 100,
    dt: 1e-4, substeps: 1, tree: 0, theta: 0.6, bodies: many, running: false
  });
  check(cap.ok && Gravity.sim.n === 16384, '16384 体（滑块上限）全部建成',
        (cap.error || '') + ' n=' + Gravity.sim.n);
  if (cap.ok && Gravity.sim.n === 16384) {
    var t2 = G.performance.now();
    for (var s2 = 0; s2 < 3; s2++) Gravity.control({ stepOnce: true });
    var ms2 = (G.performance.now() - t2) / 3;
    console.log('      wasm 16384 体：' + ms2.toFixed(1) + ' ms/步');
  }
  // ---- 下限：N = 1 ----
  var one = Gravity.setup({
    dims: 2, f: 't^2', F: '1/r^2', collision: 0, density: 100,
    dt: 1e-4, substeps: 1, tree: 0, theta: 0.6,
    bodies: [{ mass: 1, pos: [0, 0], vel: [0, 0] }], running: false
  });
  check(one.ok && Gravity.sim.n === 1, 'N=1（滑块下限）可建立',
        (one.error || '') + ' n=' + Gravity.sim.n);

  console.log('\n' + (fail === 0 ? 'ALL CHECKS PASSED' : fail + ' CHECK(S) FAILED') +
              '  (' + pass + ' passed, ' + fail + ' failed)');
  process.exit(fail === 0 ? 0 : 1);
}).catch(function (e) {
  console.error('内核初始化失败：', e);
  process.exit(1);
});
