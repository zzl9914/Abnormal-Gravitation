'use strict';
// 把 core/ 编译成 freestanding wasm32，再 base64 内联成 JS。
// 仿 Hash 2D-CAM：产物只有两个文件，浏览器里零 fetch、零依赖。
//
//   node web/tools/build-kernel.js
//
// 需要 clang 支持 --target=wasm32（MSYS2 UCRT64 自带）：
//   D:\Program Files (x86)\msys64\ucrt64\bin\clang.exe
// 可用环境变量 GRAVITY_CLANG 覆盖。

var fs = require('fs');
var path = require('path');
var child = require('child_process');

var webRoot = path.join(__dirname, '..');
var projectRoot = path.join(webRoot, '..');
var clang = process.env.GRAVITY_CLANG ||
            'D:\\Program Files (x86)\\msys64\\ucrt64\\bin\\clang.exe';

var sources = ['expr.cpp', 'sim.cpp', 'api.cpp', 'gsrt_wasm.cpp'].map(function (f) {
  return path.join(projectRoot, 'core', f);
});

var wasm = path.join(webRoot, 'gravity_core.wasm');
var out = path.join(webRoot, 'gravity_wasm.js');

// 只导出 C ABI 需要的符号；gs_pow / gs_log / ... 故意留作 import，由 JS 用 Math.* 提供。
var exportsList = [
  'gs_version', 'gs_last_error',
  'gs_create', 'gs_create_probe', 'gs_destroy',
  'gs_configure', 'gs_set_G', 'gs_get_G', 'gs_force_at', 'gs_set_tree',
  'gs_add_body', 'gs_add_body_ex', 'gs_clear_bodies',
  'gs_step', 'gs_count', 'gs_time', 'gs_energy', 'gs_dims',
  'gs_get_state', 'gs_sample_potential', 'gs_unit_ball',
  'gs_state_buf', 'gs_state_fill',
  'gs_analysis', 'gs_pot_r_ptr', 'gs_pot_v_ptr',
  'gs_ball', 'gs_ball_ptr',
  'gs_str_buf', 'gs_str_buf_size', 'gs_dscratch', 'gs_dscratch_size',
];

var args = [
  '--target=wasm32',
  '-O3',
  '-nostdlib',
  '-ffreestanding',
  '-fno-exceptions',
  '-fno-rtti',
  '-Wl,--no-entry',
  '-Wl,--allow-undefined',
  '-Wl,--export-memory',
  // 静态缓冲：2 个 Simulation 槽位（天体表 216B/体 ×16384）+ 树节点（264B ×65536）
  // + 状态帧（64B/体），实测约 26MB → 给 64MB（改 kMaxBodies/kMaxNodes 必须同步改这里）
  '-Wl,--initial-memory=67108864',
  '-Wl,--max-memory=268435456',
  // asPowerSum 的符号栈是栈上数组（~34KB），递归插入也要余量
  '-Wl,-zstack-size=1048576',
].concat(exportsList.map(function (n) { return '-Wl,--export=' + n; }))
  .concat(['-o', wasm])
  .concat(sources);

var result = child.spawnSync(clang, args, { stdio: 'inherit' });
if (result.status !== 0) process.exit(result.status || 1);

// 记录 wasm 需要 JS 提供的 import，构建时校验一次，避免运行时才发现缺函数
var mod = new WebAssembly.Module(fs.readFileSync(wasm));
var imports = WebAssembly.Module.imports(mod);
var names = imports.map(function (i) { return i.module + '.' + i.name; });
var expected = ['env.gs_pow', 'env.gs_log', 'env.gs_exp', 'env.gs_sin',
                'env.gs_cos', 'env.gs_tan', 'env.gs_atan'];
for (var i = 0; i < expected.length; i++) {
  if (names.indexOf(expected[i]) < 0) {
    console.error('内核缺少预期的 import：' + expected[i] + '（实际：' + names.join(', ') + '）');
    process.exit(1);
  }
}

var bytes = fs.readFileSync(wasm);
var b64 = bytes.toString('base64');
var js = [
  '// 由 web/tools/build-kernel.js 从 core/*.cpp 生成，不要手改。',
  '// freestanding wasm32（clang -nostdlib），base64 内联，浏览器零 fetch。',
  'var GRAVITY_WASM_B64 = ' + JSON.stringify(b64) + ';',
  "if (typeof module === 'object' && module.exports) module.exports = GRAVITY_WASM_B64;",
  '',
].join('\n');
fs.writeFileSync(out, js);
console.log('kernel wasm ' + bytes.length + ' bytes -> ' + out +
            '  (imports: ' + names.join(', ') + ')');
