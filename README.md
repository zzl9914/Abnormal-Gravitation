# Abnormal Gravitation · 自定义度量与力律的引力模拟

在 2D / 3D 中模拟多体引力。**距离度量与力律均可由表达式自定义**：

- 距离：d((x₁,y₁,z₁),(x₂,y₂,z₂)) = Fₘ( f(|Δx|) + f(|Δy|) + f(|Δz|) )，其中 f 定义在非负实数上、严格单增、f(0)=0，Fₘ = f⁻¹
- 力律：F(r)，r 为度量距离。**F > 0 = 吸引**（力指向施力对象），F < 0 = 排斥：  
  a_i = mⱼ·F(d)·(xⱼ−xᵢ)/|xⱼ−xᵢ|（力恒沿连线）
- 势能 V = ∫F dr 由内核自动导出（幂和项解析积分，其余数值积分），用于能量监视与核内正则化

**守恒律**（速度 Verlet，辛积分）：

- **总角动量对任意度量严格守恒**（成对中心力 → 每步总扭矩恒为 0，实测 1e-15）
- 总能量仅在 **f = t²**（欧氏各向同性）时严格守恒；各向异性度量下中心力的力幅随方向变化、力场非保守，报告的能量会漂移——这是「中心力 × 度量距离」模型的本征特性，不是积分器误差

## 快速开始

```bat
build_wasm.bat   :: 把 src/core/*.cpp 编译成 freestanding wasm32（需要 clang，支持 --target=wasm32）
start.bat        :: 启动静态服务器 → 浏览器打开 http://127.0.0.1:8123
```

浏览器里 **JS 直接调用 wasm 内核**：没有 Python、没有 HTTP 接口、没有 SSE、没有 JSON 帧——
状态帧是 wasm 内存上的 `Float64Array` 视图（零拷贝）。wasm 以 base64 内联在
`src/web/gravity_wasm.js` 里，页面加载时零 fetch。

原生侧（测试 / 基准 / 可选 DLL）用 `build.bat`（需要 g++）。

内核正确性测试：`bin\gs_test.exe`（**63 项**断言：圆轨道能量守恒 ~1e-12、各向异性度量下力的方向/幅值与**角动量严格守恒**（T12，1e-15）、V 从 F 的解析/数值反积分、单位球形状、并合动量守恒、弹性反弹、穿核正则化、三体 8 字、**暗物质晕（钉住天体）**、Barnes-Hut 树精度与积分健康等）。

性能基准：`bin\gs_bench.exe`（`--3d` 切 3D）——2D 下 N=2000 时 Barnes-Hut 树比精确求和快 ~14×；树在 N ≤ 32k 时近似 O(N log N)。

## 架构

```
src/core/    C++17 内核
  expr.*        表达式引擎：对偶数自动求导（值+导数一次求出），RPN 编译；
                asPowerSum：把表达式分解为幂和 Σa·r^p（用于从 F 精确反积分出 V）
  sim.*         度量反演（幂函数走解析快路径，一般 f 走 safeguarded Newton）、
                中心力（沿连线）、速度 Verlet 积分、碰撞三模式、
                派生势能（幂和解析 V=∫F；一般 F 走对数网格 Gauss-Legendre 积分
                + Hermite 插值表）、能量监控、
                Barnes-Hut 四叉/八叉树（θ 判据用度量距离 mrad < θ·d 表述；
                碰撞检测用「点到盒的度量距离下界」剪枝）
  bench_main.cpp  性能基准
  test_main.cpp   断言测试（63 项）
  api.cpp       C ABI（gs_*）：原生侧给 FFI 用，wasm 侧给 JS 用，同一套符号
  gsmath.h/.cpp 数学抽象层：gs_pow/gs_log/... 原生用 libm，wasm 由 JS 的
                Math.* 提供；sqrt/fabs/copysign 走编译器内建（wasm 原生指令）
  gsstr.h/gsrt_wasm.cpp  最小字符工具 / freestanding 运行时（memcpy 等）
  ⚠ 内核是 freestanding 的：无 STL 容器、无 std::string、无 new/delete
src/web/     前端（无构建、无 CDN），Hash 2D-CAM 风格 UI
  gravity_wasm.js  由 src/web/tools/build-kernel.js 生成：wasm 的 base64 内联版
  gravity.js       JS ↔ 内核桥接：零拷贝状态帧、setup/control/标定、
                   rAF tick 循环（每 tick 墙钟预算 70%，自适应分批）
  tools/build-kernel.js  clang --target=wasm32 -nostdlib → wasm → base64 JS
  tools/serve.js   极简静态服务器（只托管文件，没有业务逻辑）
  tools/check-bridge.js  桥接层自检（Node，无需浏览器）：21 项检查
  index.html  主页：可收起半透明「初始条件」设置浮层 + 右下运行控制台 dock +
              左上 HUD。距离度量 f 预设 5 个（t^0.5 / t / t^2 / t^3 / t^10）、
              力律 F 预设 5 个（r / 1/r / 1/r² / 汤川 e^(−r)(1+r)/r² / 1/r⁶−1/r¹²），
              选预设自动带入推荐 dt/子步数；预设下拉与表达式输入框双向同步——
              手改表达式后仍等于某个预设就回显该预设，否则显示「自定义」；
              3D 随机云在球体内真三维生成、
              近平衡初始化（包围质量圆速度）。**开局默认即 128 体随机云**
              （与「生成随机云」按钮同一条链路，滑块同步到 128）
  three.html  三体花活专页（/three.html）：同一引擎但锁定 f=t²、F=1/r²（2D）。
              内置经典场景库（figure-8、微扰 8 字、Lagrange 等边三角、
              Burrau 3-4-5、双星+环双星行星、恒星-行星-卫星），
              点击画布放置天体 + 拖拽定初速度，localStorage 中途存档/恢复 +
              JSON 导出/导入
  2D 视图（平移缩放）、软件渲染 3D 轨道视图（旋转缩放）、
  度量单位球叠加、V(r)=∫F 曲线、能量监视、拖尾（大 N 自适应缩短）；
  粒子显示尺寸同样随 N 分级收缩——N 越大越接近由密度导出的真实半径
ref/     实现说明：表达式引擎与符号分解、度量反演（幂快路径 / safeguarded
         Newton）、派生势（解析反积分 / 数值表）、Verlet、碰撞三模式、
         Barnes-Hut 不变量与剪枝、freestanding 约束与 wasm 对接、
         改动禁忌与验证方式 —— 面向要改代码的人，见 ref/implementation.md
```

## 内核是怎样被 JS 调用的

```
src/core/*.cpp ──clang --target=wasm32 -nostdlib──▶ gravity_core.wasm
                                                      │ base64
                                                      ▼
                                            src/web/gravity_wasm.js
                                                      │ new WebAssembly.Instance
                                                      ▼
   index.html / three.html ──▶ src/web/gravity.js ──▶ gs_create / gs_step / gs_state_fill …
                                     ▲
                                     └── 状态帧 = Float64Array 视图，指向 wasm 内存（零拷贝）
```

- 数学函数 `gs_pow / gs_log / gs_exp / gs_sin / gs_cos / gs_tan / gs_atan`
  在 wasm 里是 **import**，由 JS 用 `Math.*` 提供；`sqrt / fabs / copysign`
  用编译器内建直接落到 wasm 原生指令，不进 import 表。
- 字符串与 double 数组参数先写进内核的静态暂存区（`gs_str_buf()` /
  `gs_dscratch()`），再把指针传给 `gs_*`——wasm 侧没有 malloc。
- 固定容量上限：**天体数 16384**（`kMaxBodies`）、树节点 65536；
  前端 `#genN` 滑块就是 1 – 16384（对数刻度），正好落在范围内。
  ⚠ 这三个数（滑块上限 / 内核常量 / wasm `--initial-memory=64MB`）必须一起改。
- 物理 tick 仍按 30 Hz 推进（不是渲染帧率），每 tick 步数 = 子步数 × 倍率，
  并按墙钟预算 0.7×间隔分批下发；CPU 跟不上时优雅降级为「能跑多快跑多快」。

实测（同一台机器，2D，f=t²，F=1/r²，树 θ=0.6）：

| 场景             | 原生 -O2（g++） | wasm -O3（clang，浏览器/Node） |
| -------------- | ----------- | ------------------------ |
| N = 2000，树    | 23.4 ms/步   | **11.0 ms/步**            |
| N = 4000，树    | 52.2 ms/步   | —                        |

wasm 比原生 -O2 还快，主要来自 clang -O3 与固定容量数组（去掉了 STL 容器的间接层）。

## 物理约定

**力律 F(r)**（变量 r = 度量距离 d），F > 0 吸引、F < 0 排斥：

```
a_i = G · m_j · F(d) · (x_j − x_i) / |x_j − x_i|
```

**耦合常数 G（自动标定）**：保持 G·F(R/2) = 0.04（R = 随机云半径；基准即默认  
R=10、F=1/r²、G=1 的配置）。换力律或盘半径后重新「应用」即重新标定——  
任何配置下随机云半半径处的力都与默认配置相同，演化时标可比。  
随机云的圆速度随之按力律实际形状初始化：v = √(G·m_enc·F(r)·r)  
（对任意力律的圆轨道解；F=1/r² 时退化为 √(G·m_enc/r)）。  
`Gravity.calibrate()` 返回标定后的 G 与 F 的 64 点采样表，前端按对数插值求 F(r)。

检验：f=t²、F=1/r² 时严格退化为牛顿引力 a = mⱼ·Δ/d³。  
F = 1/r → V = ln r（平坦旋转曲线）；F = r → V = r²/2（胡克弹簧）；  
F = e^(−r)(1+r)/r² → V = −e^(−r)/r（汤川势，屏蔽库仑，r≫1 时指数衰减）；  
F = 1/r⁶−1/r¹² → V = −1/(5r⁵)+1/(11r¹¹)（Lennard-Jones）。  
常数力 F=c（V=cr）也是合法的。  
如需短程平滑，直接写进 F 表达式即可（如 F = 2r/(r²+ε²)²，对应软化势 V = −1/(r²+ε²)）。

- 积分器：速度 Verlet（辛积分），固定步长
- 动能仍是欧氏的 ½m|v|²；只有相互作用力由自定义度量/力律决定
- 派生势能：F 是幂和（如 1/r²、1/r⁶−1/r¹²、r、c/r）→ V 解析精确；  
  其余表达式 → 对数网格 4 点 Gauss-Legendre 积分 + 三次 Hermite 插值（能量监视用，精度 ~1e-9）

**暗物质晕（fixed 天体，紫色）**：随机云的中央天体是一个被钉住的天体——  
固定在原点、永不移动、不参与并合，但对所有质点施加全部引力（等价于外部  
静态势）。反弹模式下它表现为不可推动的弹性壁（按无穷大质量反射对方）；  
并合模式下与其接触的质点直接穿过（核内线性力正则仍然生效）。  
状态流中每体第 8 个分量是 flags，bit0 = 1 即钉住。

## 天体半径与碰撞

**半径由物质密度 ρ 导出**（重天体更大）：R = (m·Γ(d/2+1)/(ρ·π^(d/2)))^(1/d)（2D 圆盘、3D 球体）。

**核内线性力正则化**：当度量距离 d < Rᵢ+Rⱼ，力幅替换为 F(R)·d/R  
（在 d=R 处连续；对应二次势 W(x) = V(R) + F(R)(x²−R²)/(2R)）。  
力在核内随 d 线性衰减至 0，r→0 奇点被完全绕开（核内只求值 V(R) 与 F(R)）。  
穿过中心是光滑的；f=t² 时能量守恒。

**碰撞效果三选一**（d ≤ Rᵢ+Rⱼ 时触发）：

- **穿过**：不做碰撞处理（软核正则化，默认）
- **并合**：完全非弹性并合，动量守恒（能量图上可见损失）；**暗物质晕不参与并合**
- **反弹**：沿两心连线的完全弹性脉冲，动量与动能严格守恒（与中心力模型一致）；  
  钉住的天体按无穷大质量处理（只反射对方）

## 性能：Barnes-Hut 树与速度倍率

**多体近似（Barnes-Hut 树）**：N ≥ 1500 时自动启用（前端 dock 可关闭/调 θ）。  
开启判据用度量距离表述 `mrad < θ·d`（mrad = 质心到单元角点的最大度量距离）。  
θ 默认 0.6（快）；调小更准更慢（实测力误差均值：θ=0.4 约 0.3%、θ=0.3 约 0.2%、
θ=0.2 约 0.06%）。  
能量监视中的势能用更紧的 θ ≤ 0.25 估计（误差 < 0.2%）。

实测（2026-10-17，2D，f=t²，F=1/r²，单线程 g++ -O2；树在存活体数 ≥ 1500 时才启用，
所以 N ≤ 1000 两列是同一条代码路径）：

| N      | 精确求和     | 树 θ=0.6   | 加速比 |
| ------ | -------- | -------- | --- |
| 1000   | 37.1 ms/步 | 36.8 ms/步 | —   |
| 2000   | 145 ms/步  | 23.0 ms/步 | ~6× |
| 4000   | —         | 53.7 ms/步 | —   |
| 8000   | —         | 129 ms/步  | —   |
| 16000  | —         | 286 ms/步  | —   |

3D 树常数更大（八叉树，同一天实测）：N=2000 直接 176.5 ms/步、树 θ=0.6 99.3 ms/步
（约 1.8×）；树在 8000 体 659 ms/步、16000 体 1.61 s/步。

**速度倍率**：前端 dock 里 ×0.25 … ×32。每 tick 步数 = 子步数 × 倍率  
（封顶 512，且受每 tick 墙钟预算约束：CPU 跟不上时自动降为「能跑多快跑多快」，  
不会卡死界面）。小系统（几百体）×16 / ×32 可以真正提速。

**粒子数量级的现实预期**（本机单线程，树 θ=0.6）：交互式（≥10 步/秒）2D 约  
4–5 千体、3D 约 2 千体；2D 树在上限 16384 体处约 0.29 秒/步（浏览器 wasm 实测  
0.15 秒/步）。更大规模需要多线程 / GPU。

## 表达式语法

- 运算：`+ - * / ^`（`^` 右结合，支持负指数）
- 函数：`sqrt exp ln(log) abs sgn sin cos tan atan`，常量：`pi e`
- f 的变量是 `t`，F 的变量是 `r`
- f 必须通过校验：f(0)=0 且严格单增（幂函数 t^p 解析求逆；一般 f 数值求逆）
- F 必须在探测点（r = 0.5/1/2）处为有限值；近距奇点由核正则化处理

## 已知数值特性（实测，非 bug）

**各向异性度量下能量非守恒**：f ≠ t² 时中心力的力幅 F(d) 随方向变化，  
力场非保守（curl ≠ 0）→ 报告的能量随时间漂移（如 t=0.8 内 −2% 量级），  
同时**角动量仍守恒到机器精度**。这是模型的本征物理（等效于一个各向异性  
的「宇宙力场」），f=t² 时退化为严格守恒。

**核内刚性发散**：物质密度 ρ 很大（→ 半径很小）时，核内线性力的简谐频率  
ω = (G·m·|F(Rij)|/Rij)^(1/2)（Rij = 两天体半径之和）会超过积分器的可解析范围。  
Verlet 在 ω·dt > 2 时指数发散（能量爆炸；哪怕只走 1 步也救不回来）。  
**主页已内置自适应 dt**：dt_eff = min(用户 dt, 1/ω)——设定的 dt 是宽松时的上限；  
setup 之前先用探针估计 ω（按 F/ρ/维度/半径/质量谱缓存），刚性局面自动缩小 dt 并  
在警告条中提示，回到宽松力律即解除。把密度调小（半径变大）可提高上限。

**自引力盘的弛豫**：N 体盘通过两体散射把总角动量在质点间重新分配，  
轨道逐渐变扁（总 L、总 E 不变）。随机云默认近平衡初始化（包围质量圆速度、  
纯切向、±8% 弥散），但几十动力时标后仍会弛豫——真实天体物理。想抑制：  
加大 N / 用 Plummer 软化力 F = 1/(r²+ε²)（V = −atan(r/ε)/ε）。

**近距交会与树模式**：点状天体（ρ ≥ 1e6）+ 大步长时近距交会解析不了，  
能量会跳变；N ≥ 1500 自动启用的树有 ~0.4% 力噪声，都会加快弛豫。  
需要精度时减小 dt / θ。

## API 摘要（src/web/gravity.js）

没有 HTTP 接口了 —— 这些都是**同进程内的 JS 调用**：

- `await Gravity.ready()` — 实例化 wasm（页面启动时调用一次）
- `Gravity.setup(cfg)` — 全量重建：`f/F/dims/collision/density/dt/substeps/speed/tree/theta/bodies`；  
  `F` 为力律表达式，F>0 吸引；`G` 显式耦合常数，或 `autoG` = 半半径自动标定；  
  bodies 每项可带 `fixed:true` = 暗物质晕，钉住不动、不并合 → `{ok, error}`
- `Gravity.control(cmd)` — `running / dt / substeps / speed / tree / theta / G / stepOnce`
- `Gravity.reset()` — 按最近配置重建
- `Gravity.calibrate({F, r, rs})` → `{g, fs}`：满足 G·F(r)=0.04 的 G，及 F 的采样表
- `Gravity.probeForce(F, rs)` → F 在各半径处的值（不建主模拟；自适应 dt 的 ω 探测用）
- `Gravity.analysis()` → `{potR, potV, ball}`：派生 V(r) 采样 + 度量单位球
- `Gravity.sim` — 快照 `{t, n, e, e0, bodies, running, error, G}`，其中 `bodies` 是  
  **wasm 内存上的 Float64Array 视图**，每体 8 个数：`[m, x, y, z, vx, vy, vz, flags]`  
  （2D 补零；flags bit0 = 钉住的暗物质晕）
- `Gravity.onFrame = fn` / `Gravity.start()` — 每渲染帧回调（替代原来的 SSE）
- 参数：`speed` ∈ [0.25, 64]；`tree` 0=自动 1=关闭；`theta` ∈ [0.05, 1.5]

底层的 C ABI（`src/core/api.cpp`，原生 FFI 同样可用）：`gs_create / gs_create_probe /
gs_configure / gs_set_G / gs_set_tree / gs_add_body_ex / gs_step / gs_count /
gs_energy / gs_state_fill / gs_analysis / gs_ball / gs_force_at / gs_last_error`。
