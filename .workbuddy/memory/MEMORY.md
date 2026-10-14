# MEMORY.md — Simulation of Gravity（自定义度量与势能的引力模拟）

> 2026-10-14 建立。本项目规则与本层规则同效：动手前先读本文件。

## 项目结构

```
core/   C++17 内核：expr.*（对偶数表达式引擎）sim.*（模拟器，含 Barnes-Hut 树）
        api.cpp（C ABI）test_main.cpp bench_main.cpp（性能基准）
server/ core_wrap.py（ctypes）server.py（stdlib HTTP+SSE，无三方依赖）
        ~~verify_*.py / smoke_test.py 等验证脚本已于 2026-10-15 删除~~
web/    index.html（主页，Hash 2D-CAM 风格）+ three.html（三体花活页，/three.html），
        均无构建无 CDN
bin/    gravity_core.dll（-static 全静态）、gs_test.exe、gs_bench.exe（构建产物，已 gitignore）
build.bat / build.sh / run.bat（端口 8123）
.gitignore / .gitattributes（2026-10-15 新增）
GitHub：https://github.com/zzl9914/Simulation-of-Gravity （public，master，
        2026-10-15 首次推送；⚠ 未加 LICENSE）
```

## UI 约定（第 8 轮定版）

- 配色沿用 Hash 2D-CAM：bg #0c100c / panel #141a14 / line #2a382c / text #e7f0e4 /
  dim #93a394 / acc #39ff14 霓虹绿 / ink #042000。
- 布局：全屏 canvas + 左上 .hud + 右上 settings-toggle（body.is-settings 时画布 opacity .34）
  + 右下 .dock 运行控制台 + 左下 .corner-left（fit 圆钮/图例/可收起监视器）。
- 模拟对象完全拆分（不再捆绑场景预设）：F_PRESETS 5（t^0.5/t/t^2/t^3/t^10）、
  **力律 F 预设 4**（r/1/r/1/r²/1/r⁶−1/r¹²，F>0=吸引），预设带推荐 dt/sub，
  选择即 setDtSub+applySetup。前端 id 沿用 #vPreset/#vExpr（内容已是力律）。
- 3D 随机云为球体均匀采样：`rr=(0.25+0.75·cbrt(u))·R` + 球面方向 z=2u-1。
- makeCloud 为**近平衡初始化**：圆速度计入包围质量（中央体 + 均匀盘∝r² / 均匀球∝r³）、
  纯切向、±8% 速度弥散、默认旋转系数 1.0；弥散 3D 为 0.3·v0。初始 ⟨e⟩≈0.13。
- three.html：锁定 f=t²/V=−1/r/2D。经典场景库 6 个（figure-8、微扰 8 字、Lagrange、
  Burrau 3-4-5、双星+环双星行星、恒星-行星-卫星）；画布点按放天体+拖拽定初速；
  localStorage 存档（存**实时状态**，键 grav3body.scenes / .last）+ JSON 导出/导入。
- boot 必须 sim.running=true 后再 applySetup（曾只改按钮文字导致启动即暂停）。

## 物理与实现约定（改动前必读）

- **模型是「力律 F + 中心力」**（第 9 轮定版）：用户给 F(r)，**F>0=吸引（指向施力对象）**，
  F<0=排斥；力 a_i = G·m_j·F(d)·Δ/|Δ|（Δ=x_j−x_i，沿连线）。
  **V = ∫F dr 由内核导出**（F = V′！不是 −V′）：幂和走 `Expr::asPowerSum` 解析反积分
  （RPN 符号机，p=−1 项 → ln r）；其余表达式走数值表（[1e-6,1e6] 4096 点 GL4 + Hermite）。
- **耦合常数 G（第 10 轮）**：`setG(G)`（G>0 有限才生效），力与派生 V 同步乘 G；
  **自动标定 G·F(R/2)=0.04**（R=随机云半径，基准 = 默认 R=10、F=1/r²、G=1 配置，
  换力律/盘半径后演化时标可比）。server setup 收 `"G"` 显式或 `"autoG"`=半半径
  （钳 1e-6..1e6）；`POST /api/calibrate {F,r,rs?}→{g,fs?}` 返回 G 与 F 采样表。
  ⚠ **建云圆速度按力律实际形状**：v=√(G·m_enc·F(r)·r)（对数插值采样表；F=1/r² 退化
  √(G·m_enc/r)；F≤0 处兜底 0.6√G；无表退回 1/r² 假设）——曾按 √(G·mEnc/r) 隐含
  1/r² 导致 LJ 云 6 倍膨胀。
- **角动量对任意度量守恒到舍入**（中心力对瞬时总扭矩恒 0，实测 1e-15）；
  **能量仅在 f=t² 严格守恒**——各向异性度量下中心力幅随方向变、场非保守，E 漂移是
  本征物理（不是 bug，T12 演示：dL/L=1e-15 同时 dE/E=−2%）。
- **力循环不再需要度量梯度/对偶数于力的方向**（方向恒 Δ̂）；对偶数仍用于度量反演的
  Newton 与 f 的梯度（pairKernel，测试用）。f=t^p 解析快路径不变。
- 一般 f 的校验：f(0)=0 + 对数网格上**严格** `v > prev`；F 的校验：r=0.5/1/2 探测点有限。
- 并合默认**关闭**；collision ∈ {0=穿过,1=并合,2=反弹}；反弹法向 = Δ̂（f=t² 时与旧一致）。
- **核内线性力正则**（d < Rᵢ+Rⱼ）：s = F(R)·d/R（对应 W(x)=V(R)+F(R)(x²−R²)/(2R)，
  对 f=t² 与外力在边界连续）。密度大 → 半径小 → ω~√(F(R)/R)，ω·dt>2 必炸。
- 半径=密度导出：R=(m·Γ(d/2+1)/(ρ·π^(d/2)))^(1/d)；⚠ ctypes argtypes 必须与 C 签名对齐。
- 点粒子测试一律用密度 1e9（ρ=100 时 M=1000 中央核 R=1.78 会吞掉 r=1 轨道）。
- **暗物质晕（第 11 轮）**：随机云中央天体改 `fixed`（`gs_add_body_ex`，
  BodyFlag bit0）：不接收力、Verlet 跳过、**不并合**、反弹按逆质量（inv=0
  ⇒ 无穷大质量，只反射对方）。紫色 `#b06cff`（渲染按 state 的 flags&1）。
  ⚠ **getState 步长 7 → 8**（第 8 分量 = flags）：波及 test_main / core_wrap /
  index.html / three.html / verify_g.py，漏改即静默错位。
  ⚠ 晕固定后**总动量不守恒**（晕施加外力，第三定律被破坏）；角动量仍守恒
  （力始终指向原点，力矩为 0）。
- 已取消：第 11 轮曾加「两颗配平质点」抵消总动量与偏差角动量，因速度必然
  ∝√N（两颗轻质点独自扛起整盘的有限 N 抽样净动量）总是飞出视野而移除。
  结论记录备查：最小范数解形如 v=m(c+ω×r)；**3D 存在几何可行性条件
  u·(L − mid×P)=0**（矩阵 κ(|u|²I−uuᵀ) 沿 u=r₂−r₁ 恒奇异），否则无解。

## Barnes-Hut 树（2026-10-14 实现，改动前必读）

- **树节点懒创建**（state 0=空/1=叶/2=内部/3=深度上限桶），反扫传递统计
  （子节点索引恒大于父节点）。⚠ **反扫不变量：com 是「原始加权和」Σ(pos·m)**，
  父子直接相加、最后统一除质量——历史上这里有过双重加权 bug
  （子节点推 com·m 而 com 已是加权和 → com 沿树炸到 1e8 → 单元永不开启 →
  每体 1300+ 交互，比直接求和还慢；修后 2D 满盘 θ=0.6 每体 ~67 次）。
- θ 判据用**度量距离** `mrad < θ·d`（mrad=com→角点最大度量距离）。
  默认 θ=0.6（`setTree(treeMode, theta)`：0=auto(N≥1500) 1=off）。
  `kTreeThreshold=1500`、`kMaxDepth=48`、`gs_set_tree(h, mode, theta)`。
- **PE 估计器用更紧的 θ_PE=min(θ,0.25)**（能量监视精度不跟随力精度下降，
  实测误差 −1.7e-3 @ θ=0.4）。
- 碰撞检测用 `minMetricDistToBox`（可分离度量族的精确盒距离下界：各轴独立 clamp）剪枝。
- **实测基准（单线程 g++ -O2，2D 满盘，θ=0.6）**：N=2000 直接 247ms vs 树 21.6ms（11×），
  N=4000 988ms vs 47ms（21×），32k 树 499ms/步（近线性）。3D 常数大 4×（N=2000 3.3×）。
- 实测力误差（环形域 N=3000）：θ=0.4 → maxErr/aMax 2.1e-2 / mean 3.4e-3；
  θ=0.3 → 8.4e-3 / 1.7e-3；θ=0.2 → 2.2e-3 / 6.1e-4。
- **误差指标教训**：按单体相对误差会被近零合力天体发散带偏（假阳性），
  必须按力尺度 aMax 归一化。
- 力正确但单元不开启的常见原因不是遍历而是 com 统计——先查反扫不变量。

## 服务器/前端性能约定（2026-10-14）

- **速度倍率** `speed` ∈ [0.25, 64]：每 tick 步数 = substeps×speed（≤512）。
  ⚠ **每 tick 墙钟预算 70%**：单步可能超预算 → 高倍率重系统自动降为"能跑多快跑多快"，
  否则 tick 循环持锁 1s+ 会饿死 HTTP/SSE（实测 2500 体 ×16 曾把 API 饿死）。
- **energy() 在大 N 很贵**（树重建 + PE 遍历 ~30ms @2500）→ snapshot 降频缓存
  （n<1000 每 2 帧，大系统每 15 帧）。
- **SSE 帧率随 N 分级**：<1500→30Hz，<3000→12Hz，≥3000→8Hz（JSON 7×N doubles）。
- **前端拖尾自适应**：N>600→40，N>1500→16，N>3000→0（2500 体×240 拖尾曾把
  headless 浏览器渲染吃到 97fps 拖垮服务器；真实 GPU 浏览器压力小得多）。
- genN 上限 4000（getState 缓冲 4096 体）。

## 已知数值特性（实测，非 bug）

- **角动量 L 对任意度量被 Verlet 精确守恒**（中心力对，实测 1e-15，即使 E 爆炸也不动；
  第 9 轮后 f≠t² 也守恒——T12 演示 dL/L=1e-15 同时 dE/E=−2%）。
- **各向异性度量下能量非守恒**（第 9 轮后的本征特性）：f≠t² 时中心力幅随方向变、
  场非保守（curl≠0），E 漂移是模型物理不是积分器误差；f=t² 时严格守恒
  （圆轨道 ~1e-13，dt=2e-4）。旧「L¹ 能量漂移随 dt 线性」结论随语义变更作废
  （drift_test.py 已切到新语义，历史数字不再适用）。
- **核内刚性发散**：密度 ρ 大（半径小）→ 核内线性力简谐频率 ω~√(F(R)/R)，
  ω·dt>2 时 Verlet 指数发散。ρ=1e9/dt=5e-4 必炸，ρ=1e4 安全；调大密度必须同步调小 dt。
- **自引力盘弛豫**：轨道变扁是真实二体散射（L、E 守恒下的重新分配），随机云已
  近平衡初始化（初始 ⟨e⟩≈0.13）但几十动力时标后仍会弛豫；
  诊断脚本 `server/eccentricity_test.py`。
- 树模式下长时 N 体积分会形成硬双星/三体弹弓（真实物理），dt 解析不了时能量失控——
  测试里用弱 G + 冷盘 + minSep 采样规避（t≪t_ff）。

## 验证工具

- `bin/gs_test.exe`：**57 项**断言（含 T2/T3 中心力方向/幅值、T12 各向异性角动量、
  T11 树精度/PE/冷盘积分、T14 钉住晕 5 项），全过为准。
- `bin/gs_bench.exe`：直接 vs 树基准（`--3d`、`--direct-max`、`--steps`）。
  第 9 轮后力循环更快：N=2000 直接 22ms/步、树 1.8ms/步；N=4000 树 4.3ms/步；32k 树 57ms。

### ⚠ 2026-10-15 订正：Python 验证脚本已全部删除
用户要求「移除测试脚本」后按选定档位删除，以下条目**在磁盘上已不存在**，
需要时须重写（原记录保留备查，勿再当可用工具引用）：
- ~~`server/smoke_test.py`、`server/eccentricity_test.py`、`server/drift_test.py`~~
- ~~`server/verify_pages.py`~~（两页 headless 冒烟 + 截图）
- ~~`server/verify_g.py`~~（G 标定 E2E：默认 G=1 / R=20→G=4 / F=r / LJ 近平衡）
- ~~`server/verify_halo.py`~~（纯晕验证 10/10：flags=1 / m 不变 / |x|=|v|=0 /
  体数=N+1 / 盘 L_z 非零 / 演化后不动 / 并合模式下质量不变）
- ~~`server/verify_speed_tree.py`~~（加速 + 树 E2E）
- ~~`server/screenshot.py`~~（本就过时，引用第 8 轮前的 #preset 控件）
同时删除全部 `shot_*.png`（server/ 10 个 + bin/ 11 个 + web/ 8 个）。
**当前仓库内唯一可跑的验证是 `bin/gs_test.exe`（需先 build.bat）。**
（脚本均需 venv playwright：`...\envs\default\Scripts\python.exe`）

### 重写这些脚本时必须记住的坑
- ⚠ 截图/联调**不能并行跑**：多个浏览器共享同一服务器模拟实例，会互相抢配置。
  第 12 轮因此误报过一次「LJ 云 drift 170%」，串行重跑即 0.1%。
- ⚠ 写浏览器端诊断时**按 `s.n` 遍历，不要拿 `i+k < a.length` 判界**——stride 8 下
  写成 `i+2` 会读到数组外，看到一片 `undefined` 却误判成「SSE 数组有洞」。
- 第 8 轮 UI 改成设置抽屉后：必须先 `click("#settingsToggle")`，否则 #genBtn
  不可见 → 30s 超时。
- **`#genN` 是 log10 刻度**（value=log10 N，范围 0.3–3.602），直接赋 2500 会被
  钳到 4000，要赋 `Math.log10(2500)`。
- E2E 暂停要用前端 `sim.running` 判位后点 #playBtn（不是按钮文本、也不是
  直接 POST /api/control，后者会被 applySetup 的 running:true 覆盖）。
