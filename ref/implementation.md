# 实现说明（Implementation Notes）

面向**要改代码的人**。README 讲的是「它是什么、怎么跑」，本文件讲的是
「它怎么实现、改动时哪些地方会静默出错」。所有常量、函数与不变量均取自
`core/` 当前源码。

---

## 0. 总体形态：JS 直连 wasm 内核

```
core/*.cpp ──clang --target=wasm32 -nostdlib──▶ gravity_core.wasm（~54 KB）
                                                      │ base64 内联
                                                      ▼
                                          web/gravity_wasm.js（零 fetch）
                                                      │ WebAssembly.Instance
                                                      ▼
  index.html / three.html ──▶ web/gravity.js ──▶ gs_create / gs_step / gs_state_fill …
                                     ▲
                                     └─ 状态帧 = Float64Array 视图，指向 wasm 内存（零拷贝）
```

没有 Python、没有 HTTP 接口、没有 SSE、没有 JSON 帧。仿真完全在浏览器里跑，
静态服务器（`web/tools/serve.js`）只负责把文件托管出去。

这条路线决定了内核必须满足 **freestanding 约束**，见第 1 节。

## 1. 代码地图

| 文件 | 行数 | 职责 |
| --- | --- | --- |
| `core/expr.hpp` / `expr.cpp` | ~80 / ~560 | 表达式引擎：词法 → 递归下降**直接发射 RPN**；对偶数前向自动微分；符号分解 `asPowerOfVar` / `asPowerSum` |
| `core/sim.hpp` / `sim.cpp` | ~210 / ~900 | 模拟器：度量核、中心力、派生势、速度 Verlet、碰撞三模式、Barnes-Hut 树、能量监控 |
| `core/api.cpp` | ~200 | C ABI（`gs_*`）：原生 FFI 与 wasm 共用同一套符号；含 wasm 零拷贝辅助函数 |
| `core/gsmath.h` / `gsmath.cpp` | — | 数学抽象层：`gs_pow/gs_log/...` 原生用 libm，**wasm 由 JS 的 `Math.*` 提供** |
| `core/gsstr.h` | — | 最小字符工具（`gs_strtod` / `gs_strcat*` / `isdigit`…），替代 `<cstring>/<cctype>/<cstdlib>` |
| `core/gsrt_wasm.cpp` | — | freestanding 运行时：自己实现 `memcpy / memset / memmove`（仅 wasm 构建编译） |
| `core/test_main.cpp` | 741 | 断言测试（当前 **63** 项，全过为准） |
| `core/bench_main.cpp` | 97 | 直接求和 vs 树的性能基准 |
| `web/gravity.js` | ~330 | JS ↔ 内核桥接：零拷贝状态帧、setup/control/标定、rAF tick 循环 |
| `web/tools/build-kernel.js` | — | clang → wasm → base64 JS，并校验 import 表 |
| `web/tools/check-bridge.js` | — | 桥接层自检（Node，无需浏览器）：**19** 项检查 |
| `web/tools/serve.js` | — | 极简静态服务器（无业务逻辑） |
| `web/index.html` | 837 | 主页：度量 f 预设 5、力律 F 预设 4、设置抽屉、运行控制台 |
| `web/three.html` | 716 | 三体专页：锁定 f=t²、F=1/r²（2D）、场景库、交互布置 |

### 1.1 freestanding 约束（改代码前必读）

内核要能被 `clang --target=wasm32 -nostdlib -ffreestanding -fno-exceptions -fno-rtti`
编译，因此：

1. **不用 STL 容器**：没有 `std::vector` / `std::string` / `std::pair`。
   一律固定容量数组 + 计数（`kMaxBodies = 4096`、`kMaxNodes = 16384`、
   `kMaxExtra = 2048`、`kMaxContacts = 4096`、`kTabN = 4096`）。
2. **不用 `new` / `delete`**：`Simulation` 放在 `api.cpp` 的静态槽位数组里
   （slot 0 = 主模拟，slot 1 = 一次性探针）；句柄 = 槽位索引 + 1。
3. **默认构造必须平凡**：所有成员都**没有**默认初始化器（否则会生成需要运行时
   执行的动态初始化，而 freestanding 没有 `__wasm_call_ctors`）。语义默认值
   （`density=100`、`G=1`、`theta=0.6`…）一律在 `init()` 里设置。静态实例靠
   `.bss` 清零，所以固定容量大数组**不占 wasm 文件体积**。
4. **数学走 `gsmath.h`**：`gs_pow / gs_log / gs_exp / gs_sin / gs_cos / gs_tan /
   gs_atan` 是 `extern "C"`，wasm 侧由 JS 提供（`Math.*`）；
   `sqrt / fabs / copysign` 用编译器内建，直接落到 wasm 原生指令；
   `gs_tgamma` 用 Lanczos 近似自己实现（两侧同一份代码，避免 import）。
5. **只 include clang 自带的 C 头**：`<stdint.h>` / `<stddef.h>`，**不能**用
   `<cstdint>` / `<cstddef>`（wasm32 目标下找不到 libstdc++ 的 C++ 头）。
6. 表达式引擎**不建 AST**，递归下降时直接发射 RPN（省掉 `new`）。

共享静态工作缓冲（同一时刻只有一个树在活动，每次都重建）：

- `g_nodes[kMaxNodes]` — 树节点（放在文件作用域，因为两个 Simulation 槽位
  各存一份会太大）
- `g_extraBody / g_extraNext` — bucket 的溢出链表（`TreeNode::extra` 是头索引，
  取代原来的 `std::vector<int>`）
- `g_cA / g_cB` — 收集到的接触对

### 1.2 与 JS 的数据交换

- **字符串参数**：JS 用 `TextEncoder` 写进 `gs_str_buf()`（1 KB），再传指针。
- **double 数组参数**（pos/vel）：写进 `gs_dscratch()`（64 个 double）。
- **状态帧**：`gs_state_fill()` 把 `[m,x,y,z,vx,vy,vz,flags]×N` 写进 `g_state`，
  返回体数；`gs_state_buf()` 给指针；JS 直接
  `new Float64Array(memory.buffer, ptr, n*8)` —— **零拷贝**。
- **错误信息**：`gs_last_error()` 返回指针，JS 用 `TextDecoder` 解 UTF-8。
- wasm 侧没有 malloc，JS 也不在 wasm 堆上分配，全部走上面这些静态缓冲。

构建参数要点（`web/tools/build-kernel.js`）：

```
--target=wasm32 -O3 -nostdlib -ffreestanding -fno-exceptions -fno-rtti
-Wl,--no-entry -Wl,--allow-undefined -Wl,--export-memory
-Wl,--initial-memory=8388608 -Wl,--max-memory=268435456 -Wl,-zstack-size=1048576
-Wl,--export=gs_create ...
```

⚠ 栈要 1 MB：`asPowerSum` 的符号栈是栈上数组（~34 KB），插入又是递归的。
⚠ 原生侧反过来：`Simulation` 约 1 MB，测试里是**栈对象**，所以
`gs_test.exe` / `gs_bench.exe` 链接时必须 `-Wl,--stack,67108864`（否则一启动就崩）。

关键常量（`sim.cpp` / `sim.hpp`）：

- `kMaxDepth = 48` — 树深上限，超限转 bucket（state 3）
- `kTreeThreshold = 1500` — auto 模式启用树的存活体数阈值
- `kGradClamp = 1e6` — 度量梯度分量钳位
- `kForceClamp = 1e12` — 力幅钳位
- 维度上限 8（前端只开放 2D / 3D）

---

## 2. 表达式引擎（`expr.*`）

### 2.1 词法与语法

`lex()` 支持十进制小数、科学计数法（`1e-3`、`.5`）、标识符与运算符。递归下降按
（优先级由低到高）：

```
expr  := term (('+'|'-') term)*
term  := unary (('*'|'/') unary)*
unary := ('-'|'+') unary | power
power := atom ('^' unary)?        // 右结合；指数可带符号，如 t^-2
atom  := number | const | var | func '(' expr ')' | '(' expr ')'
```

内置函数：`sqrt exp ln(log) abs sgn sin cos tan atan`；常量 `pi e`。
错误一律非抛出：构造后 `ok() == false`，`error()` 给出中文原因（含出错位置）。

### 2.2 RPN 编译

`Compiler` 把 AST 走成 `std::vector<Insn>` 的栈式程序，同时跟踪栈深：

- `push1` / `net0` / `binary` 三类指令分别对应净 +1 / 净 0 / 净 −1
- 编译结束要求 `depth == 1 && maxDepth > 0`，且 `maxDepth ≤ 128`，否则报
  「表达式结构无效 / 嵌套过深」

⚠ **幂的特化**：只有「指数是字面常数」时编译为 `PowC`（指数折叠进 `insn.num`），
否则为 `PowGen`（底数与指数都走对偶数）。这个分岔同时被 `asPowerOfVar`（识别
`[PushVar]` 或 `[PushVar, PowC]`）和 `asPowerSum` 依赖，改动编译顺序会连带破坏
度量快路径与解析反积分。

### 2.3 对偶数求值（前向 AD）

`eval(x)` 用固定大小的栈（`Dual st[130]`，与 maxDepth 128 对齐）执行 RPN，
每个算子同时推进值与导数。边界处理集中在 `eval` 内：

- `Div`：除零 → `0/0` 得 NaN，非零除零得带符号 Inf
- `PowC`：`a != 0` 走常规幂律；`a == 0` 按 c 与 1 的关系决定 0 / ±Inf
- `PowGen`：`a > 0` 时用 `u^w·(w' ln u + w·u'/u)`；`0^0` 定义为 1；其余 NaN
- `Sqrt / Ln`：负宗量 → NaN；`Sqrt(0)` 的导数在 `a.d != 0` 时为 ±Inf

### 2.4 符号分解 `asPowerSum`

在 RPN 上跑一个符号栈 `PolySum`（term 列表 `(p, a)` 表示 `a·r^p`），保持「幂和」
这个表示在运算下**闭合**，一旦不闭合就返回 false：

- `Add / Sub`：逐项合并同类幂，系数为 0 的项删除
- `Mul`：笛卡尔积，指数相加、系数相乘，再合并同类项（term 数超过 `kCap = 64` 视为失败）
- `Div`：要求分母是**单项**，逐项指数相减、系数相除
- `PowC` / `Sqrt`：要求**单项**，`(a r^p)^c = a^c r^{pc}`（`Sqrt` 即 `c = 0.5`）
- `Abs`：要求至多一项，取系数绝对值
- `PowGen / Exp / Ln / Sgn / Sin / Cos / Tan / Atan`：直接失败 → 走数值表路径

---

## 3. 度量核（`sim.cpp`）

### 3.1 度量距离

```
d(Δ) = F_m( Σ_k f(|Δ_k|) ),   F_m = f⁻¹
```

要求 f 在 `[0,∞)` 上严格单增且 `f(0) = 0`（保证重合点距离为 0、正定性）。

### 3.2 `invertF(S)`：求 f(x) = S

1. `S <= 0` → 0
2. **幂快路径**：若 f = t^p（由 `asPowerOfVar` 识别），直接 `S^(1/p)`
3. 一般 f：**倍增定界**（hi 从 1 起翻倍，上限 200 次；f 有界时视为距离无界）→
   **safeguarded Newton**（导数由 `f_.eval(x)` 的对偶数给出；迭代不在 `(lo,hi)`
   内则退回二分点）
4. 收敛判据：残差 `≤ 1e-14·max(1,S)`，或区间长度 `≤ 1e-15·max(1,hi)`；最多 80 次迭代

### 3.3 梯度 `pairKernel`

- 幂 f = t^p：`grad_k = (|Δ_k|/d)^(p−1) · sgn(Δ_k)`
- 一般 f：`grad_k = f'(|Δ_k|)/f'(d) · sgn(Δ_k)`（链式法则，`fpd` 做了 1e-300 下限保护）

均钳位到 `±1e6`。⚠ 力的**方向**不用这个梯度：见下节。

### 3.4 可分离性 → 精确的盒距离下界

度量族 `d = F_m(Σ f(|Δ_k|))` 是**逐轴可分离**的，因此「点到轴对齐盒的最小度量距离」
可以逐轴 clamp 后精确求得（`minMetricDistToBox`）。碰撞检测用它剪枝，
所以**碰撞检测是精确的，不是近似**（与树近似的力计算不同）。

---

## 4. 力律与派生势

### 4.1 中心力

```
a_i = G · m_j · F(d) · Δ̂ ,   Δ = x_j − x_i ,  Δ̂ = Δ/|Δ|
```

`F > 0` = 吸引（指向施力对象），`F < 0` = 排斥。力**恒沿连线**，
方向用欧氏 Δ̂ 而非度量梯度——这正是角动量对任意度量精确守恒的原因。
构造时只在 `r = 0.5, 1, 2` 三个探测点检查 F 有限；其余奇点由运行时钳位与
核正则化兜底。

### 4.2 派生势 V：解析路径

`V = ∫F dr`，即 **V′ = F**（注意符号约定：不是 −V′）。

若 `F_.asPowerSum()` 成功（`F = Σ a_i r^{p_i}`），则逐项解析积分：

- `p_i = −1` → 累加到 `vLn_`，贡献 `vLn_·ln r`
- 否则 → `vPow_ += (p_i+1, a_i/(p_i+1))`

### 4.3 派生势 V：数值表路径

`buildPotentialTable()`：

- 网格：`[1e-6, 1e6]` 上 **4096** 个对数节点
- 每段在 ln r 空间用 **4 点 Gauss–Legendre** 积分（`dr = r·dln r`）
- 节点同时存 `r`、`V`、`F`
- 查询：ln r 上的**三次 Hermite 插值**，节点斜率取精确值 `dV/dln r = F(r)·r`
- 网格外：用端点斜率线性外推
- 溢出保护：`|V| > 1e300` 或 NaN 时钳到 `±1e300`

### 4.4 核内正则化 `pairV`

`d < R = R_i + R_j` 时，力幅替换为二次势的导数，二者在 `d = R` 处连续：

```
W(x) = V(R) + G·F(R)·(x² − R²) / (2R)
```

力在核内随 d 线性衰减到 0，**r → 0 的奇点永不被求值**（只需求 `V(R)` 与 `F(R)`）。

### 4.5 耦合常数 G

`setG(G)` 仅在 `G > 0` 且有限时生效；`potentialAt` 在返回前统一乘 `G_`。
自动标定目标 `G·F(R/2) = 0.04`（`kGRef`，基准为默认 R=10、F=1/r²、G=1 的配置），
换力律或盘半径后重新「应用」即可，保证演化时标可比。

---

## 5. 积分器

`step(dt, substeps)` 是标准 kick–drift–kick（速度 Verlet，辛积分）：

```
if (!accValid_) computeAccel();
for s in 1..substeps:
    v += 0.5·h·a ;  x += h·v        // half kick + drift
    computeAccel()                    // new a
    v += 0.5·h·a                      // half kick
    handleCollisions()
    t += h
```

- `accValid_` 让「刚算过加速度」这一步免于重复计算；`addBody` / `clearBodies` 会置 false
- 死体（`!alive`）与钉住体（`fixed`）在所有循环中跳过

---

## 6. 碰撞

**半径由物质密度导出**（2D 圆盘 / 3D 球体）：

```
R = ( m·Γ(d/2+1) / (ρ·π^(d/2)) )^(1/d)
```

`setCollision` 改密度会即时重算所有半径。接触判据用**度量距离** `d ≤ R_i + R_j`。

三模式（`Collision::{Pass, Merge, Bounce}`）：

- **Pass**：直接返回，靠核正则化平滑穿过（默认）
- **Merge**：合并 j → i。位置/速度按质量加权（动量守恒），质量相加、半径重算，
  `B.alive = false`，随后 `computeAccel()` 重算（树模式下会连带重建树）
- **Bounce**：法向取**两心连线方向**（与中心力模型一致；f=t² 时与度量梯度法向相同）。
  接近速度 `w = (v_A − v_B)·n̂`，`w ≤ 0`（已在分离）则不动；以逆质量分配冲量：
  `invA = A.fixed ? 0 : 1/m_A`，`tot = invA + invB`，`e = 1`，
  `v_A −= (1+e)(invA/tot)w·n̂`，`v_B += (1+e)(invB/tot)w·n̂`。两者都钉住则直接返回

⚠ **钉住天体（暗物质晕）**：`fixed` 体在 `applyOneSided` 开头就 return（施力不受力）、
在 `treeForce` 与 Verlet 循环中被跳过、**不参与并合**（Merge 分支里 `A.fixed || B.fixed`
直接返回，退化成穿过）、反弹时按**逆质量 0**（无穷大惯性，只反射对方）。
副作用：**总动量不再守恒**（晕施加外力，第三定律被破坏），但角动量仍守恒
（力始终指向原点，力矩为 0）。

---

## 7. Barnes-Hut 树

### 7.1 节点与状态机

`TreeNode` 固定 8 个子指针（2D 用前 4 个，3D 用 8 个）：

| state | 含义 |
| --- | --- |
| 0 | 空 |
| 1 | 叶（恰好一个天体，`body` 有效） |
| 2 | 内部节点（`child[]` 有效） |
| 3 | bucket — 深度达 `kMaxDepth` 仍重合，用 `extra[]` 兜住 |

字段：`mass`、`com[]`、`ctr[]`、`half[]`、`mrad`、`maxRad`、`parent`、`state`、`extra`。

### 7.2 建树（`buildTree`）

1. 求所有存活体的包围盒；`half = 0.5·(mx−mn)·1.000001 + 1e-12`
   （略微外扩，避免边界点反复落在格子线上）
2. **子节点懒创建**：只在真正要下沉时才 `newNode` 并挂到 `child[k]`
3. 插入：空 → 叶；叶 → 分裂并把旧天体下推；内部 → 按各轴中点（`>= ctr` 取 1）选象限位

**反扫累积**（⚠ 最容易写错的地方）：

- 因 `newNode` 追加，子节点索引恒大于父节点，故**倒序遍历**即可保证「子先于父结算」
- **不变量：扫过程中 `nd.com` 保存的是原始加权和 `Σ(pos·m)`，不是质心**
- 父子之间直接把该加权和原样相加，全部累积完毕后才统一除以质量
- 历史 bug：子节点把「已算好的质心」再乘一次质量推给父节点 → 双重加权 →
  `com` 沿树炸到 1e8 → 单元永远打不开 → 每体 1300+ 次交互，比直接求和还慢

最后计算 `mrad`：遍历 `2^d` 个角点，取 `max d(com, corner)`（度量距离）。

### 7.3 力遍历（`treeForce`）

显式栈（`int stack[512]`，入栈时检查 `sp < 511`）：

- 叶 / bucket：跳过自身（`nd.body != i`）后逐体调用 `applyOneSided`
- 内部节点：判据 **`mrad < θ·d`** 成立 → 用质心整体近似；否则压入所有非空子节点

`applyOneSided(i, posj, mj, radj)`：

```
delta = posj − x_i ; rho2 = |delta|²
if (rho2 <= 0) return          // 重合，无方向
if (bodies_[i].fixed) return   // 钉住体不受力
d = metricDistance(delta)
s = (d < Rij) ? G·F(Rij)·d/Rij : G·F(d)
s 钳到 ±1e12；w = s/√rho2
acc += mj · w · delta
```

### 7.4 势能估计器（`treePotential`）

沿用同样的遍历，但用**更紧的判据** `θ_PE = min(θ, 0.25)`，因为 PE 每帧只算一次，
而力每子步都算——让能量监视的精度不跟随力精度下降。`U = 0.5·Σ_i U_i`。

### 7.5 接触收集（`collectContacts`）

对每个 i 走树，内部节点用
`minMetricDistToBox(x_i, nd) > R_i + nd.maxRad` 剪枝（该下界是精确的，故不漏对）。
树模式下每次合并都会触发 `computeAccel()`（重建树），因此收集到的 `(i,j)` 对在
`resolveContact` 内部还会重新校验 `alive` 与距离。

### 7.6 已知代价

- 单极近似 → **动量仅近似守恒**（与直接求和不同）
- 树遍历有 ~0.4% 量级的力噪声（θ=0.6），会加速 N 体盘的弛豫
- `treeActive()` 每次都 O(N) 数存活体；`computeAccel` 每步都重建树

---

## 8. 运行调度（`web/gravity.js` 的 tick 循环）

调度逻辑从原来的 Python 服务端搬到了浏览器里，行为保持一致：

- **物理 tick 固定 30 Hz**（不是渲染帧率）：`requestAnimationFrame` 驱动渲染，
  tick 按 `1000/30` ms 累积推进；单帧最多补 4 个 tick（切回标签页时不会一次性
  补几百毫秒）。
- **每 tick 步数** = `substeps × speed`，夹到 `[1, 512]`。
- **墙钟预算**：`budget = 0.7 × interval`。批次以自适应 chunk 下发：
  chunk 从 1 起，若这批耗时 `< 0.25 × budget` 就 ×4（上限 32），
  每批结束检查预算，超了就停 → 高倍率重系统自动降级为「能跑多快跑多快」。
  ⚠ 这个预算不能省：同步步进会阻塞 rAF，界面直接卡死。
- **能量降频**：`n < 1000` 每 2 帧算一次，否则每 15 帧（树模式下 `energy()`
  要重建树并全量遍历 PE）。
- **零拷贝帧**：`sim.bodies` 每帧新建一个指向 wasm 内存的 `Float64Array`
  子视图（不复制数据），渲染代码按 stride 8 读即可。

`Gravity.setup/control/reset/calibrate/analysis` 与旧的
`/api/setup`、`/api/control`、`/api/reset`、`/api/calibrate`、`/api/analysis`
一一对应，语义不变，只是从 HTTP + JSON 变成了进程内函数调用。
其中 `calibrate` 返回标定后的 `g` 与可选采样表 `fs`，前端据此按
`v = √(G·m_enc·F(r)·r)` 做**任意力律**的近平衡建云。

---

## 9. 前端

无构建、无 CDN、原生 canvas 2D。前端只做两件事：把 UI 状态 `POST /api/setup`，
然后消费 SSE 帧。

- **index.html**：度量 f 预设 5 个（`t^0.5 / t / t^2 / t^3 / t^10`）、
  力律 F 预设 4 个（`r / 1/r / 1/r² / 1/r⁶−1/r¹²`），选预设会同时带入推荐 dt/子步数；
  左下角 fit 圆钮、图例、可收起监视器；右下 dock；3D 随机云在球体内做均匀采样。
- **three.html**：锁定 `f=t²`、`F=1/r²`、2D。经典场景库 6 个 + 点击放置天体 +
  拖拽定初速 + `localStorage`（键 `grav3body.scenes` / `grav3body.last`）+ JSON 导入导出。
- **拖尾自适应**：N > 600 → 40，N > 1500 → 16，N > 3000 → 0
  （2500 体 × 240 拖尾曾把 headless 浏览器的渲染开销吃到 97 fps，反过来拖垮服务器）。
- **genN 上限 4000**（内核 `getState` 缓冲 4096 体）。

---

## 10. 改动禁忌（每一条都是实测踩出来的）

1. **`getState` 步长是 8**，第 8 个分量是 `flags`。改动波及 `test_main.cpp`、
   `web/gravity.js`、`index.html`、`three.html`——漏改一处就是**静默错位**，不报错。
2. **新增内核代码必须遵守 freestanding 约束**（§1.1）：不用 STL 容器、
   不用 `new`、不给成员加默认初始化器、不用 `<cstdint>`。
   一旦违反，wasm 构建会当场失败（或更糟：生成需要运行时初始化的段）。
3. **需要新的数学函数时**，加到 `gsmath.h`：原生在 `gsmath.cpp` 里用 `std::*`
   实现，wasm 侧则要同步在 `web/gravity.js` 的 importObject 里提供，
   否则 `WebAssembly.instantiate` 直接报缺 import。
4. **树反扫的 `com` 不变量**：扫过程中是 `Σ(pos·m)`，最后一趟才除质量。
   力正确但单元打不开，先查这里，别去查遍历逻辑。
5. **`kMaxBodies = 4096` 是硬上限**，`addBody` 超限返回 −1 并设置错误。
   基准脚本已尊重该上限（曾出现 N=32000 只建出 4096 个天体的假加速）。
6. **力的误差指标要按力尺度 `aMax` 归一化**。按单体相对误差统计会被近零合力的
   天体发散带偏，得到假阳性。
7. **原生测试/基准必须 `-Wl,--stack,67108864`**：`Simulation` 是 ~1 MB 的栈对象，
   栈不够时进程一启动就死（表现为无输出、退出码 127，很难查）。
8. 在浏览器端遍历状态数组时**按 `sim.n` 判界**，不要写 `i+k < a.length`——stride 8 下
   写成 `i+2` 会读到数组之外。
9. `#genN` 是 **log10 刻度**（value = log10 N，范围 0.3–3.602），直接赋 2500 会被钳到 4000。
10. **boot 顺序**：`sim.running = true` → `await Gravity.ready()` →
    `connectStream()` → `applySetup()`。内核没 ready 就 setup 会拿到
    「内核尚未初始化」，只改按钮文字则会启动即暂停。

---

## 11. 验证与基准

```bat
build.bat            :: 原生：内核 DLL + gs_test.exe + gs_bench.exe（需要 g++）
bin\gs_test.exe      :: 断言测试，输出 "ALL TESTS PASSED (0 failures)"
bin\gs_bench.exe [--3d] [--direct-max N] [--steps N]

build_wasm.bat       :: wasm：core/*.cpp → web/gravity_wasm.js（需要 clang + wasm32）
node web\tools\check-bridge.js   :: 桥接层自检（19 项，含圆轨道守恒/标定/单位球/错误路径）
start.bat            :: 静态服务器 → http://127.0.0.1:8123
```

`check-bridge.js` 是 wasm 路径最快的回归手段：它覆盖 setup（含 fixed 晕）、
圆轨道能量守恒、状态帧 stride/flags、G 标定、F 采样表、单位球形状、
reset、非法 f/F 的拒绝，以及 2000 体的树模式耗时。

测试覆盖面（按 `test_main.cpp` 分组）：T1 圆轨道能量守恒；T2/T3 各向异性度量下
中心力方向与幅值；T2b 一般 F 的「力 = −∇(∫F)」有限差分校验；T4 度量梯度对拍；
T5 单位球形状（L¹ 为菱形）与 V 的解析反积分（含 ln 项与 LJ）；T6/T9 并合与反弹的
动量守恒；T7 三体 8 字；T8 非法 f / F 的拒绝；**T11 树精度（maxErr/aMax、误差分布）、
PE 估计器误差、冷盘积分健康**；**T12 各向异性度量下角动量守恒（1e-15）**；
**T13 G 的缩放与标定不变量**；**T14 钉住晕的 12 项行为**。

⚠ 基准耗时参考（本机单线程）：`--direct-max 2000 --steps 10` 约 4 分钟，
完整一轮（默认 `--direct-max 4000`）约 5 分钟起步。
⚠ bench 里 N=500/1000 的「tree」行其实是直接求和（auto 阈值 1500，树未开）。
