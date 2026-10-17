// C ABI 边界：原生侧给 ctypes 用，wasm 侧给 JS 用（同一套 gs_* 符号）。
//
// ⚠ freestanding（wasm）约束：没有 new / delete，Simulation 放在**静态槽位数组**里。
//   slot 0 = 主模拟；slot 1 = 一次性探针（标定 G 时临时建一个力律不同的实例）。
//   ⚠ 同一个句柄上的调用不做线程保护：调用方（Python 服务器 / JS）负责串行化。

#include "sim.hpp"

#include "gsstr.h"

#if defined(_WIN32)
#define GS_API __declspec(dllexport)
#else
#define GS_API __attribute__((visibility("default")))
#endif

using grav::Simulation;

enum { kSlots = 2 };

static Simulation g_slots[kSlots];
static char g_err[256];

// wasm 侧的零拷贝输出缓冲：JS 直接在这段内存上建 Float64Array 视图
enum { kStateDoubles = grav::kMaxBodies * 8 };
enum { kAnaMax = 512 };
enum { kBallMax = 2048 };  // 1024 个点 × 2

static double g_state[kStateDoubles];
static double g_potR[kAnaMax];
static double g_potV[kAnaMax];
static double g_ball[kBallMax];

// JS → wasm 的输入暂存区（wasm 侧没有 malloc，JS 也不好自己分配）：
// 字符串参数先 UTF-8 写进 g_strbuf，再传指针；double 数组写进 g_dscratch。
static char g_strbuf[1024];
static double g_dscratch[64];

static void set_err(const char* s) {
    g_err[0] = '\0';
    if (s) gs_strcat(g_err, (int)sizeof(g_err), s);
}

// 句柄 = slot 索引 + 1（0 / nullptr 表示无效）
static Simulation* simOf(void* h) {
    unsigned long long v = (unsigned long long)(size_t)h;
    if (v < 1 || v > kSlots) return nullptr;
    return &g_slots[v - 1];
}

extern "C" {

GS_API const char* gs_version() { return "gravity_core 2.0 (freestanding-capable)"; }

GS_API const char* gs_last_error() { return g_err; }

GS_API void* gs_create(int dims, const char* fExpr, const char* forceExpr) {
    set_err("");
    if (!fExpr || !forceExpr) {
        set_err("表达式指针为空");
        return nullptr;
    }
    // 优先复用 slot 0（前端同时只需要一个模拟；slot 1 留给一次性探针）
    Simulation* s = &g_slots[0];
    if (!s->init(dims, fExpr, forceExpr)) {
        set_err(s->error());
        return nullptr;
    }
    return (void*)(size_t)1;
}

// 一次性探针：在 slot 1 上建一个实例（标定 / 采样用），不影响 slot 0
GS_API void* gs_create_probe(int dims, const char* fExpr, const char* forceExpr) {
    set_err("");
    if (!fExpr || !forceExpr) {
        set_err("表达式指针为空");
        return nullptr;
    }
    Simulation* s = &g_slots[1];
    if (!s->init(dims, fExpr, forceExpr)) {
        set_err(s->error());
        return nullptr;
    }
    return (void*)(size_t)2;
}

GS_API void gs_destroy(void* h) {
    Simulation* s = simOf(h);
    if (s) s->clearBodies();
}

// collisionMode: 0 = 穿过（核内正则）1 = 并合 2 = 反弹
// density: 导出天体半径用的物质密度
GS_API int gs_configure(void* h, int collisionMode, double density) {
    Simulation* s = simOf(h);
    if (!s) return -1;
    int m = collisionMode;
    if (m < 0 || m > 2) m = 0;
    s->setCollision(static_cast<grav::Collision>(m), density);
    return 0;
}

// 力律耦合常数（同时缩放力与派生势）。G 必须 > 0 且有限，非法值被忽略。
GS_API int gs_set_G(void* h, double G) {
    Simulation* s = simOf(h);
    if (!s) return -1;
    s->setG(G);
    return 0;
}

GS_API double gs_get_G(void* h) {
    Simulation* s = simOf(h);
    return s ? s->g() : 0.0;
}

// 原始力律值 F(r)（不含 G）
GS_API double gs_force_at(void* h, double r) {
    Simulation* s = simOf(h);
    if (!s) return 0.0;
    return s->forceAt(r);
}

// treeMode: 0 = auto（alive >= 1500 用 Barnes-Hut）1 = 直接求和
GS_API int gs_set_tree(void* h, int treeMode, double theta) {
    Simulation* s = simOf(h);
    if (!s) return -1;
    s->setTree(treeMode, theta);
    return 0;
}

// fixed != 0 钉住该天体（暗物质晕）：施力但不动、不并合，反弹时按无穷大质量
GS_API int gs_add_body_ex(void* h, double mass, const double* pos, const double* vel,
                          int fixed) {
    Simulation* s = simOf(h);
    if (!s || !pos || !vel) return -1;
    set_err("");
    int idx = s->addBody(mass, pos, vel, fixed != 0);
    if (idx < 0) set_err("质点参数无效（质量必须为正且坐标/速度为有限值），或已达上限");
    return idx;
}

GS_API int gs_add_body(void* h, double mass, const double* pos, const double* vel) {
    return gs_add_body_ex(h, mass, pos, vel, 0);
}

GS_API void gs_clear_bodies(void* h) {
    Simulation* s = simOf(h);
    if (s) s->clearBodies();
}

GS_API int gs_step(void* h, double dt, int substeps) {
    Simulation* s = simOf(h);
    if (!s) return -1;
    set_err("");
    if (!s->step(dt, substeps)) {
        set_err("步进失败：dt 必须为正的有限值");
        return -1;
    }
    return 0;
}

GS_API int gs_count(void* h) {
    Simulation* s = simOf(h);
    return s ? s->aliveCount() : 0;
}

GS_API double gs_time(void* h) {
    Simulation* s = simOf(h);
    return s ? s->time() : 0.0;
}

GS_API double gs_energy(void* h) {
    Simulation* s = simOf(h);
    return s ? s->energy() : 0.0;
}

GS_API int gs_dims(void* h) {
    Simulation* s = simOf(h);
    return s ? s->dims() : 0;
}

// out 需要 8 * maxBodies 个 double：[mass, x, y, z, vx, vy, vz, flags]
GS_API int gs_get_state(void* h, double* out, int maxBodies) {
    Simulation* s = simOf(h);
    if (!s || !out || maxBodies <= 0) return 0;
    return s->getState(out, maxBodies);
}

GS_API int gs_sample_potential(void* h, int n, double rmin, double rmax, double* outR,
                               double* outV) {
    Simulation* s = simOf(h);
    if (!s || !outR || !outV || n <= 0) return 0;
    s->samplePotential(n, rmin, rmax, outR, outV);
    return n;
}

GS_API int gs_unit_ball(void* h, int n, double* outXY) {
    Simulation* s = simOf(h);
    if (!s || !outXY || n <= 0) return 0;
    return s->unitBall2D(n, outXY);
}

// ------------------------------------------------------- wasm 零拷贝辅助
// JS 侧不自己分配内存：这些函数把结果写进模块内的静态缓冲并回传指针，
// JS 直接在 memory.buffer 上建视图（Float64Array）读取。

GS_API double* gs_state_buf() { return g_state; }

// 把当前状态写进 g_state，返回写入的体数（每体 8 个 double）
GS_API int gs_state_fill(void* h) {
    Simulation* s = simOf(h);
    if (!s) return 0;
    return s->getState(g_state, grav::kMaxBodies);
}

// 采样派生势 V(r)，写入静态缓冲，返回点数
GS_API int gs_analysis(int n, double rmin, double rmax, void* h) {
    Simulation* s = simOf(h);
    if (!s) return 0;
    if (n > kAnaMax) n = kAnaMax;
    if (n < 1) n = 1;
    s->samplePotential(n, rmin, rmax, g_potR, g_potV);
    return n;
}

GS_API double* gs_pot_r_ptr() { return g_potR; }
GS_API double* gs_pot_v_ptr() { return g_potV; }

// 度量单位球，写入静态缓冲，返回点数
GS_API int gs_ball(int n, void* h) {
    Simulation* s = simOf(h);
    if (!s) return 0;
    int pts = n;
    if (pts * 2 > kBallMax) pts = kBallMax / 2;
    if (pts < 3) return 0;
    return s->unitBall2D(pts, g_ball);
}

GS_API double* gs_ball_ptr() { return g_ball; }

// ---- 输入暂存区 ----
GS_API char* gs_str_buf() { return g_strbuf; }
GS_API int gs_str_buf_size() { return (int)sizeof(g_strbuf); }
GS_API double* gs_dscratch() { return g_dscratch; }
GS_API int gs_dscratch_size() { return (int)(sizeof(g_dscratch) / sizeof(g_dscratch[0])); }

}  // extern "C"
