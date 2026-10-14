// C API boundary for ctypes / FFI consumers.
// All handles are opaque; error text is available via gs_last_error().
// Calls are NOT thread-safe per handle: the Python server serializes them
// behind a lock.

#include "sim.hpp"

#include <cstring>
#include <string>

#if defined(_WIN32)
#define GS_API __declspec(dllexport)
#else
#define GS_API __attribute__((visibility("default")))
#endif

using grav::Simulation;

static char g_err[1024] = "";

static void set_err(const std::string& s) {
    std::strncpy(g_err, s.c_str(), sizeof(g_err) - 1);
    g_err[sizeof(g_err) - 1] = '\0';
}

extern "C" {

GS_API const char* gs_version() {
    return "gravity_core 1.0";
}

GS_API const char* gs_last_error() {
    return g_err;
}

GS_API void* gs_create(int dims, const char* fExpr, const char* vExpr) {
    set_err("");
    if (!fExpr || !vExpr) {
        set_err("表达式指针为空");
        return nullptr;
    }
    Simulation* s = new Simulation(dims, fExpr, vExpr);
    if (!s->ok()) {
        set_err(s->error());
        delete s;
        return nullptr;
    }
    return s;
}

GS_API void gs_destroy(void* h) {
    delete static_cast<Simulation*>(h);
}

// collisionMode: 0 = pass through (core-regularized), 1 = merge, 2 = bounce.
// density: matter density used to derive body radii.
GS_API int gs_configure(void* h, int collisionMode, double density) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s) return -1;
    int m = collisionMode;
    if (m < 0 || m > 2) m = 0;
    s->setCollision(static_cast<grav::Collision>(m), density);
    return 0;
}

// Coupling constant of the force law (scales the force and the derived
// potential).  G must be > 0 and finite; invalid values are ignored.
GS_API int gs_set_G(void* h, double G) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s) return -1;
    s->setG(G);
    return 0;
}

// Raw force-law value F(r) at metric distance r (no G factor).
GS_API double gs_force_at(void* h, double r) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s) return 0.0;
    return s->forceAt(r);
}

// treeMode: 0 = auto (Barnes-Hut when alive >= 1500), 1 = direct summation.
// theta: opening angle, smaller = more accurate (typical 0.3..0.8).
GS_API int gs_set_tree(void* h, int treeMode, double theta) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s) return -1;
    s->setTree(treeMode, theta);
    return 0;
}

GS_API int gs_add_body(void* h, double mass, const double* pos, const double* vel) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s || !pos || !vel) return -1;
    set_err("");
    int idx = s->addBody(mass, pos, vel);
    if (idx < 0) set_err("质点参数无效（质量必须为正且坐标/速度为有限值）");
    return idx;
}

// `fixed` != 0 pins the body (dark-matter halo): it exerts force but never
// moves, never merges, and reflects others elastically in Bounce mode.
GS_API int gs_add_body_ex(void* h, double mass, const double* pos, const double* vel,
                          int fixed) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s || !pos || !vel) return -1;
    set_err("");
    int idx = s->addBody(mass, pos, vel, fixed != 0);
    if (idx < 0) set_err("质点参数无效（质量必须为正且坐标/速度为有限值）");
    return idx;
}

GS_API void gs_clear_bodies(void* h) {
    Simulation* s = static_cast<Simulation*>(h);
    if (s) s->clearBodies();
}

GS_API int gs_step(void* h, double dt, int substeps) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s) return -1;
    set_err("");
    if (!s->step(dt, substeps)) {
        set_err("步进失败：dt 必须为正的有限值");
        return -1;
    }
    return 0;
}

GS_API int gs_count(void* h) {
    Simulation* s = static_cast<Simulation*>(h);
    return s ? s->aliveCount() : 0;
}

GS_API double gs_time(void* h) {
    Simulation* s = static_cast<Simulation*>(h);
    return s ? s->time() : 0.0;
}

GS_API double gs_energy(void* h) {
    Simulation* s = static_cast<Simulation*>(h);
    return s ? s->energy() : 0.0;
}

GS_API int gs_dims(void* h) {
    Simulation* s = static_cast<Simulation*>(h);
    return s ? s->dims() : 0;
}

// out must have room for 8 * maxBodies doubles
// ([mass, x, y, z, vx, vy, vz, flags] per alive body).
// Returns number of bodies written; each record: [mass, x, y, z, vx, vy, vz].
GS_API int gs_get_state(void* h, double* out, int maxBodies) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s || !out || maxBodies <= 0) return 0;
    return s->getState(out, maxBodies);
}

// n samples, log-spaced r in [rmin, rmax]. outR/outV: n doubles each.
GS_API int gs_sample_potential(void* h, int n, double rmin, double rmax,
                               double* outR, double* outV) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s || !outR || !outV || n <= 0) return 0;
    s->samplePotential(n, rmin, rmax, outR, outV);
    return n;
}

// n points on the metric unit ball (first two axes). outXY: 2n doubles.
GS_API int gs_unit_ball(void* h, int n, double* outXY) {
    Simulation* s = static_cast<Simulation*>(h);
    if (!s || !outXY || n <= 0) return 0;
    return s->unitBall2D(n, outXY);
}

}  // extern "C"
