#pragma once
// 数学抽象层。
//
// 目的：让 core/ 同时能编译成「原生 DLL」与「freestanding wasm32」。
//
//   - 原生构建：gs_pow / gs_log / ... 由 gsmath.cpp 用 libm 实现。
//   - wasm 构建（-nostdlib）：这些符号保持未定义，链接时用
//     -Wl,--allow-undefined 变成 wasm import，由 JS 侧用 Math.* 提供。
//
// 能映射到 wasm 原生指令的（sqrt / fabs / copysign）一律用编译器内建，
// 不进 import 表，避免 JS↔wasm 调用开销。

extern "C" {
double gs_pow(double a, double b);
double gs_log(double x);
double gs_exp(double x);
double gs_sin(double x);
double gs_cos(double x);
double gs_tan(double x);
double gs_atan(double x);
}

static const double GS_PI = 3.14159265358979323846;
static const double GS_SQRT_2PI = 2.50662827463100050242;

// ---- 内建：直接落到 wasm 指令 ----

inline double gs_sqrt(double x) { return __builtin_sqrt(x); }
inline double gs_fabs(double x) { return __builtin_fabs(x); }
inline double gs_copysign(double x, double y) { return __builtin_copysign(x, y); }
inline bool gs_isfinite(double x) {
    return x == x && x != __builtin_inf() && x != -__builtin_inf();
}
inline double gs_min2(double a, double b) { return (b < a) ? b : a; }
inline double gs_max2(double a, double b) { return (a < b) ? b : a; }
inline double gs_clamp(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

// ---- Γ(x)：Lanczos 近似 ----
// 两个构建共用同一份代码，所以 wasm 侧不需要为它开 import，
// 原生与 wasm 的结果也完全一致（半径公式只用到 x = d/2+1，d ∈ [1,8]）。
inline double gs_tgamma(double x) {
    if (x < 0.5) {
        // 反射公式。本工程用不到（dims ≥ 1 ⇒ x ≥ 1.5），仅为完备性保留。
        return GS_PI / (gs_sin(GS_PI * x) * gs_tgamma(1.0 - x));
    }
    static const double kG[9] = {
        0.99999999999980993,      676.5203681218851,     -1259.1392167224028,
        771.32342877765313,       -176.61502916214059,   12.507343278686905,
        -0.13857109526572012,     9.9843695780195716e-6, 1.5056327351493116e-7};
    x -= 1.0;
    double a = kG[0];
    for (int i = 1; i < 9; i++) a += kG[i] / (x + (double)i);
    double t = x + 7.5;
    return GS_SQRT_2PI * gs_pow(t, x + 0.5) * gs_exp(-t) * a;
}
