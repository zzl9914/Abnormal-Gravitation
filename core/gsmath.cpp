// 原生构建专用的数学实现（wasm 构建不编译本文件：同名符号由 JS import 提供）。
#include "gsmath.h"

#include <cmath>

extern "C" {

double gs_pow(double a, double b) { return std::pow(a, b); }
double gs_log(double x) { return std::log(x); }
double gs_exp(double x) { return std::exp(x); }
double gs_sin(double x) { return std::sin(x); }
double gs_cos(double x) { return std::cos(x); }
double gs_tan(double x) { return std::tan(x); }
double gs_atan(double x) { return std::atan(x); }

}  // extern "C"
