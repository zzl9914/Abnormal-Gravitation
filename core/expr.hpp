#pragma once
// 对偶数表达式引擎（freestanding 友好）。
//
// 解析用户公式（变量为 "r" 或 "t"）成 RPN 程序，可同时求值与一阶精确导数
// （前向模式自动微分）。
//
// 语法（优先级由低到高）：
//   expr  := term (('+'|'-') term)*
//   term  := unary (('*'|'/') unary)*
//   unary := ('-'|'+') unary | power
//   power := atom ('^' unary)?          // 右结合，指数可带符号
//   atom  := number | const | var | func '(' expr ')' | '(' expr ')'
//
// 内置函数：sqrt exp ln log abs sgn sin cos tan atan；常量：pi e
//
// ⚠ freestanding 约束：不使用 std::string / std::vector / new / 异常。
//   token 与指令都是固定容量数组；解析时**直接发射 RPN**（不再建 AST）；
//   默认构造必须保持平凡（静态存储靠 .bss 清零，没有运行时初始化器）。

#include <stdint.h>

#include "gsmath.h"

namespace grav {

// 定长容量（超出即报错，不动态扩容）
enum {
    kMaxInsn = 256,     // RPN 指令数
    kMaxTok = 256,      // token 数
    kMaxExprStack = 130,
    kMaxPowTerms = 16,  // 幂和项数
    kMaxErr = 192       // 错误信息字符数
};

struct Dual {
    double v;  // 值
    double d;  // 对变量的导数
};

enum class Op : uint8_t {
    PushVar,
    PushConst,
    Add, Sub, Mul, Div, Neg,
    PowC,      // u ^ 常数   (insn.num = 指数)
    PowGen,    // u ^ w      （两边都是对偶数）
    Sqrt, Exp, Ln, Abs, Sgn, Sin, Cos, Tan, Atan
};

struct Insn {
    Op op;
    double num;  // PushConst 的值 / PowC 的指数
};

class Expr {
public:
    Expr() {}  // 平凡：静态对象由 .bss 清零
    Expr(const char* text, const char* varName) { compile(text, varName); }

    bool ok() const { return ok_; }
    const char* error() const { return err_; }

    // 只求值
    double value(double x) const { return eval(x).v; }
    // 值 + 精确一阶导数
    Dual eval(double x) const;

    // 快路径识别：表达式恰好是 var^p（p != 0）或 var（p = 1）
    bool asPowerOfVar(double& p) const;

    // 分解成幂和 sum a_i * var^p_i（可用于从 F 解析反积分出 V）。
    // 成功返回 true，项数写入 n（<= cap）。
    bool asPowerSum(double* outP, double* outA, int cap, int& n) const;

    bool triviallyConstant() const;  // 程序只有一条 PushConst

    // 编译（可重复调用，覆盖上一次的程序）。Simulation::init() 会调用它。
    bool compile(const char* text, const char* varName);

private:

    Insn prog_[kMaxInsn];
    int progLen_;
    int maxDepth_;
    bool ok_;
    char err_[kMaxErr];
};

}  // namespace grav
