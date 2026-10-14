#pragma once
// Dual-number expression engine.
// Parses a user supplied formula (variable: "r" or "t") into an RPN program
// that can be evaluated for value and exact first derivative simultaneously
// (forward-mode automatic differentiation).
//
// Grammar (precedence low -> high):
//   expr  := term (('+'|'-') term)*
//   term  := unary (('*'|'/') unary)*
//   unary := ('-'|'+') unary | power
//   power := atom ('^' unary)?          // right associative, exponent may be signed
//   atom  := number | const | var | func '(' expr ')' | '(' expr ')'
//
// Built-in functions: sqrt exp ln log abs sgn sin cos tan atan
// Constants: pi e

#include <cstdint>
#include <string>
#include <vector>

namespace grav {

struct Dual {
    double v = 0.0;
    double d = 0.0;  // derivative w.r.t. the variable
    Dual() = default;
    Dual(double val, double der = 0.0) : v(val), d(der) {}
};

enum class Op : uint8_t {
    PushVar,
    PushConst,
    Add, Sub, Mul, Div, Neg,
    PowC,      // u ^ const          (insn.num = exponent)
    PowGen,    // u ^ w  (both operands dual)
    Sqrt, Exp, Ln, Abs, Sgn, Sin, Cos, Tan, Atan
};

struct Insn {
    Op op;
    double num = 0.0;  // PushConst value / PowC exponent
};

class Expr {
public:
    Expr() = default;
    // Non-throwing: on failure ok() == false and error() holds the reason.
    Expr(const std::string& text, const std::string& varName);

    bool ok() const { return ok_; }
    const std::string& error() const { return err_; }

    // Value only.
    double value(double x) const { return eval(x).v; }
    // Value + exact derivative w.r.t. the variable.
    Dual eval(double x) const;

    // Fast-path detection: expression is exactly  var^p  (p != 0) or  var.
    // Returns true and sets p (1.0 for a bare variable).
    bool asPowerOfVar(double& p) const;

    // Decompose into a sum of power terms  sum_i a_i * var^p_i  if the
    // expression is one (affine combinations of powers, products/quotients
    // and integer roots of single terms).  Returns true with terms
    // (p_i, a_i); zero coefficients are dropped.  Used to derive an exact
    // antiderivative of a force law.
    bool asPowerSum(std::vector<std::pair<double, double>>& terms) const;

    bool triviallyConstant() const;  // program is a single PushConst

private:
    bool compile(const std::string& text, const std::string& varName);

    std::vector<Insn> prog_;
    size_t depth_ = 0;
    bool ok_ = false;
    std::string err_;
};

}  // namespace grav
