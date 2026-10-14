#include "expr.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>

namespace grav {

namespace {

// ---------------------------------------------------------------- tokenizer

struct Tok {
    enum Kind { Num, Name, Op, LP, RP, End } kind;
    double num = 0.0;
    std::string text;
    int pos = 0;
};

std::vector<Tok> lex(const std::string& s, std::string& err) {
    std::vector<Tok> out;
    size_t i = 0;
    const size_t n = s.size();
    while (i < n) {
        unsigned char c = (unsigned char)s[i];
        if (std::isspace(c)) { i++; continue; }
        bool digit = std::isdigit(c) != 0;
        bool dot = c == '.' && i + 1 < n && std::isdigit((unsigned char)s[i + 1]);
        if (digit || dot) {
            size_t j = i;
            while (j < n && (std::isdigit((unsigned char)s[j]) || s[j] == '.')) j++;
            if (j < n && (s[j] == 'e' || s[j] == 'E')) {
                size_t k = j + 1;
                if (k < n && (s[k] == '+' || s[k] == '-')) k++;
                if (k < n && std::isdigit((unsigned char)s[k])) {
                    while (k < n && std::isdigit((unsigned char)s[k])) k++;
                    j = k;
                }
            }
            Tok t;
            t.kind = Tok::Num;
            t.text = s.substr(i, j - i);
            t.num = std::strtod(t.text.c_str(), nullptr);
            t.pos = (int)i;
            if (!std::isfinite(t.num)) {
                err = "数字无效: " + t.text;
                return {};
            }
            out.push_back(t);
            i = j;
            continue;
        }
        if (std::isalpha(c) || c == '_') {
            size_t j = i;
            while (j < n && (std::isalnum((unsigned char)s[j]) || s[j] == '_')) j++;
            Tok t;
            t.kind = Tok::Name;
            t.text = s.substr(i, j - i);
            t.pos = (int)i;
            out.push_back(t);
            i = j;
            continue;
        }
        if (std::strchr("+-*/^()", c)) {
            Tok t;
            t.kind = (c == '(') ? Tok::LP : (c == ')') ? Tok::RP : Tok::Op;
            t.text.assign(1, (char)c);
            t.pos = (int)i;
            out.push_back(t);
            i++;
            continue;
        }
        err = std::string("无法识别的字符 '") + (char)c + "'";
        return {};
    }
    Tok end;
    end.kind = Tok::End;
    out.push_back(end);
    return out;
}

// ---------------------------------------------------------------------- AST

struct Node {
    enum Kind { Var, Num, Neg, Bin, Fn } kind;
    double num = 0.0;
    char op = 0;
    std::string fn;
    std::unique_ptr<Node> a, b;
};

using NodeP = std::unique_ptr<Node>;

NodeP makeNum(double v) {
    NodeP p(new Node);
    p->kind = Node::Num;
    p->num = v;
    return p;
}

NodeP makeVar() {
    NodeP p(new Node);
    p->kind = Node::Var;
    return p;
}

struct Parser {
    const std::vector<Tok>& toks;
    const std::string& varName;
    std::string err;
    size_t pos = 0;

    const Tok& peek() const { return toks[pos]; }
    const Tok& next() { return toks[pos++]; }

    NodeP parseExpr() {
        NodeP left = parseTerm();
        if (!left) return nullptr;
        while (peek().kind == Tok::Op && (peek().text == "+" || peek().text == "-")) {
            char op = next().text[0];
            NodeP right = parseTerm();
            if (!right) return nullptr;
            NodeP p(new Node);
            p->kind = Node::Bin;
            p->op = op;
            p->a = std::move(left);
            p->b = std::move(right);
            left = std::move(p);
        }
        return left;
    }

    NodeP parseTerm() {
        NodeP left = parseUnary();
        if (!left) return nullptr;
        while (peek().kind == Tok::Op && (peek().text == "*" || peek().text == "/")) {
            char op = next().text[0];
            NodeP right = parseUnary();
            if (!right) return nullptr;
            NodeP p(new Node);
            p->kind = Node::Bin;
            p->op = op;
            p->a = std::move(left);
            p->b = std::move(right);
            left = std::move(p);
        }
        return left;
    }

    NodeP parseUnary() {
        if (peek().kind == Tok::Op && peek().text == "-") {
            next();
            NodeP inner = parseUnary();
            if (!inner) return nullptr;
            if (inner->kind == Node::Num) return makeNum(-inner->num);
            NodeP p(new Node);
            p->kind = Node::Neg;
            p->a = std::move(inner);
            return p;
        }
        if (peek().kind == Tok::Op && peek().text == "+") {
            next();
            return parseUnary();
        }
        return parsePower();
    }

    NodeP parsePower() {
        NodeP base = parseAtom();
        if (!base) return nullptr;
        if (peek().kind == Tok::Op && peek().text == "^") {
            next();
            NodeP expo = parseUnary();  // right associative; allows r^-2
            if (!expo) return nullptr;
            NodeP p(new Node);
            p->kind = Node::Bin;
            p->op = '^';
            p->a = std::move(base);
            p->b = std::move(expo);
            return p;
        }
        return base;
    }

    NodeP parseAtom() {
        const Tok& t = peek();
        if (t.kind == Tok::Num) {
            next();
            return makeNum(t.num);
        }
        if (t.kind == Tok::LP) {
            next();
            NodeP e = parseExpr();
            if (!e) return nullptr;
            if (peek().kind != Tok::RP) {
                err = "缺少右括号";
                return nullptr;
            }
            next();
            return e;
        }
        if (t.kind == Tok::Name) {
            next();
            const std::string& name = t.text;
            if (name == varName) return makeVar();
            if (name == "pi") return makeNum(3.14159265358979323846);
            if (name == "e") return makeNum(2.71828182845904523536);
            static const char* kFns[] = {"sqrt", "exp", "ln", "log", "abs", "sgn",
                                         "sin", "cos", "tan", "atan"};
            for (const char* f : kFns) {
                if (name == f) {
                    if (peek().kind != Tok::LP) {
                        err = "函数 " + name + " 后需要 '('";
                        return nullptr;
                    }
                    next();
                    NodeP arg = parseExpr();
                    if (!arg) return nullptr;
                    if (peek().kind != Tok::RP) {
                        err = "函数 " + name + " 缺少右括号";
                        return nullptr;
                    }
                    next();
                    NodeP p(new Node);
                    p->kind = Node::Fn;
                    p->fn = name;
                    p->a = std::move(arg);
                    return p;
                }
            }
            err = "未知名称 '" + name + "'（变量应为 " + varName + "）";
            return nullptr;
        }
        err = "此处应为数字、变量、函数或 '('";
        return nullptr;
    }
};

// ----------------------------------------------------------------- compiler
// `depth` = number of values currently on the evaluation stack.
// Every compiled sub-expression leaves exactly one value on the stack.

struct Compiler {
    std::vector<Insn> prog;
    size_t depth = 0;
    size_t maxDepth = 0;

    void push1(Op op, double num = 0.0) {  // net +1
        Insn in;
        in.op = op;
        in.num = num;
        prog.push_back(in);
        depth++;
        if (depth > maxDepth) maxDepth = depth;
    }
    void net0(Op op, double num = 0.0) {  // consumes 1, pushes 1
        Insn in;
        in.op = op;
        in.num = num;
        prog.push_back(in);
    }
    void binary(Op op) {  // consumes 2, pushes 1
        Insn in;
        in.op = op;
        prog.push_back(in);
        depth--;
    }

    void run(const Node* nd) {
        switch (nd->kind) {
            case Node::Var:
                push1(Op::PushVar);
                break;
            case Node::Num:
                push1(Op::PushConst, nd->num);
                break;
            case Node::Neg:
                run(nd->a.get());
                net0(Op::Neg);
                break;
            case Node::Bin:
                if (nd->op == '^' && nd->b->kind == Node::Num) {
                    run(nd->a.get());
                    net0(Op::PowC, nd->b->num);
                } else {
                    run(nd->a.get());
                    run(nd->b.get());
                    switch (nd->op) {
                        case '+': binary(Op::Add); break;
                        case '-': binary(Op::Sub); break;
                        case '*': binary(Op::Mul); break;
                        case '/': binary(Op::Div); break;
                        case '^': binary(Op::PowGen); break;
                    }
                }
                break;
            case Node::Fn: {
                run(nd->a.get());
                Op op;
                if (nd->fn == "sqrt") op = Op::Sqrt;
                else if (nd->fn == "exp") op = Op::Exp;
                else if (nd->fn == "ln" || nd->fn == "log") op = Op::Ln;
                else if (nd->fn == "abs") op = Op::Abs;
                else if (nd->fn == "sgn") op = Op::Sgn;
                else if (nd->fn == "sin") op = Op::Sin;
                else if (nd->fn == "cos") op = Op::Cos;
                else if (nd->fn == "tan") op = Op::Tan;
                else op = Op::Atan;
                net0(op);
                break;
            }
        }
    }
};

}  // namespace

Expr::Expr(const std::string& text, const std::string& varName) {
    compile(text, varName);
}

bool Expr::compile(const std::string& text, const std::string& varName) {
    if (text.empty()) {
        err_ = "表达式为空";
        return false;
    }
    std::string lexErr;
    std::vector<Tok> toks = lex(text, lexErr);
    if (!lexErr.empty()) {
        err_ = lexErr;
        return false;
    }
    Parser p{toks, varName, {}, 0};
    NodeP root = p.parseExpr();
    if (!root) {
        err_ = p.err.empty() ? "解析失败" : p.err;
        return false;
    }
    if (p.peek().kind != Tok::End) {
        err_ = "表达式末尾有多余内容";
        return false;
    }
    Compiler c;
    c.run(root.get());
    if (c.depth != 1 || c.maxDepth == 0) {
        err_ = "表达式结构无效";
        return false;
    }
    if (c.maxDepth > 128) {
        err_ = "表达式嵌套过深";
        return false;
    }
    prog_ = std::move(c.prog);
    depth_ = c.maxDepth;
    ok_ = true;
    return true;
}

Dual Expr::eval(double x) const {
    Dual st[130];
    int sp = 0;
    for (const Insn& in : prog_) {
        switch (in.op) {
            case Op::PushVar: st[sp++] = Dual(x, 1.0); break;
            case Op::PushConst: st[sp++] = Dual(in.num, 0.0); break;
            case Op::Add: {
                Dual b = st[--sp], a = st[sp - 1];
                st[sp - 1] = Dual(a.v + b.v, a.d + b.d);
                break;
            }
            case Op::Sub: {
                Dual b = st[--sp], a = st[sp - 1];
                st[sp - 1] = Dual(a.v - b.v, a.d - b.d);
                break;
            }
            case Op::Mul: {
                Dual b = st[--sp], a = st[sp - 1];
                st[sp - 1] = Dual(a.v * b.v, a.d * b.v + a.v * b.d);
                break;
            }
            case Op::Div: {
                Dual b = st[--sp], a = st[sp - 1];
                if (b.v == 0.0) {
                    if (a.v == 0.0) st[sp - 1] = Dual(NAN, NAN);
                    else st[sp - 1] = Dual(std::copysign(INFINITY, a.v),
                                           a.d == 0.0 ? 0.0 : std::copysign(INFINITY, a.d));
                } else {
                    st[sp - 1] = Dual(a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v));
                }
                break;
            }
            case Op::Neg: {
                Dual& a = st[sp - 1];
                a = Dual(-a.v, -a.d);
                break;
            }
            case Op::PowC: {
                Dual a = st[sp - 1];
                const double c = in.num;
                double val = std::pow(a.v, c);
                double der;
                if (a.v != 0.0) {
                    der = c * std::pow(a.v, c - 1.0) * a.d;
                } else if (c == 1.0) {
                    der = a.d;
                } else if (c > 1.0) {
                    der = 0.0;
                } else {
                    der = (a.d == 0.0) ? 0.0 : std::copysign(INFINITY, a.d);
                }
                st[sp - 1] = Dual(val, der);
                break;
            }
            case Op::PowGen: {
                Dual w = st[--sp], a = st[sp - 1];
                if (a.v > 0.0) {
                    double val = std::pow(a.v, w.v);
                    st[sp - 1] = Dual(val, val * (w.d * std::log(a.v) + w.v * a.d / a.v));
                } else if (a.v == 0.0 && w.v == 0.0) {
                    st[sp - 1] = Dual(1.0, 0.0);
                } else {
                    st[sp - 1] = Dual(NAN, NAN);
                }
                break;
            }
            case Op::Sqrt: {
                Dual a = st[sp - 1];
                if (a.v > 0.0) {
                    double s = std::sqrt(a.v);
                    st[sp - 1] = Dual(s, a.d / (2.0 * s));
                } else if (a.v == 0.0) {
                    st[sp - 1] = Dual(0.0, a.d == 0.0 ? 0.0 : std::copysign(INFINITY, a.d));
                } else {
                    st[sp - 1] = Dual(NAN, NAN);
                }
                break;
            }
            case Op::Exp: {
                Dual a = st[sp - 1];
                double ev = std::exp(a.v);
                st[sp - 1] = Dual(ev, ev * a.d);
                break;
            }
            case Op::Ln: {
                Dual a = st[sp - 1];
                if (a.v > 0.0) st[sp - 1] = Dual(std::log(a.v), a.d / a.v);
                else st[sp - 1] = Dual(NAN, NAN);
                break;
            }
            case Op::Abs: {
                Dual a = st[sp - 1];
                double s = (a.v > 0.0) - (a.v < 0.0);
                st[sp - 1] = Dual(std::fabs(a.v), a.d * s);
                break;
            }
            case Op::Sgn: {
                Dual a = st[sp - 1];
                st[sp - 1] = Dual((a.v > 0.0) - (a.v < 0.0), 0.0);
                break;
            }
            case Op::Sin: {
                Dual a = st[sp - 1];
                st[sp - 1] = Dual(std::sin(a.v), std::cos(a.v) * a.d);
                break;
            }
            case Op::Cos: {
                Dual a = st[sp - 1];
                st[sp - 1] = Dual(std::cos(a.v), -std::sin(a.v) * a.d);
                break;
            }
            case Op::Tan: {
                Dual a = st[sp - 1];
                double c = std::cos(a.v);
                if (std::fabs(c) < 1e-300) st[sp - 1] = Dual(NAN, NAN);
                else st[sp - 1] = Dual(std::tan(a.v), a.d / (c * c));
                break;
            }
            case Op::Atan: {
                Dual a = st[sp - 1];
                st[sp - 1] = Dual(std::atan(a.v), a.d / (1.0 + a.v * a.v));
                break;
            }
        }
    }
    return st[0];
}

bool Expr::asPowerOfVar(double& p) const {
    if (!ok_) return false;
    if (prog_.size() == 1 && prog_[0].op == Op::PushVar) {
        p = 1.0;
        return true;
    }
    // Compiler emits [PushVar, PowC(exp)] (exponent folded into the insn).
    if (prog_.size() == 2 && prog_[0].op == Op::PushVar && prog_[1].op == Op::PowC) {
        p = prog_[1].num;
        return p != 0.0;
    }
    return false;
}

// ---------------------------------------------------------------- power sums

namespace {

// Symbolic "sum of a*r^p" carried through the RPN program.  Operations keep
// the representation closed whenever the result is again a power sum and
// fail (return false) otherwise.
struct PolySum {
    std::vector<std::pair<double, double>> t;  // (p, a): sum a*r^p

    static PolySum cst(double c) {
        PolySum s;
        if (c != 0.0) s.t.push_back({0.0, c});
        return s;
    }
    static PolySum var() {
        PolySum s;
        s.t.push_back({1.0, 1.0});
        return s;
    }

    void insert(double p, double a) {
        if (a == 0.0) return;
        for (size_t i = 0; i < t.size(); i++) {
            if (t[i].first == p) {
                t[i].second += a;
                if (t[i].second == 0.0) t.erase(t.begin() + i);
                return;
            }
        }
        t.push_back({p, a});
    }
    void add(const PolySum& o, double sign) {
        for (const auto& e : o.t) insert(e.first, sign * e.second);
    }
    void negate() {
        for (auto& e : t) e.second = -e.second;
    }
    bool mul(const PolySum& o, size_t cap) {  // cartesian product stays closed
        std::vector<std::pair<double, double>> r;
        for (const auto& e1 : t) {
            for (const auto& e2 : o.t) {
                double p = e1.first + e2.first, a = e1.second * e2.second;
                bool merged = false;
                for (auto& x : r)
                    if (x.first == p) { x.second += a; merged = true; break; }
                if (!merged) r.push_back({p, a});
            }
            if (r.size() > cap) return false;
        }
        t = std::move(r);
        return true;
    }
    bool div(const PolySum& o) {  // denominator must be a single term
        if (o.t.size() != 1) return false;
        double q = o.t[0].first, b = o.t[0].second;
        if (b == 0.0) return false;
        for (auto& e : t) {
            e.first -= q;
            e.second /= b;
        }
        return true;
    }
    bool powC(double c) {  // single term only: (a r^p)^c = a^c r^{pc}
        if (t.empty()) {           // zero^c
            if (c == 0.0) t.push_back({0.0, 1.0});
            return c >= 0.0;
        }
        if (t.size() != 1) return false;
        if (t[0].second == 0.0) return false;  // (0 r^p): keep as zero
        t[0].second = std::pow(t[0].second, c);
        t[0].first *= c;
        return true;
    }
};

}  // namespace

bool Expr::asPowerSum(std::vector<std::pair<double, double>>& terms) const {
    if (!ok_) return false;
    static const size_t kCap = 64;
    std::vector<PolySum> st;
    auto pop = [&]() { PolySum v = std::move(st.back()); st.pop_back(); return v; };
    for (const Insn& in : prog_) {
        switch (in.op) {
            case Op::PushVar: st.push_back(PolySum::var()); break;
            case Op::PushConst: st.push_back(PolySum::cst(in.num)); break;
            case Op::Neg: { PolySum a = pop(); a.negate(); st.push_back(std::move(a)); break; }
            case Op::Add: { PolySum b = pop(), a = pop(); a.add(b, +1.0); st.push_back(std::move(a)); break; }
            case Op::Sub: { PolySum b = pop(), a = pop(); a.add(b, -1.0); st.push_back(std::move(a)); break; }
            case Op::Mul: { PolySum b = pop(), a = pop(); if (!a.mul(b, kCap)) return false; st.push_back(std::move(a)); break; }
            case Op::Div: { PolySum b = pop(), a = pop(); if (!a.div(b)) return false; st.push_back(std::move(a)); break; }
            case Op::PowC: { PolySum a = pop(); if (!a.powC(in.num)) return false; st.push_back(std::move(a)); break; }
            case Op::Sqrt: { PolySum a = pop(); if (!a.powC(0.5)) return false; st.push_back(std::move(a)); break; }
            case Op::Abs: {
                PolySum a = pop();
                if (a.t.size() > 1) return false;
                if (a.t.size() == 1) {
                    if (a.t[0].second < 0.0) a.t[0].second = -a.t[0].second;
                }
                st.push_back(std::move(a));
                break;
            }
            default:  // PowGen, Exp, Ln, Sgn, Sin, Cos, Tan, Atan
                return false;
        }
    }
    if (st.size() != 1) return false;
    terms = st[0].t;
    for (const auto& e : terms)
        if (!std::isfinite(e.first) || !std::isfinite(e.second)) return false;
    return true;
}

bool Expr::triviallyConstant() const {
    return ok_ && prog_.size() == 1 && prog_[0].op == Op::PushConst;
}

}  // namespace grav
