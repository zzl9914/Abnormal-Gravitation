#include "expr.hpp"

#include "gsstr.h"

namespace grav {

namespace {

struct Tok {
    enum Kind { Num, Name, Op, LP, RP, End } kind;
    double num = 0.0;
    int pos = 0;
    int len = 0;  // 文本长度（配合源串定位），避免保存字符串
};

// ---------------------------------------------------------------- tokenizer

struct Lexer {
    const char* s;
    Tok toks[kMaxTok];
    int n = 0;
    char err[kMaxErr];

    bool fail(const char* msg) {
        err[0] = '\0';
        gs_strcat(err, kMaxErr, msg);
        return false;
    }
    bool room() {
        if (n + 2 > kMaxTok) return fail("表达式过长");
        return true;
    }

    bool run() {
        int i = 0;
        while (s[i]) {
            char c = s[i];
            if (gs_isspace(c)) {
                i++;
                continue;
            }
            if (!room()) return false;
            if (gs_isdigit(c) ||
                (c == '.' && gs_isdigit(s[i + 1]))) {
                const char* start = s + i;
                const char* end = nullptr;
                double v = gs_strtod(start, &end);
                if (end == start) return fail("无法解析的数字");
                if (!gs_isfinite(v)) return fail("数字无效");
                Tok& t = toks[n++];
                t.kind = Tok::Num;
                t.num = v;
                t.pos = i;
                t.len = (int)(end - start);
                i += t.len;
                continue;
            }
            if (gs_isalpha(c)) {
                int j = i;
                while (s[j] && gs_isalnum(s[j])) j++;
                Tok& t = toks[n++];
                t.kind = Tok::Name;
                t.pos = i;
                t.len = j - i;
                i = j;
                continue;
            }
            if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^' ||
                c == '(' || c == ')') {
                Tok& t = toks[n++];
                t.kind = (c == '(') ? Tok::LP : (c == ')') ? Tok::RP : Tok::Op;
                t.pos = i;
                t.len = 1;
                i++;
                continue;
            }
            err[0] = '\0';
            gs_strcat(err, kMaxErr, "无法识别的字符 '");
            gs_strcat_char(err, kMaxErr, c);
            gs_strcat(err, kMaxErr, "'");
            return false;
        }
        Tok& end = toks[n++];
        end.kind = Tok::End;
        end.pos = i;
        end.len = 0;
        return true;
    }
};

// ------------------------------------------------------------------ parser
// 递归下降，**直接发射 RPN**（freestanding 下没有 new，不再构造 AST）。

struct Parser {
    const Tok* toks;
    const char* src;
    const char* varName;
    Insn* prog;
    int progLen = 0;
    int depth = 0;
    int maxDepth = 0;
    char err[kMaxErr];

    const Tok& peek() const { return toks[pos]; }
    int pos = 0;

    bool isOp(char c) const {
        const Tok& t = toks[pos];
        return t.kind == Tok::Op && src[t.pos] == c;
    }
    bool tokIs(const Tok& t, const char* lit) const {
        int L = gs_strlen(lit);
        if (t.len != L) return false;
        for (int i = 0; i < L; i++)
            if (src[t.pos + i] != lit[i]) return false;
        return true;
    }
    bool fail(const char* msg) {
        err[0] = '\0';
        gs_strcat(err, kMaxErr, msg);
        return false;
    }
    bool room(int k) {
        if (progLen + k > kMaxInsn) return fail("表达式过长");
        return true;
    }

    bool emitPush(Op op, double num = 0.0) {
        if (!room(1)) return false;
        prog[progLen].op = op;
        prog[progLen].num = num;
        progLen++;
        depth++;
        if (depth > maxDepth) maxDepth = depth;
        return true;
    }
    bool emitNet0(Op op, double num = 0.0) {  // 消耗 1 压入 1
        if (!room(1)) return false;
        prog[progLen].op = op;
        prog[progLen].num = num;
        progLen++;
        return true;
    }
    bool emitBin(Op op) {  // 消耗 2 压入 1
        if (!room(1)) return false;
        prog[progLen].op = op;
        prog[progLen].num = 0.0;
        progLen++;
        depth--;
        return true;
    }

    bool parseExpr() {
        if (!parseTerm()) return false;
        while (isOp('+') || isOp('-')) {
            char c = src[toks[pos].pos];
            pos++;
            if (!parseTerm()) return false;
            if (!emitBin(c == '+' ? Op::Add : Op::Sub)) return false;
        }
        return true;
    }

    bool parseTerm() {
        if (!parseUnary()) return false;
        while (isOp('*') || isOp('/')) {
            char c = src[toks[pos].pos];
            pos++;
            if (!parseUnary()) return false;
            if (!emitBin(c == '*' ? Op::Mul : Op::Div)) return false;
        }
        return true;
    }

    bool parseUnary() {
        if (isOp('-')) {
            pos++;
            int mark = progLen;
            if (!parseUnary()) return false;
            // 常量折叠：-2 → PushConst(-2)。这样 r^-2 也能折叠成 PowC，
            // 从而让 F = 1/r^2 之类的表达式走**解析反积分**而不是数值表。
            if (progLen == mark + 1 && prog[mark].op == Op::PushConst) {
                prog[mark].num = -prog[mark].num;
                return true;
            }
            return emitNet0(Op::Neg);
        }
        if (isOp('+')) {
            pos++;
            return parseUnary();
        }
        return parsePower();
    }

    bool parsePower() {
        if (!parseAtom()) return false;
        if (isOp('^')) {
            pos++;
            int mark = progLen;
            if (!parseUnary()) return false;  // 右结合；允许 r^-2
            if (progLen == mark + 1 && prog[mark].op == Op::PushConst) {
                double c = prog[mark].num;
                progLen = mark;  // 弹出常数，换成 PowC
                // 深度修正：[常数, PowGen] 净 −1；[PowC] 净 0
                depth--;
                return emitNet0(Op::PowC, c);
            }
            return emitBin(Op::PowGen);
        }
        return true;
    }

    bool parseAtom() {
        const Tok& t = peek();
        if (t.kind == Tok::Num) {
            pos++;
            return emitPush(Op::PushConst, t.num);
        }
        if (t.kind == Tok::LP) {
            pos++;
            if (!parseExpr()) return false;
            if (peek().kind != Tok::RP) return fail("缺少右括号");
            pos++;
            return true;
        }
        if (t.kind == Tok::Name) {
            if (tokIs(t, varName)) {
                pos++;
                return emitPush(Op::PushVar);
            }
            if (tokIs(t, "pi")) {
                pos++;
                return emitPush(Op::PushConst, GS_PI);
            }
            if (tokIs(t, "e")) {
                pos++;
                return emitPush(Op::PushConst, 2.71828182845904523536);
            }
            static const char* kFns[] = {"sqrt", "exp", "ln",  "log", "abs",
                                         "sgn",  "sin", "cos", "tan", "atan"};
            static const Op kOps[] = {Op::Sqrt, Op::Exp, Op::Ln,  Op::Ln, Op::Abs,
                                      Op::Sgn,  Op::Sin, Op::Cos, Op::Tan, Op::Atan};
            for (int i = 0; i < 10; i++) {
                if (!tokIs(t, kFns[i])) continue;
                pos++;
                if (peek().kind != Tok::LP) {
                    err[0] = '\0';
                    gs_strcat(err, kMaxErr, "函数 ");
                    gs_strcat(err, kMaxErr, kFns[i]);
                    gs_strcat(err, kMaxErr, " 后需要 '('");
                    return false;
                }
                pos++;
                if (!parseExpr()) return false;
                if (peek().kind != Tok::RP) {
                    err[0] = '\0';
                    gs_strcat(err, kMaxErr, "函数 ");
                    gs_strcat(err, kMaxErr, kFns[i]);
                    gs_strcat(err, kMaxErr, " 缺少右括号");
                    return false;
                }
                pos++;
                return emitNet0(kOps[i]);
            }
            err[0] = '\0';
            gs_strcat(err, kMaxErr, "未知名称 '");
            for (int i = 0; i < t.len && i < 24; i++)
                gs_strcat_char(err, kMaxErr, src[t.pos + i]);
            gs_strcat(err, kMaxErr, "'（变量应为 ");
            gs_strcat(err, kMaxErr, varName);
            gs_strcat(err, kMaxErr, "）");
            return false;
        }
        return fail("此处应为数字、变量、函数或 '('");
    }
};

// ------------------------------------------------------------ 幂和符号栈

struct PolySum {
    double p[kMaxPowTerms];
    double a[kMaxPowTerms];
    int n = 0;
    bool over = false;  // 项数溢出：调用方按失败处理（走数值表路径）

    void clear() {
        n = 0;
        over = false;
    }
    int find(double pp) const {
        for (int i = 0; i < n; i++)
            if (p[i] == pp) return i;
        return -1;
    }
    void insert(double pp, double aa) {
        if (aa == 0.0) return;
        int k = find(pp);
        if (k >= 0) {
            a[k] += aa;
            if (a[k] == 0.0) {  // 完全抵消：删除该幂
                for (int j = k; j + 1 < n; j++) {
                    p[j] = p[j + 1];
                    a[j] = a[j + 1];
                }
                n--;
            }
            return;
        }
        if (n >= kMaxPowTerms) {
            over = true;
            return;
        }
        p[n] = pp;
        a[n] = aa;
        n++;
    }
    void add(const PolySum& o, double sign) {
        for (int i = 0; i < o.n; i++) insert(o.p[i], sign * o.a[i]);
        if (o.over) over = true;
    }
    void negate() {
        for (int i = 0; i < n; i++) a[i] = -a[i];
    }
    bool mul(const PolySum& o) {  // 笛卡尔积，指数相加系数相乘
        PolySum r;
        r.clear();
        for (int i = 0; i < n; i++)
            for (int j = 0; j < o.n; j++) {
                r.insert(p[i] + o.p[j], a[i] * o.a[j]);
                if (r.over) return false;
            }
        *this = r;
        return true;
    }
    bool div(const PolySum& o) {  // 分母必须是单项
        if (o.n != 1) return false;
        if (o.a[0] == 0.0) return false;
        for (int i = 0; i < n; i++) {
            p[i] -= o.p[0];
            a[i] /= o.a[0];
        }
        return true;
    }
    bool powC(double c) {  // 只允许单项：(a r^p)^c = a^c r^{pc}
        if (n == 0) {
            if (c == 0.0) insert(0.0, 1.0);
            return c >= 0.0;
        }
        if (n != 1) return false;
        if (a[0] == 0.0) return false;
        a[0] = gs_pow(a[0], c);
        p[0] *= c;
        return true;
    }
};

}  // namespace

// ------------------------------------------------------------------ compile

bool Expr::compile(const char* text, const char* varName) {
    progLen_ = 0;
    maxDepth_ = 0;
    ok_ = false;
    err_[0] = '\0';
    if (!text || !text[0]) {
        gs_strcat(err_, kMaxErr, "表达式为空");
        return false;
    }
    Lexer lx;
    lx.s = text;
    lx.n = 0;
    lx.err[0] = '\0';
    if (!lx.run()) {
        gs_strcat(err_, kMaxErr, lx.err);
        return false;
    }
    Parser p;
    p.toks = lx.toks;
    p.src = text;
    p.varName = varName;
    p.prog = prog_;
    p.pos = 0;
    p.err[0] = '\0';
    if (!p.parseExpr()) {
        gs_strcat(err_, kMaxErr, p.err);
        return false;
    }
    if (p.peek().kind != Tok::End) {
        gs_strcat(err_, kMaxErr, "表达式末尾有多余内容");
        return false;
    }
    if (p.depth != 1 || p.maxDepth == 0) {
        gs_strcat(err_, kMaxErr, "表达式结构无效");
        return false;
    }
    if (p.maxDepth > kMaxExprStack - 8) {
        gs_strcat(err_, kMaxErr, "表达式嵌套过深");
        return false;
    }
    progLen_ = p.progLen;
    maxDepth_ = p.maxDepth;
    ok_ = true;
    return true;
}

// --------------------------------------------------------------------- eval

// ⚠ 整数指数快路径（实测动机）：wasm 侧 gs_pow 是 JS import（Math.pow），
//   指数不是 ±1/±2 这类小整数时会掉进 fdlibm——N=8192 实测 1/r^12 每步 257 ms，
//   而 1/r^2 只要 55 ms（同样是一次幂运算），LJ（两个幂项）更是 437 ms。
//   改成乘法平方链后：① 没有 import 调用；② 力循环后续可做 SIMD128 向量化。
static inline bool integralExp(double c, int& n){
    // 不用 floor（gsmath 没有，也不想为一个判断再加一个 wasm import）：
    // |c| ≤ 64 时 int↔double 往返是精确的，直接比 (double)(int)c == c
    if (c < -64.0 || c > 64.0) return false;
    int k = (int)c;
    if ((double)k != c) return false;
    n = k;
    return true;
}
static inline double ipowChain(double x, int n){
    unsigned m = (n < 0) ? (unsigned)(-n) : (unsigned)n;
    double r = 1.0, b = x;
    while (m) { if (m & 1u) r *= b; b *= b; m >>= 1; }
    return (n < 0) ? 1.0 / r : r;
}

Dual Expr::eval(double x) const {
    Dual st[kMaxExprStack];
    int sp = 0;
    const double kNaN = __builtin_nan("");
    for (int i = 0; i < progLen_; i++) {
        const Insn& in = prog_[i];
        switch (in.op) {
            case Op::PushVar: st[sp].v = x; st[sp].d = 1.0; sp++; break;
            case Op::PushConst: st[sp].v = in.num; st[sp].d = 0.0; sp++; break;
            case Op::Add: {
                Dual b = st[sp - 1], a = st[sp - 2];
                st[sp - 2].v = a.v + b.v;
                st[sp - 2].d = a.d + b.d;
                sp--;
                break;
            }
            case Op::Sub: {
                Dual b = st[sp - 1], a = st[sp - 2];
                st[sp - 2].v = a.v - b.v;
                st[sp - 2].d = a.d - b.d;
                sp--;
                break;
            }
            case Op::Mul: {
                Dual b = st[sp - 1], a = st[sp - 2];
                st[sp - 2].v = a.v * b.v;
                st[sp - 2].d = a.d * b.v + a.v * b.d;
                sp--;
                break;
            }
            case Op::Div: {
                Dual b = st[sp - 1], a = st[sp - 2];
                if (b.v == 0.0) {
                    if (a.v == 0.0) {
                        st[sp - 2].v = kNaN;
                        st[sp - 2].d = kNaN;
                    } else {
                        st[sp - 2].v = gs_copysign(__builtin_inf(), a.v);
                        st[sp - 2].d = (a.d == 0.0) ? 0.0
                                                    : gs_copysign(__builtin_inf(), a.d);
                    }
                } else {
                    st[sp - 2].v = a.v / b.v;
                    st[sp - 2].d = (a.d * b.v - a.v * b.d) / (b.v * b.v);
                }
                sp--;
                break;
            }
            case Op::Neg: {
                st[sp - 1].v = -st[sp - 1].v;
                st[sp - 1].d = -st[sp - 1].d;
                break;
            }
            case Op::PowC: {
                Dual a = st[sp - 1];
                const double c = in.num;
                double val, dvdx;                 // dvdx = c·x^(c−1)
                int n;
                if (integralExp(c, n)) {          // 整数指数：乘法链，不调 gs_pow
                    val = ipowChain(a.v, n);
                    dvdx = c * ipowChain(a.v, n - 1);
                } else {
                    val = gs_pow(a.v, c);
                    dvdx = c * gs_pow(a.v, c - 1.0);
                }
                double der;
                if (a.v != 0.0) {
                    der = dvdx * a.d;
                } else if (c == 1.0) {
                    der = a.d;
                } else if (c > 1.0) {
                    der = 0.0;
                } else {
                    der = (a.d == 0.0) ? 0.0 : gs_copysign(__builtin_inf(), a.d);
                }
                st[sp - 1].v = val;
                st[sp - 1].d = der;
                break;
            }
            case Op::PowGen: {
                Dual w = st[sp - 1], a = st[sp - 2];
                if (a.v > 0.0) {
                    double val = gs_pow(a.v, w.v);
                    st[sp - 2].v = val;
                    st[sp - 2].d = val * (w.d * gs_log(a.v) + w.v * a.d / a.v);
                } else if (a.v == 0.0 && w.v == 0.0) {
                    st[sp - 2].v = 1.0;
                    st[sp - 2].d = 0.0;
                } else {
                    st[sp - 2].v = kNaN;
                    st[sp - 2].d = kNaN;
                }
                sp--;
                break;
            }
            case Op::Sqrt: {
                Dual a = st[sp - 1];
                if (a.v > 0.0) {
                    double s = gs_sqrt(a.v);
                    st[sp - 1].v = s;
                    st[sp - 1].d = a.d / (2.0 * s);
                } else if (a.v == 0.0) {
                    st[sp - 1].v = 0.0;
                    st[sp - 1].d =
                        (a.d == 0.0) ? 0.0 : gs_copysign(__builtin_inf(), a.d);
                } else {
                    st[sp - 1].v = kNaN;
                    st[sp - 1].d = kNaN;
                }
                break;
            }
            case Op::Exp: {
                Dual a = st[sp - 1];
                double ev = gs_exp(a.v);
                st[sp - 1].v = ev;
                st[sp - 1].d = ev * a.d;
                break;
            }
            case Op::Ln: {
                Dual a = st[sp - 1];
                if (a.v > 0.0) {
                    st[sp - 1].v = gs_log(a.v);
                    st[sp - 1].d = a.d / a.v;
                } else {
                    st[sp - 1].v = kNaN;
                    st[sp - 1].d = kNaN;
                }
                break;
            }
            case Op::Abs: {
                Dual a = st[sp - 1];
                double s = (a.v > 0.0) - (a.v < 0.0);
                st[sp - 1].v = gs_fabs(a.v);
                st[sp - 1].d = a.d * s;
                break;
            }
            case Op::Sgn: {
                Dual a = st[sp - 1];
                st[sp - 1].v = (a.v > 0.0) - (a.v < 0.0);
                st[sp - 1].d = 0.0;
                break;
            }
            case Op::Sin: {
                Dual a = st[sp - 1];
                st[sp - 1].v = gs_sin(a.v);
                st[sp - 1].d = gs_cos(a.v) * a.d;
                break;
            }
            case Op::Cos: {
                Dual a = st[sp - 1];
                st[sp - 1].v = gs_cos(a.v);
                st[sp - 1].d = -gs_sin(a.v) * a.d;
                break;
            }
            case Op::Tan: {
                Dual a = st[sp - 1];
                double c = gs_cos(a.v);
                if (gs_fabs(c) < 1e-300) {
                    st[sp - 1].v = kNaN;
                    st[sp - 1].d = kNaN;
                } else {
                    st[sp - 1].v = gs_tan(a.v);
                    st[sp - 1].d = a.d / (c * c);
                }
                break;
            }
            case Op::Atan: {
                Dual a = st[sp - 1];
                st[sp - 1].v = gs_atan(a.v);
                st[sp - 1].d = a.d / (1.0 + a.v * a.v);
                break;
            }
        }
    }
    return st[0];
}

// 只求值的走法：**不算导数**。力循环每一步都要调 F(r)，而 eval() 里一半算术
// 是在算用不到的 .d（对偶数导数）。语义与 eval(x).v 完全一致（除了整数指数走
// 同一条乘法链快路径），只有舍入次序的差别。
double Expr::valueOnly(double x) const {
    double st[kMaxExprStack];
    int sp = 0;
    const double kNaN = __builtin_nan("");
    for (int i = 0; i < progLen_; i++) {
        const Insn& in = prog_[i];
        switch (in.op) {
            case Op::PushVar:   st[sp++] = x; break;
            case Op::PushConst: st[sp++] = in.num; break;
            case Op::Neg:       st[sp - 1] = -st[sp - 1]; break;
            case Op::Add: { double b = st[sp - 1], a = st[sp - 2]; st[sp - 2] = a + b; sp--; break; }
            case Op::Sub: { double b = st[sp - 1], a = st[sp - 2]; st[sp - 2] = a - b; sp--; break; }
            case Op::Mul: { double b = st[sp - 1], a = st[sp - 2]; st[sp - 2] = a * b; sp--; break; }
            case Op::Div: { double b = st[sp - 1], a = st[sp - 2]; st[sp - 2] = a / b; sp--; break; }
            case Op::PowC: {
                int n; double c = in.num;
                st[sp - 1] = integralExp(c, n) ? ipowChain(st[sp - 1], n) : gs_pow(st[sp - 1], c);
                break;
            }
            case Op::PowGen: {
                double w = st[sp - 1], a = st[sp - 2];
                if (a > 0.0) st[sp - 2] = gs_pow(a, w);
                else if (a == 0.0 && w == 0.0) st[sp - 2] = 1.0;
                else st[sp - 2] = kNaN;
                sp--;
                break;
            }
            case Op::Sqrt: {
                double a = st[sp - 1];
                st[sp - 1] = (a > 0.0) ? gs_sqrt(a) : (a == 0.0 ? 0.0 : kNaN);
                break;
            }
            case Op::Exp:  st[sp - 1] = gs_exp(st[sp - 1]); break;
            case Op::Ln: {
                double a = st[sp - 1];
                st[sp - 1] = (a > 0.0) ? gs_log(a) : kNaN;
                break;
            }
            case Op::Abs:  st[sp - 1] = gs_fabs(st[sp - 1]); break;
            case Op::Sgn:  st[sp - 1] = (double)((st[sp - 1] > 0.0) - (st[sp - 1] < 0.0)); break;
            case Op::Sin:  st[sp - 1] = gs_sin(st[sp - 1]); break;
            case Op::Cos:  st[sp - 1] = gs_cos(st[sp - 1]); break;
            case Op::Tan: {
                double c = gs_cos(st[sp - 1]);
                st[sp - 1] = (gs_fabs(c) < 1e-300) ? kNaN : gs_tan(st[sp - 1]);
                break;
            }
            case Op::Atan: st[sp - 1] = gs_atan(st[sp - 1]); break;
        }
    }
    return st[0];
}

bool Expr::asPowerOfVar(double& p) const {
    if (!ok_) return false;
    if (progLen_ == 1 && prog_[0].op == Op::PushVar) {
        p = 1.0;
        return true;
    }
    // 编译器对「指数为字面常数」发射 [PushVar, PowC(指数)]
    if (progLen_ == 2 && prog_[0].op == Op::PushVar && prog_[1].op == Op::PowC) {
        p = prog_[1].num;
        return p != 0.0;
    }
    return false;
}

bool Expr::asPowerSum(double* outP, double* outA, int cap, int& n) const {
    n = 0;
    if (!ok_) return false;
    PolySum st[kMaxExprStack];
    int sp = 0;
    for (int i = 0; i < progLen_; i++) {
        const Insn& in = prog_[i];
        switch (in.op) {
            case Op::PushVar: {
                st[sp].clear();
                st[sp].insert(1.0, 1.0);
                sp++;
                break;
            }
            case Op::PushConst: {
                st[sp].clear();
                st[sp].insert(0.0, in.num);
                sp++;
                break;
            }
            case Op::Neg: {
                st[sp - 1].negate();
                break;
            }
            case Op::Add: {
                PolySum b = st[sp - 1];
                st[sp - 2].add(b, +1.0);
                sp--;
                break;
            }
            case Op::Sub: {
                PolySum b = st[sp - 1];
                st[sp - 2].add(b, -1.0);
                sp--;
                break;
            }
            case Op::Mul: {
                PolySum b = st[sp - 1];
                if (!st[sp - 2].mul(b)) return false;
                sp--;
                break;
            }
            case Op::Div: {
                PolySum b = st[sp - 1];
                if (!st[sp - 2].div(b)) return false;
                sp--;
                break;
            }
            case Op::PowC: {
                if (!st[sp - 1].powC(in.num)) return false;
                break;
            }
            case Op::Sqrt: {
                if (!st[sp - 1].powC(0.5)) return false;
                break;
            }
            case Op::Abs: {
                PolySum& a = st[sp - 1];
                if (a.n > 1) return false;
                if (a.n == 1 && a.a[0] < 0.0) a.a[0] = -a.a[0];
                break;
            }
            default:  // PowGen, Exp, Ln, Sgn, Sin, Cos, Tan, Atan
                return false;
        }
    }
    if (sp != 1) return false;
    PolySum& r = st[0];
    if (r.n > cap) return false;
    for (int i = 0; i < r.n; i++) {
        if (!gs_isfinite(r.p[i]) || !gs_isfinite(r.a[i])) return false;
        outP[i] = r.p[i];
        outA[i] = r.a[i];
    }
    n = r.n;
    return true;
}

bool Expr::triviallyConstant() const {
    return ok_ && progLen_ == 1 && prog_[0].op == Op::PushConst;
}

}  // namespace grav
