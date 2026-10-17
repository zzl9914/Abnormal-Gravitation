#include "sim.hpp"

#include "gsmath.h"
#include "gsstr.h"

namespace grav {

static const double kGradClamp = 1e6;    // |d(d)/d(dx_k)| 钳位
static const double kForceClamp = 1e12;  // 力幅钳位
static const double GS_LN10 = 2.30258509299404568402;

// ------------------------------------------------- 共享静态工作缓冲
// freestanding 下没有 new；树节点 / 溢出槽位 / 接触对都用固定容量的静态数组。
// 同一时刻只有一个树在活动（computeAccel / energy / collectContacts 每次重建），
// 因此共享是安全的。

static TreeNode g_nodes[kMaxNodes];
static int g_nodeCount = 0;
static int g_root = -1;

static int g_extraBody[kMaxExtra];
static int g_extraNext[kMaxExtra];
static int g_extraUsed = 0;

static int g_cA[kMaxContacts];
static int g_cB[kMaxContacts];
static int g_cN = 0;

// 溢出链表：头插，返回新的头索引
static int extraAlloc(int body, int head) {
    if (g_extraUsed >= kMaxExtra) return head;  // 满了就丢弃（极罕见）
    int idx = g_extraUsed++;
    g_extraBody[idx] = body;
    g_extraNext[idx] = head;
    return idx;
}

// ------------------------------------------------------------------ 构造

bool Simulation::init(int dims, const char* fExpr, const char* forceExpr) {
    // 语义默认值（静态实例靠 .bss 清零，这里显式设置）
    dim_ = dims < 1 ? 1 : dims;
    collisionMode_ = Collision::Pass;
    density_ = 100.0;
    G_ = 1.0;
    treeMode_ = 0;
    theta_ = 0.6;
    t_ = 0.0;
    nbodies_ = 0;
    accValid_ = false;
    nPow_ = 0;
    vLn_ = 0.0;
    vAnalytic_ = false;
    fIsPower_ = false;
    fExp_ = 1.0;
    ok_ = false;
    err_[0] = '\0';

    if (dim_ > 8) {
        gs_strcat(err_, kMaxErr, "内核最多支持 8 维（前端只开放 2D/3D）");
        return false;
    }
    if (!fExpr || !forceExpr) {
        gs_strcat(err_, kMaxErr, "表达式为空");
        return false;
    }
    if (!f_.compile(fExpr, "t")) {
        gs_strcat(err_, kMaxErr, "f(t) 表达式错误: ");
        gs_strcat(err_, kMaxErr, f_.error());
        return false;
    }
    if (!F_.compile(forceExpr, "r")) {
        gs_strcat(err_, kMaxErr, "F(r) 表达式错误: ");
        gs_strcat(err_, kMaxErr, F_.error());
        return false;
    }
    // 力律体检：几个探测距离上必须是有限值（其余极点交给运行时钳位）
    static const double kProbe[3] = {0.5, 1.0, 2.0};
    for (int i = 0; i < 3; i++) {
        if (!gs_isfinite(F_.value(kProbe[i]))) {
            gs_strcat(err_, kMaxErr, "F(r) 在 r = ");
            gs_strcat_num(err_, kMaxErr, kProbe[i]);
            gs_strcat(err_, kMaxErr, " 处不是有限值");
            return false;
        }
    }

    double p = 0.0;
    if (f_.asPowerOfVar(p)) {
        if (!(p > 0.0) || !gs_isfinite(p)) {
            gs_strcat(err_, kMaxErr, "f(t) = t^p 要求 p > 0 才能求逆");
            return false;
        }
        fIsPower_ = true;
        fExp_ = p;
    } else {
        // 校验：f(0) == 0 且在采样网格上严格递增
        if (gs_fabs(f_.value(0.0)) > 1e-9) {
            gs_strcat(err_, kMaxErr, "f 必须满足 f(0) = 0（否则两点重合时距离不为 0）");
            return false;
        }
        double prev = 0.0;
        for (int i = 0; i < 80; i++) {
            double t = gs_exp(-6.0 * GS_LN10 + 10.0 * GS_LN10 * i / 79.0);  // 1e-6 .. 1e4
            double v = f_.value(t);
            if (!gs_isfinite(v)) {
                gs_strcat(err_, kMaxErr, "f 在 t = ");
                gs_strcat_num(err_, kMaxErr, t);
                gs_strcat(err_, kMaxErr, " 处不是有限值");
                return false;
            }
            if (!(v > prev)) {  // 严格递增
                gs_strcat(err_, kMaxErr, "f 必须在 [0, +inf) 上严格递增（在 t = ");
                gs_strcat_num(err_, kMaxErr, t);
                gs_strcat(err_, kMaxErr, " 附近不满足）");
                return false;
            }
            prev = v;
        }
    }

    // 派生势 V(r) = ∫F dr（V' = F）
    double tp[kMaxPowTermsSim], ta[kMaxPowTermsSim];
    int nt = 0;
    if (F_.asPowerSum(tp, ta, kMaxPowTermsSim, nt)) {
        vAnalytic_ = true;  // 精确原函数
        nPow_ = 0;
        vLn_ = 0.0;
        for (int i = 0; i < nt; i++) {
            double pp = tp[i], a = ta[i];
            if (pp == -1.0) {
                vLn_ += a;  // ∫ r^-1 = ln r
            } else if (nPow_ < kMaxPowTermsSim) {
                vPowP_[nPow_] = pp + 1.0;
                vPowA_[nPow_] = a / (pp + 1.0);
                nPow_++;
            }
        }
    } else {
        buildPotentialTable();  // 数值积分
    }
    ok_ = true;
    return true;
}

// -------------------------------------------------------- 派生势 V

// 数值兜底：V 在 [1e-6, 1e6] 的对数网格上采样，每段用 4 点 Gauss-Legendre
// 在 ln 空间积分（dr = r·dln r）；插值用 ln r 上的三次 Hermite，节点斜率取
// 精确值 dV/dln r = F(r)·r；网格外用端点斜率线性外推。
void Simulation::buildPotentialTable() {
    const int N = kTabN;
    static const double xg[4] = {-0.8611363115940526, -0.3399810435848563,
                                 0.3399810435848563, 0.8611363115940526};
    static const double wg[4] = {0.34785484513745385, 0.6521451548625461,
                                 0.6521451548625461, 0.34785484513745385};
    const double lmin = -6.0 * GS_LN10;  // ln(1e-6)
    const double dl = (12.0 * GS_LN10) / (N - 1);
    double V = 0.0;
    for (int i = 0; i < N; i++) {
        double l = lmin + dl * i;
        double r = gs_exp(l);
        double fr = F_.value(r);
        vTabR_[i] = r;
        vTabV_[i] = V;
        vTabF_[i] = gs_isfinite(fr) ? fr : 0.0;
        if (i + 1 < N) {
            double lm = l + 0.5 * dl;
            double seg = 0.0;
            for (int g = 0; g < 4; g++) {
                double lg = lm + 0.5 * dl * xg[g];
                double rg = gs_exp(lg);
                double fg = F_.value(rg);
                if (gs_isfinite(fg)) seg += wg[g] * 0.5 * dl * fg * rg;
            }
            V += seg;
            if (!gs_isfinite(V)) V = V > 0 ? 1e300 : -1e300;
            if (gs_fabs(V) > 1e300) V = V > 0 ? 1e300 : -1e300;
        }
    }
}

double Simulation::potentialAt(double r) const {
    if (!(r > 0.0) || !gs_isfinite(r)) return 0.0;
    double v;
    if (vAnalytic_) {
        v = vLn_ * gs_log(r);
        for (int i = 0; i < nPow_; i++) v += vPowA_[i] * gs_pow(r, vPowP_[i]);
    } else {
        const int N = kTabN;
        const double lmin = gs_log(vTabR_[0]);
        const double dl = gs_log(vTabR_[1]) - lmin;
        const double lr = gs_log(r);
        if (lr <= lmin) {  // 网格下方：按首端斜率线性外推
            v = vTabV_[0] + vTabF_[0] * (r - vTabR_[0]);
        } else {
            const double lmax = lmin + dl * (N - 1);
            if (lr >= lmax) {  // 网格上方：按末端斜率线性外推
                v = vTabV_[N - 1] + vTabF_[N - 1] * (r - vTabR_[N - 1]);
            } else {
                double fi = (lr - lmin) / dl;
                int i = (int)fi;
                if (i > N - 2) i = N - 2;
                double u = fi - i;
                // 三次 Hermite，斜率 dV/dln r = F(r)·r
                double s0 = vTabF_[i] * vTabR_[i];
                double s1 = vTabF_[i + 1] * vTabR_[i + 1];
                double h00 = (1 + 2 * u) * (1 - u) * (1 - u);
                double h10 = u * (1 - u) * (1 - u);
                double h01 = u * u * (3 - 2 * u);
                double h11 = u * u * (u - 1);
                v = h00 * vTabV_[i] + h10 * dl * s0 + h01 * vTabV_[i + 1] +
                    h11 * dl * s1;
            }
        }
    }
    v = G_ * v;  // 耦合常数同样缩放派生势
    return gs_isfinite(v) ? v : 0.0;
}

// 度量距离 d 处的对势，含核内正则（R = R_i + R_j）
double Simulation::pairV(double d, double R) const {
    if (R > 0.0 && d < R) {
        double FR = G_ * F_.value(R);
        if (!gs_isfinite(FR)) FR = 0.0;
        return potentialAt(R) + FR * (d * d - R * R) / (2.0 * R);
    }
    return potentialAt(d);
}

void Simulation::setCollision(Collision mode, double density) {
    collisionMode_ = mode;
    density_ = density > 1e-12 ? density : 1e-12;
    for (int i = 0; i < nbodies_; i++)  // 密度变了半径要重算
        bodies_[i].radius = radiusFromMass(bodies_[i].mass);
}

void Simulation::setG(double G) {
    if (G > 0.0 && gs_isfinite(G)) G_ = G;
}

double Simulation::radiusFromMass(double m) const {
    // R = ( m·Γ(d/2+1) / (ρ·π^(d/2)) )^(1/d)
    if (!(m > 0.0)) return 0.0;
    double g = gs_tgamma(dim_ / 2.0 + 1.0);
    return gs_pow(m * g / (density_ * gs_pow(GS_PI, dim_ / 2.0)), 1.0 / dim_);
}

// ---------------------------------------------------------- Barnes-Hut 树

void Simulation::setTree(int treeMode, double theta) {
    treeMode_ = (treeMode == 1) ? 1 : 0;
    theta_ = theta > 0.05 ? theta : 0.05;
}

bool Simulation::treeActive() const {
    if (treeMode_ != 0) return false;
    return aliveCount() >= kTreeThreshold;
}

int Simulation::newNode(int parent, const double* ctr, const double* half) const {
    if (g_nodeCount >= kMaxNodes) return -1;  // 节点耗尽：调用方退化成 bucket
    TreeNode& nd = g_nodes[g_nodeCount];
    nd.parent = parent;
    nd.mass = 0.0;
    nd.mrad = 0.0;
    nd.maxRad = 0.0;
    nd.body = -1;
    nd.extra = -1;
    nd.state = 0;
    for (int a = 0; a < 8; a++) {
        nd.child[a] = -1;
        nd.com[a] = 0.0;
        nd.ctr[a] = ctr[a];
        nd.half[a] = half[a];
    }
    return g_nodeCount++;
}

void Simulation::insertBody(int ni, int b, int depth) const {
    uint8_t st = g_nodes[ni].state;
    if (st == 0) {  // 空单元
        g_nodes[ni].state = 1;
        g_nodes[ni].body = b;
        return;
    }
    if (st == 3) {  // 深度封顶的 bucket
        g_nodes[ni].extra = extraAlloc(b, g_nodes[ni].extra);
        return;
    }
    if (st == 1) {  // 被占用的叶 → 细分
        if (depth >= kMaxDepth) {
            g_nodes[ni].state = 3;
            g_nodes[ni].extra = extraAlloc(b, g_nodes[ni].extra);
            return;
        }
        int old = g_nodes[ni].body;
        g_nodes[ni].body = -1;
        g_nodes[ni].state = 2;
        insertBody(ni, old, depth + 1);
        insertBody(ni, b, depth + 1);
        return;
    }
    // 内部节点：下沉到包含该天体的子单元
    int k = 0;
    for (int a = 0; a < dim_; a++)
        if (bodies_[b].pos[a] >= g_nodes[ni].ctr[a]) k |= (1 << a);
    int c = g_nodes[ni].child[k];
    if (c < 0) {
        double ctr[8], half[8];
        for (int a = 0; a < dim_; a++) {
            double h2 = g_nodes[ni].half[a] * 0.5;
            ctr[a] = g_nodes[ni].ctr[a] + (((k >> a) & 1) ? h2 : -h2);
            half[a] = h2;
        }
        c = newNode(ni, ctr, half);
        if (c < 0) {  // 节点耗尽：退化成 bucket（力仍精确，只是失去剪枝）
            g_nodes[ni].state = 3;
            g_nodes[ni].extra = extraAlloc(b, g_nodes[ni].extra);
            return;
        }
        g_nodes[ni].child[k] = c;
    }
    insertBody(c, b, depth + 1);
}

void Simulation::buildTree() const {
    g_nodeCount = 0;
    g_root = -1;
    g_extraUsed = 0;
    double mn[8], mx[8];
    for (int k = 0; k < dim_; k++) {
        mn[k] = __builtin_inf();
        mx[k] = -__builtin_inf();
    }
    bool any = false;
    for (int i = 0; i < nbodies_; i++) {
        const Body& b = bodies_[i];
        if (!b.alive) continue;
        any = true;
        for (int k = 0; k < dim_; k++) {
            mn[k] = gs_min2(mn[k], b.pos[k]);
            mx[k] = gs_max2(mx[k], b.pos[k]);
        }
    }
    if (!any) return;
    double ctr[8], half[8];
    for (int k = 0; k < dim_; k++) {
        ctr[k] = 0.5 * (mn[k] + mx[k]);
        half[k] = 0.5 * (mx[k] - mn[k]) * 1.000001 + 1e-12;
        if (!(half[k] > 0)) half[k] = 1e-6;
    }
    g_root = newNode(-1, ctr, half);
    if (g_root < 0) return;
    for (int i = 0; i < nbodies_; i++)
        if (bodies_[i].alive) insertBody(g_root, i, 0);

    // 自底向上：质量 / 质心 / 最大半径。
    // 子节点索引恒大于父节点，所以倒序遍历保证「子先于父结算」。
    // ⚠ 不变量：扫过程中 nd.com 保存的是**原始加权和** Σ(pos·m)，不是质心；
    //   父子之间原样相加，最后一趟才统一除以质量。
    for (int n = g_nodeCount - 1; n >= 0; n--) {
        TreeNode& nd = g_nodes[n];
        if (nd.state == 1) {
            const Body& b = bodies_[nd.body];
            nd.mass = b.mass;
            for (int k = 0; k < dim_; k++) nd.com[k] = b.pos[k] * b.mass;
            nd.maxRad = b.radius;
        } else if (nd.state == 3) {
            double m = 0.0;
            double c[8];
            for (int k = 0; k < 8; k++) c[k] = 0.0;
            for (int e = nd.extra; e >= 0; e = g_extraNext[e]) {
                const Body& b = bodies_[g_extraBody[e]];
                m += b.mass;
                for (int k = 0; k < dim_; k++) c[k] += b.pos[k] * b.mass;
                nd.maxRad = gs_max2(nd.maxRad, b.radius);
            }
            if (nd.body >= 0) {
                const Body& b = bodies_[nd.body];
                m += b.mass;
                for (int k = 0; k < dim_; k++) c[k] += b.pos[k] * b.mass;
                nd.maxRad = gs_max2(nd.maxRad, b.radius);
            }
            nd.mass = m;
            for (int k = 0; k < dim_; k++) nd.com[k] = c[k];
        }
        if (n != g_root && nd.mass > 0) {
            TreeNode& p = g_nodes[nd.parent];
            for (int k = 0; k < dim_; k++) p.com[k] += nd.com[k];
            p.mass += nd.mass;
            p.maxRad = gs_max2(p.maxRad, nd.maxRad);
        } else if (n != g_root) {
            TreeNode& p = g_nodes[nd.parent];
            p.maxRad = gs_max2(p.maxRad, nd.maxRad);
        }
    }
    // 收尾：加权和除以质量得到质心；计算度量半径界 mrad
    const int nc = 1 << dim_;
    for (int i = 0; i < g_nodeCount; i++) {
        TreeNode& nd = g_nodes[i];
        if (nd.mass > 0) {
            for (int k = 0; k < dim_; k++) nd.com[k] /= nd.mass;
        } else {
            for (int k = 0; k < dim_; k++) nd.com[k] = nd.ctr[k];
        }
        double mr = 0.0;
        for (int c = 0; c < nc; c++) {
            double dd[8];
            for (int k = 0; k < dim_; k++)
                dd[k] = (nd.ctr[k] + (((c >> k) & 1) ? nd.half[k] : -nd.half[k])) -
                        nd.com[k];
            mr = gs_max2(mr, metricDistance(dd));
        }
        nd.mrad = mr;
    }
}

// 点 xi 到盒内任意点的度量距离**精确下界**。
// 对可分离族 d = F(Σ f(|dx_k|))，目标 Σ f(|dx_k|) 各轴独立最小，即逐轴 clamp。
double Simulation::minMetricDistToBox(const double* xi, const TreeNode& nd) const {
    double dd[8];
    for (int k = 0; k < dim_; k++) {
        double lo = nd.ctr[k] - nd.half[k], hi = nd.ctr[k] + nd.half[k];
        double t = xi[k] < lo ? lo : (xi[k] > hi ? hi : xi[k]);
        dd[k] = t - xi[k];
    }
    return metricDistance(dd);
}

// 中心力：F > 0 把 i 拉向 posj 处的源。
//   a_i += m_j · s · (posj − x_i)/|posj − x_i|
//   s = F(d)（核外）或 F(R)·d/R（核内线性正则）
void Simulation::applyOneSided(int i, const double* posj, double mj, double radj) {
    double delta[8];
    double rho2 = 0.0;
    for (int k = 0; k < dim_; k++) {
        delta[k] = posj[k] - bodies_[i].pos[k];
        rho2 += delta[k] * delta[k];
    }
    if (!(rho2 > 0.0)) return;      // 重合：没有方向，力为 0
    if (bodies_[i].fixed) return;   // 钉住天体施力但不受力
    double d = metricDistance(delta);
    if (!(d > 0.0) || !gs_isfinite(d)) return;
    const double Rij = bodies_[i].radius + radj;
    double s;
    if (Rij > 0.0 && d < Rij) {
        s = G_ * F_.value(Rij) * (d / Rij);
    } else {
        s = G_ * F_.value(d);
    }
    if (!gs_isfinite(s)) return;
    if (gs_fabs(s) > kForceClamp) s = s > 0 ? kForceClamp : -kForceClamp;
    if (s == 0.0) return;
    const double w = s / gs_sqrt(rho2);  // 每单位欧氏长度的力
    for (int k = 0; k < dim_; k++) bodies_[i].acc[k] += mj * w * delta[k];
}

double Simulation::pairUOneSided(int i, const double* posj, double mj, double radj) const {
    double delta[8];
    for (int k = 0; k < dim_; k++) delta[k] = posj[k] - bodies_[i].pos[k];
    double d = metricDistance(delta);
    if (!gs_isfinite(d)) return 0.0;
    const double Rij = bodies_[i].radius + radj;
    double val = pairV(d, Rij);
    if (!gs_isfinite(val)) return 0.0;
    return bodies_[i].mass * mj * val;
}

void Simulation::treeForce() {
    for (int i = 0; i < nbodies_; i++)
        for (int k = 0; k < dim_; k++) bodies_[i].acc[k] = 0.0;
    if (g_root < 0) return;
    int stack[512];
    const int nc = 1 << dim_;
    for (int i = 0; i < nbodies_; i++) {
        if (!bodies_[i].alive) continue;
        if (bodies_[i].fixed) continue;  // 钉住天体没有属于自己的力
        double xi[8];
        for (int k = 0; k < dim_; k++) xi[k] = bodies_[i].pos[k];
        int sp = 0;
        stack[sp++] = g_root;
        while (sp > 0) {
            int ni = stack[--sp];
            const TreeNode& nd = g_nodes[ni];
            if (nd.mass <= 0.0) continue;
            if (nd.state == 1) {
                if (nd.body != i)
                    applyOneSided(i, bodies_[nd.body].pos, bodies_[nd.body].mass,
                                  bodies_[nd.body].radius);
            } else if (nd.state == 3) {
                if (nd.body >= 0 && nd.body != i)
                    applyOneSided(i, bodies_[nd.body].pos, bodies_[nd.body].mass,
                                  bodies_[nd.body].radius);
                for (int e = nd.extra; e >= 0; e = g_extraNext[e]) {
                    int j = g_extraBody[e];
                    if (j != i)
                        applyOneSided(i, bodies_[j].pos, bodies_[j].mass, bodies_[j].radius);
                }
            } else {
                double delta[8];
                for (int k = 0; k < dim_; k++) delta[k] = nd.com[k] - xi[k];
                double d = metricDistance(delta);
                if (d > 0.0 && nd.mrad < theta_ * d) {
                    applyOneSided(i, nd.com, nd.mass, 0.0);
                } else {
                    for (int k = 0; k < nc; k++) {
                        int c = nd.child[k];
                        if (c >= 0 && g_nodes[c].mass > 0 && sp < 511) stack[sp++] = c;
                    }
                }
            }
        }
    }
}

double Simulation::treePotential() const {
    if (g_root < 0) return 0.0;
    // 势能容忍比力更紧的开角判据：θ_PE = 0.25 时估计误差 ~0.3%
    // （θ = 0.6 时约 −7%）；PE 每帧只算一次，力每子步都算。
    const double thPE = gs_min2(theta_, 0.25);
    double U = 0.0;
    int stack[512];
    const int nc = 1 << dim_;
    for (int i = 0; i < nbodies_; i++) {
        if (!bodies_[i].alive) continue;
        double xi[8];
        for (int k = 0; k < dim_; k++) xi[k] = bodies_[i].pos[k];
        double Ui = 0.0;
        int sp = 0;
        stack[sp++] = g_root;
        while (sp > 0) {
            int ni = stack[--sp];
            const TreeNode& nd = g_nodes[ni];
            if (nd.mass <= 0.0) continue;
            if (nd.state == 1) {
                if (nd.body != i)
                    Ui += pairUOneSided(i, bodies_[nd.body].pos, bodies_[nd.body].mass,
                                        bodies_[nd.body].radius);
            } else if (nd.state == 3) {
                if (nd.body >= 0 && nd.body != i)
                    Ui += pairUOneSided(i, bodies_[nd.body].pos, bodies_[nd.body].mass,
                                        bodies_[nd.body].radius);
                for (int e = nd.extra; e >= 0; e = g_extraNext[e]) {
                    int j = g_extraBody[e];
                    if (j != i)
                        Ui += pairUOneSided(i, bodies_[j].pos, bodies_[j].mass,
                                            bodies_[j].radius);
                }
            } else {
                double delta[8];
                for (int k = 0; k < dim_; k++) delta[k] = nd.com[k] - xi[k];
                double d = metricDistance(delta);
                if (d > 0.0 && nd.mrad < thPE * d) {
                    Ui += pairUOneSided(i, nd.com, nd.mass, 0.0);
                } else {
                    for (int k = 0; k < nc; k++) {
                        int c = nd.child[k];
                        if (c >= 0 && g_nodes[c].mass > 0 && sp < 511) stack[sp++] = c;
                    }
                }
            }
        }
        U += Ui;
    }
    return 0.5 * U;
}

void Simulation::collectContacts() const {
    g_cN = 0;
    if (g_root < 0) return;
    int stack[512];
    const int nc = 1 << dim_;
    for (int i = 0; i < nbodies_; i++) {
        if (!bodies_[i].alive) continue;
        const double Ri = bodies_[i].radius;
        double xi[8];
        for (int k = 0; k < dim_; k++) xi[k] = bodies_[i].pos[k];
        int sp = 0;
        stack[sp++] = g_root;
        while (sp > 0) {
            int ni = stack[--sp];
            const TreeNode& nd = g_nodes[ni];
            if (nd.mass <= 0.0) continue;
            if (nd.state == 1 || nd.state == 3) {
                // 叶：只测 nd.body；bucket：先测 extra 链表，再测 nd.body
                if (nd.state == 3) {
                    for (int e = nd.extra; e >= 0 && g_cN < kMaxContacts;
                         e = g_extraNext[e]) {
                        int j = g_extraBody[e];
                        if (j <= i || !bodies_[j].alive) continue;
                        double delta[8];
                        for (int k = 0; k < dim_; k++)
                            delta[k] = bodies_[j].pos[k] - xi[k];
                        double d = metricDistance(delta);
                        if (gs_isfinite(d) && d <= Ri + bodies_[j].radius) {
                            g_cA[g_cN] = i;
                            g_cB[g_cN] = j;
                            g_cN++;
                        }
                    }
                }
                if (nd.body >= 0 && g_cN < kMaxContacts) {
                    int j = nd.body;
                    if (j > i && bodies_[j].alive) {
                        double delta[8];
                        for (int k = 0; k < dim_; k++)
                            delta[k] = bodies_[j].pos[k] - xi[k];
                        double d = metricDistance(delta);
                        if (gs_isfinite(d) && d <= Ri + bodies_[j].radius) {
                            g_cA[g_cN] = i;
                            g_cB[g_cN] = j;
                            g_cN++;
                        }
                    }
                }
            } else {
                if (minMetricDistToBox(xi, nd) > Ri + nd.maxRad) continue;  // 精确剪枝
                for (int k = 0; k < nc; k++) {
                    int c = nd.child[k];
                    if (c >= 0 && g_nodes[c].mass > 0 && sp < 511) stack[sp++] = c;
                }
            }
        }
    }
}

// -------------------------------------------------------------------- 天体

int Simulation::addBody(double mass, const double* pos, const double* vel, bool fixed) {
    if (!ok_ || mass <= 0 || !gs_isfinite(mass)) return -1;
    if (nbodies_ >= kMaxBodies) return -1;
    Body& b = bodies_[nbodies_];
    b.mass = mass;
    b.radius = radiusFromMass(mass);
    b.alive = true;
    b.fixed = fixed;
    for (int k = 0; k < dim_; k++) {
        b.pos[k] = pos[k];
        b.vel[k] = fixed ? 0.0 : vel[k];
        b.acc[k] = 0.0;
        if (!gs_isfinite(b.pos[k])) return -1;
        if (!gs_isfinite(b.vel[k])) return -1;
    }
    nbodies_++;
    accValid_ = false;
    return nbodies_ - 1;
}

void Simulation::clearBodies() {
    nbodies_ = 0;
    t_ = 0.0;
    accValid_ = false;
}

// -------------------------------------------------------------- 度量核

// 解 f(x) = S（x >= 0），f 严格递增
double Simulation::invertF(double S) const {
    if (!(S > 0.0)) return 0.0;
    if (fIsPower_) return gs_pow(S, 1.0 / fExp_);

    double lo = 0.0, hi = 1.0;
    double fhi = f_.value(hi);
    int expand = 0;
    while (fhi < S) {
        lo = hi;
        hi *= 2.0;
        fhi = f_.value(hi);
        if (++expand > 200) return hi;  // f 有界：距离视为无界
    }
    if (f_.value(lo) > S) return 0.0;

    // 带保护的 Newton（导数来自对偶数求值）
    double x = 0.5 * (lo + hi);
    for (int it = 0; it < 80; it++) {
        Dual fx = f_.eval(x);
        if (!gs_isfinite(fx.v)) break;
        if (gs_fabs(fx.v - S) <= 1e-14 * gs_max2(1.0, S)) return x;
        if (fx.v < S) lo = x;
        else hi = x;
        if (hi - lo <= 1e-15 * gs_max2(1.0, hi)) return 0.5 * (lo + hi);
        double xn;
        if (fx.d > 1e-300) xn = x - (fx.v - S) / fx.d;
        else xn = 0.5 * (lo + hi);
        if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
        x = xn;
    }
    return 0.5 * (lo + hi);
}

double Simulation::metricDistance(const double* delta) const {
    if (fIsPower_) {
        const double p = fExp_;
        double S = 0.0;
        for (int k = 0; k < dim_; k++) S += gs_pow(gs_fabs(delta[k]), p);
        return S > 0.0 ? gs_pow(S, 1.0 / p) : 0.0;
    }
    double S = 0.0;
    for (int k = 0; k < dim_; k++) S += f_.value(gs_fabs(delta[k]));
    return invertF(S);
}

bool Simulation::pairKernel(const double* delta, double& dOut, double* grad) const {
    double d;
    if (fIsPower_) {
        const double p = fExp_;
        double S = 0.0;
        for (int k = 0; k < dim_; k++) S += gs_pow(gs_fabs(delta[k]), p);
        d = S > 0.0 ? gs_pow(S, 1.0 / p) : 0.0;
        if (!(d > 0.0) || !gs_isfinite(d)) {
            dOut = d > 0 ? d : 0.0;
            for (int k = 0; k < dim_; k++) grad[k] = 0.0;
            return d > 0.0 && gs_isfinite(d);
        }
        // grad[k] = (|dx_k|/d)^(p−1) · sgn(dx_k)
        for (int k = 0; k < dim_; k++) {
            double ratio = gs_fabs(delta[k]) / d;  // ∈ [0,1]
            double g = gs_pow(ratio, p - 1.0);
            if (!gs_isfinite(g))
                g = (p < 1.0 && ratio > 0.0) ? kGradClamp : (ratio == 0.0 ? 0.0 : 1.0);
            g *= (delta[k] > 0.0) - (delta[k] < 0.0);
            grad[k] = gs_clamp(g, -kGradClamp, kGradClamp);
        }
    } else {
        double S = 0.0;
        double fp[8];
        for (int k = 0; k < dim_; k++) {
            Dual fk = f_.eval(gs_fabs(delta[k]));
            S += fk.v;
            fp[k] = fk.d;
        }
        d = invertF(S);
        if (!gs_isfinite(d)) {
            dOut = 0.0;
            for (int k = 0; k < dim_; k++) grad[k] = 0.0;
            return false;
        }
        double fpd = f_.eval(d).d;
        if (gs_fabs(fpd) < 1e-300) fpd = fpd >= 0 ? 1e-300 : -1e-300;
        for (int k = 0; k < dim_; k++) {
            double g = fp[k] / fpd;
            g *= (delta[k] > 0.0) - (delta[k] < 0.0);
            grad[k] = gs_clamp(g, -kGradClamp, kGradClamp);
        }
    }
    dOut = d;
    return true;
}

// ------------------------------------------------------------------ 物理

void Simulation::computeAccel() {
    if (treeActive()) {
        buildTree();
        treeForce();
        accValid_ = true;
        return;
    }
    for (int i = 0; i < nbodies_; i++)
        for (int k = 0; k < dim_; k++) bodies_[i].acc[k] = 0.0;
    double delta[8];
    for (int i = 0; i < nbodies_; i++) {
        if (!bodies_[i].alive) continue;
        for (int j = i + 1; j < nbodies_; j++) {
            if (!bodies_[j].alive) continue;
            double rho2 = 0.0;
            for (int k = 0; k < dim_; k++) {
                delta[k] = bodies_[j].pos[k] - bodies_[i].pos[k];
                rho2 += delta[k] * delta[k];
            }
            if (!(rho2 > 0.0)) continue;  // 重合
            double d = metricDistance(delta);
            if (!(d > 0.0) || !gs_isfinite(d)) continue;
            const double Rij = bodies_[i].radius + bodies_[j].radius;
            // 核外 G·F(d)；核内线性正则 G·F(R)·(d/R)（在 d = R 处力连续）
            double s;
            if (Rij > 0.0 && d < Rij) {
                s = G_ * F_.value(Rij) * (d / Rij);
            } else {
                s = G_ * F_.value(d);
            }
            if (!gs_isfinite(s)) continue;
            if (gs_fabs(s) > kForceClamp) s = s > 0 ? kForceClamp : -kForceClamp;
            if (s == 0.0) continue;
            const double mi = bodies_[i].mass, mj = bodies_[j].mass;
            const double w = s / gs_sqrt(rho2);
            // 中心力：等大反向，沿连线（钉住天体施力但不受力）
            for (int k = 0; k < dim_; k++) {
                double f = w * delta[k];
                if (!bodies_[i].fixed) bodies_[i].acc[k] += mj * f;
                if (!bodies_[j].fixed) bodies_[j].acc[k] -= mi * f;
            }
        }
    }
    accValid_ = true;
}

void Simulation::resolveContact(int i, int j) {
    Body& A = bodies_[i];
    Body& B = bodies_[j];
    if (!A.alive || !B.alive) return;
    double delta[8];
    for (int k = 0; k < dim_; k++) delta[k] = B.pos[k] - A.pos[k];
    double d = metricDistance(delta);
    if (!gs_isfinite(d)) return;
    const double Rij = A.radius + B.radius;
    if (d > Rij) return;

    if (collisionMode_ == Collision::Merge) {
        // 钉住天体永不参与并合（晕必须保持完整）：退化成穿过（核内正则仍生效）
        if (A.fixed || B.fixed) return;
        // 把 j 并进 i（完全非弹性，动量守恒）
        double M = A.mass + B.mass;
        for (int k = 0; k < dim_; k++) {
            A.pos[k] = (A.pos[k] * A.mass + B.pos[k] * B.mass) / M;
            A.vel[k] = (A.vel[k] * A.mass + B.vel[k] * B.mass) / M;
            A.acc[k] = 0.0;
        }
        A.mass = M;
        A.radius = radiusFromMass(M);
        B.alive = false;
        computeAccel();
    } else {  // Collision::Bounce
        // 沿接触法向的完全弹性脉冲。与中心力模型一致，法向即两心连线
        // （f = t^2 时与度量梯度法向完全相同）。
        double gn2 = 0.0;
        for (int k = 0; k < dim_; k++) gn2 += delta[k] * delta[k];
        double gn = gs_sqrt(gn2);
        if (gn < 1e-300) return;
        double nrm[8];
        double w = 0.0;  // 沿法向的接近速度
        for (int k = 0; k < dim_; k++) {
            nrm[k] = delta[k] / gn;  // A → B
            w += (A.vel[k] - B.vel[k]) * nrm[k];
        }
        if (w <= 0.0) return;  // 已经在分离
        // 逆质量形式的冲量：钉住天体逆质量为 0 ⇒ 无穷大惯性，只反射对方
        const double invA = A.fixed ? 0.0 : 1.0 / A.mass;
        const double invB = B.fixed ? 0.0 : 1.0 / B.mass;
        const double tot = invA + invB;
        if (tot <= 0.0) return;  // 两者都被钉住
        const double e = 1.0;    // 完全弹性
        for (int k = 0; k < dim_; k++) {
            A.vel[k] -= (1.0 + e) * (invA / tot) * w * nrm[k];
            B.vel[k] += (1.0 + e) * (invB / tot) * w * nrm[k];
        }
    }
}

void Simulation::handleCollisions() {
    if (collisionMode_ == Collision::Pass) return;
    if (treeActive()) {
        // 本子步 computeAccel() 建的树仍然有效（kick 不移动位置）；
        // resolveContact() 内部的并合会自己重算加速度。
        collectContacts();
        for (int c = 0; c < g_cN; c++) resolveContact(g_cA[c], g_cB[c]);
        return;
    }
    for (int i = 0; i < nbodies_; i++) {
        if (!bodies_[i].alive) continue;
        for (int j = i + 1; j < nbodies_; j++) {
            if (!bodies_[j].alive) continue;
            resolveContact(i, j);
        }
    }
}

bool Simulation::step(double dt, int substeps) {
    if (!ok_) return false;
    if (substeps < 1) substeps = 1;
    if (!(dt > 0.0) || !gs_isfinite(dt)) return false;
    if (!accValid_) computeAccel();
    const double h = dt;
    for (int s = 0; s < substeps; s++) {
        // 半 kick + 漂移
        for (int i = 0; i < nbodies_; i++) {
            Body& b = bodies_[i];
            if (!b.alive || b.fixed) continue;  // 钉住天体永不移动
            for (int k = 0; k < dim_; k++) {
                b.vel[k] += 0.5 * h * b.acc[k];
                b.pos[k] += b.vel[k] * h;
            }
        }
        computeAccel();  // 新加速度
        for (int i = 0; i < nbodies_; i++) {
            Body& b = bodies_[i];
            if (!b.alive || b.fixed) continue;
            for (int k = 0; k < dim_; k++) b.vel[k] += 0.5 * h * b.acc[k];
        }
        handleCollisions();
        t_ += h;
    }
    return true;
}

int Simulation::aliveCount() const {
    int c = 0;
    for (int i = 0; i < nbodies_; i++)
        if (bodies_[i].alive) c++;
    return c;
}

double Simulation::potentialEnergy() const {
    if (treeActive()) {
        buildTree();
        return treePotential();
    }
    double U = 0.0;
    double delta[8];
    for (int i = 0; i < nbodies_; i++) {
        if (!bodies_[i].alive) continue;
        for (int j = i + 1; j < nbodies_; j++) {
            if (!bodies_[j].alive) continue;
            for (int k = 0; k < dim_; k++) delta[k] = bodies_[j].pos[k] - bodies_[i].pos[k];
            double d = metricDistance(delta);
            if (!gs_isfinite(d)) continue;
            const double Rij = bodies_[i].radius + bodies_[j].radius;
            double val = pairV(d, Rij);
            if (gs_isfinite(val)) U += bodies_[i].mass * bodies_[j].mass * val;
        }
    }
    return U;
}

double Simulation::energy() const {
    double KE = 0.0;
    for (int i = 0; i < nbodies_; i++) {
        const Body& b = bodies_[i];
        if (!b.alive) continue;
        double v2 = 0.0;
        for (int k = 0; k < dim_; k++) v2 += b.vel[k] * b.vel[k];
        KE += 0.5 * b.mass * v2;
    }
    return KE + potentialEnergy();
}

int Simulation::getState(double* out, int maxBodies) const {
    int w = 0;
    for (int i = 0; i < nbodies_; i++) {
        const Body& b = bodies_[i];
        if (!b.alive) continue;
        if (w >= maxBodies) break;
        out[0] = b.mass;
        for (int k = 0; k < 3; k++) out[1 + k] = k < dim_ ? b.pos[k] : 0.0;
        for (int k = 0; k < 3; k++) out[4 + k] = k < dim_ ? b.vel[k] : 0.0;
        out[7] = b.fixed ? (double)kFlagFixed : 0.0;
        out += 8;
        w++;
    }
    return w;
}

void Simulation::samplePotential(int n, double rmin, double rmax, double* outR,
                                 double* outV) const {
    if (n < 1) return;
    if (!(rmin > 0.0)) rmin = 1e-6;
    if (rmax <= rmin) rmax = rmin * 1e6;
    for (int i = 0; i < n; i++) {
        double frac = (n == 1) ? 0.0 : (double)i / (n - 1);
        double r = gs_exp(gs_log(rmin) + frac * (gs_log(rmax) - gs_log(rmin)));
        double v = potentialAt(r);  // 派生势 V = ∫F
        outR[i] = r;
        outV[i] = gs_isfinite(v) ? v : __builtin_nan("");
    }
}

int Simulation::unitBall2D(int n, double* outXY) const {
    if (n < 3) return 0;
    for (int i = 0; i < n; i++) {
        double th = 2.0 * GS_PI * i / n;
        double ux = gs_cos(th), uy = gs_sin(th);
        // 二分求解 metricDistance(rho·u) == 1
        double lo = 0.0, hi = 1.0;
        int expand = 0;
        double d[8];
        for (int k = 0; k < 8; k++) d[k] = 0.0;
        d[0] = ux;
        d[1] = uy;
        // 度量对 rho 单调：用缩放后的向量求距离
        auto dist = [&](double rho) {
            double v[8];
            for (int k = 0; k < 8; k++) v[k] = 0.0;
            v[0] = rho * ux;
            v[1] = rho * uy;
            return metricDistance(v);
        };
        while (dist(hi) < 1.0 && expand < 60) {
            hi *= 2.0;
            expand++;
        }
        if (dist(hi) < 1.0) {  // 不可达（f 有界）
            outXY[2 * i] = __builtin_nan("");
            outXY[2 * i + 1] = __builtin_nan("");
            continue;
        }
        for (int it = 0; it < 60; it++) {
            double mid = 0.5 * (lo + hi);
            if (dist(mid) < 1.0) lo = mid;
            else hi = mid;
        }
        double rho = 0.5 * (lo + hi);
        outXY[2 * i] = rho * ux;
        outXY[2 * i + 1] = rho * uy;
    }
    return n;
}

}  // namespace grav
