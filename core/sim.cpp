#include "sim.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>

namespace grav {

static const double kGradClamp = 1e6;   // clamp on |d(d)/d(dx_k)|
static const double kForceClamp = 1e12; // clamp on force scale
static const double kPi = 3.14159265358979323846;

Simulation::Simulation(int dims, const std::string& fExpr, const std::string& forceExpr)
    : dim_(dims < 1 ? 1 : dims), f_(fExpr, "t"), F_(forceExpr, "r") {
    if (dim_ > 8) {
        err_ = "内核最多支持 8 维（前端只开放 2D/3D）";
        return;
    }
    if (!f_.ok()) {
        err_ = "f(t) 表达式错误: " + f_.error();
        return;
    }
    if (!F_.ok()) {
        err_ = "F(r) 表达式错误: " + F_.error();
        return;
    }
    // Force law sanity: must be finite at a few probe distances.  (A pole
    // elsewhere, e.g. F = 1/(r-1), is handled by the runtime clamps.)
    for (double r : {0.5, 1.0, 2.0}) {
        if (!std::isfinite(F_.value(r))) {
            err_ = "F(r) 在 r = " + std::to_string(r) + " 处不是有限值";
            return;
        }
    }

    double p = 0.0;
    if (f_.asPowerOfVar(p)) {
        if (!(p > 0.0) || !std::isfinite(p)) {
            err_ = "f(t) = t^p 要求 p > 0 才能求逆";
            return;
        }
        fIsPower_ = true;
        fExp_ = p;
    } else {
        // Validate: f(0) == 0 and f strictly increasing on a sampled grid.
        if (std::fabs(f_.value(0.0)) > 1e-9) {
            err_ = "f 必须满足 f(0) = 0（否则两点重合时距离不为 0）";
            return;
        }
        double prev = 0.0;
        for (int i = 0; i < 80; i++) {
            double t = std::pow(10.0, -6.0 + 10.0 * i / 79.0);  // 1e-6 .. 1e4
            double v = f_.value(t);
            if (!std::isfinite(v)) {
                err_ = "f 在 t = " + std::to_string(t) + " 处不是有限值";
                return;
            }
            if (!(v > prev)) {  // strict increase
                err_ = "f 必须在 [0, +inf) 上严格递增（在 t = " + std::to_string(t) + " 附近不满足）";
                return;
            }
            prev = v;
        }
    }

    // Derive the potential V(r) = integral of F dr  (so that V' = F).
    std::vector<std::pair<double, double>> terms;  // F = sum a_i r^{p_i}
    if (F_.asPowerSum(terms)) {
        vAnalytic_ = true;  // exact antiderivative
        for (const auto& tm : terms) {
            double pp = tm.first, a = tm.second;
            if (pp == -1.0) {
                vLn_ += a;                      // int r^-1 = ln r
            } else {
                vPow_.push_back({pp + 1.0, a / (pp + 1.0)});
            }
        }
    } else {
        buildPotentialTable();                  // numeric quadrature
    }
    ok_ = true;
}

// ------------------------------------------------- derived potential V

// Numeric fallback: V sampled on a log grid over [1e-6, 1e6], integrated
// with 4-point Gauss-Legendre per segment in log space (dr = r dln r).
// Interpolation is cubic Hermite in ln r with the exact slope dV/dln r =
// F(r)*r at the nodes; outside the grid V is extended linearly with the
// end slope F.
void Simulation::buildPotentialTable() {
    const int N = 4096;
    const double lmin = std::log(1e-6), lmax = std::log(1e6);
    static const double xg[4] = {-0.8611363115940526, -0.3399810435848563,
                                 0.3399810435848563, 0.8611363115940526};
    static const double wg[4] = {0.34785484513745385, 0.6521451548625461,
                                 0.6521451548625461, 0.34785484513745385};
    vTabR_.assign(N, 0.0);
    vTabV_.assign(N, 0.0);
    vTabF_.assign(N, 0.0);
    const double dl = (lmax - lmin) / (N - 1);
    double V = 0.0;
    for (int i = 0; i < N; i++) {
        double l = lmin + dl * i;
        double r = std::exp(l);
        double fr = F_.value(r);
        vTabR_[i] = r;
        vTabV_[i] = V;
        vTabF_[i] = std::isfinite(fr) ? fr : 0.0;
        if (i + 1 < N) {
            double lm = l + 0.5 * dl;
            double seg = 0.0;
            for (int g = 0; g < 4; g++) {
                double lg = lm + 0.5 * dl * xg[g];
                double rg = std::exp(lg);
                double fg = F_.value(rg);
                if (std::isfinite(fg)) seg += wg[g] * 0.5 * dl * fg * rg;
            }
            V += seg;
            if (!std::isfinite(V)) V = V > 0 ? 1e300 : -1e300;
            if (std::fabs(V) > 1e300) V = V > 0 ? 1e300 : -1e300;
        }
    }
}

double Simulation::potentialAt(double r) const {
    if (!(r > 0.0) || !std::isfinite(r)) return 0.0;
    double v;
    if (vAnalytic_) {
        v = vLn_ * std::log(r);
        for (const auto& tm : vPow_) v += tm.second * std::pow(r, tm.first);
    } else {
        const int N = (int)vTabR_.size();
        if (N < 2) return 0.0;
        const double lmin = std::log(vTabR_[0]);
        const double dl = std::log(vTabR_[1]) - lmin;
        const double lr = std::log(r);
        if (lr <= lmin) {  // below the grid: linear with the first slope
            v = vTabV_[0] + vTabF_[0] * (r - vTabR_[0]);
        } else {
            const double lmax = lmin + dl * (N - 1);
            if (lr >= lmax) {  // above the grid: linear with the last slope
                v = vTabV_[N - 1] + vTabF_[N - 1] * (r - vTabR_[N - 1]);
            } else {
                double fi = (lr - lmin) / dl;
                int i = (int)fi;
                if (i > N - 2) i = N - 2;
                double u = fi - i;
                // cubic Hermite with slopes dV/dln r = F(r) * r
                double s0 = vTabF_[i] * vTabR_[i];
                double s1 = vTabF_[i + 1] * vTabR_[i + 1];
                double h00 = (1 + 2 * u) * (1 - u) * (1 - u);
                double h10 = u * (1 - u) * (1 - u);
                double h01 = u * u * (3 - 2 * u);
                double h11 = u * u * (u - 1);
                v = h00 * vTabV_[i] + h10 * dl * s0 + h01 * vTabV_[i + 1] + h11 * dl * s1;
            }
        }
    }
    v = G_ * v;  // coupling constant scales the derived potential
    return std::isfinite(v) ? v : 0.0;
}

double Simulation::forceAt(double r) const { return F_.value(r); }

// Pair potential at metric distance d with core radius R = R_i + R_j.
// Inside the core the quadratic matches W(R) = V(R), W'(R) = G*F(R).
double Simulation::pairV(double d, double R) const {
    if (R > 0.0 && d < R) {
        double FR = G_ * F_.value(R);
        if (!std::isfinite(FR)) FR = 0.0;
        return potentialAt(R) + FR * (d * d - R * R) / (2.0 * R);
    }
    return potentialAt(d);
}

void Simulation::setCollision(Collision mode, double density) {
    collisionMode_ = mode;
    density_ = density > 1e-12 ? density : 1e-12;
    // density changes body radii
    for (auto& b : bodies_) b.radius = radiusFromMass(b.mass);
}

void Simulation::setG(double G) {
    if (G > 0.0 && std::isfinite(G)) G_ = G;
}

double Simulation::radiusFromMass(double m) const {
    // R = ( m * Gamma(d/2+1) / (rho * pi^(d/2)) )^(1/d)
    if (!(m > 0.0)) return 0.0;
    double g = std::tgamma(dim_ / 2.0 + 1.0);
    return std::pow(m * g / (density_ * std::pow(kPi, dim_ / 2.0)), 1.0 / dim_);
}

// ---------------------------------------------------------- Barnes-Hut tree

void Simulation::setTree(int treeMode, double theta) {
    treeMode_ = (treeMode == 1) ? 1 : 0;
    theta_ = theta > 0.05 ? theta : 0.05;
}

bool Simulation::treeActive() const {
    if (treeMode_ != 0) return false;
    int alive = 0;
    for (const auto& b : bodies_)
        if (b.alive) alive++;
    return alive >= kTreeThreshold;
}

int Simulation::newNode(int parent, const double* ctr, const double* half) const {
    TreeNode nd;
    nd.parent = parent;
    for (int a = 0; a < dim_; a++) {
        nd.ctr[a] = ctr[a];
        nd.half[a] = half[a];
    }
    nodes_.push_back(std::move(nd));
    return (int)nodes_.size() - 1;
}

void Simulation::insertBody(int ni, int b, int depth) const {
    uint8_t st = nodes_[ni].state;
    if (st == 0) {  // empty cell
        nodes_[ni].state = 1;
        nodes_[ni].body = b;
        return;
    }
    if (st == 3) {  // depth-capped bucket
        nodes_[ni].extra.push_back(b);
        return;
    }
    if (st == 1) {  // occupied leaf -> subdivide
        if (depth >= kMaxDepth) {
            nodes_[ni].state = 3;
            nodes_[ni].extra.push_back(b);
            return;
        }
        int old = nodes_[ni].body;
        nodes_[ni].body = -1;
        nodes_[ni].state = 2;
        insertBody(ni, old, depth + 1);
        insertBody(ni, b, depth + 1);
        return;
    }
    // internal: descend into the child cell containing the body
    int k = 0;
    for (int a = 0; a < dim_; a++)
        if (bodies_[b].pos[a] >= nodes_[ni].ctr[a]) k |= (1 << a);
    int c = nodes_[ni].child[k];
    if (c < 0) {
        double ctr[8], half[8];
        for (int a = 0; a < dim_; a++) {
            double h2 = nodes_[ni].half[a] * 0.5;
            ctr[a] = nodes_[ni].ctr[a] + (((k >> a) & 1) ? h2 : -h2);
            half[a] = h2;
        }
        c = newNode(ni, ctr, half);
        nodes_[ni].child[k] = c;
    }
    insertBody(c, b, depth + 1);
}

void Simulation::buildTree() const {
    nodes_.clear();
    root_ = -1;
    double mn[8], mx[8];
    for (int k = 0; k < dim_; k++) { mn[k] = INFINITY; mx[k] = -INFINITY; }
    for (const auto& b : bodies_) {
        if (!b.alive) continue;
        for (int k = 0; k < dim_; k++) {
            mn[k] = std::min(mn[k], b.pos[k]);
            mx[k] = std::max(mx[k], b.pos[k]);
        }
    }
    if (!std::isfinite(mn[0])) return;  // no alive bodies
    double ctr[8], half[8];
    for (int k = 0; k < dim_; k++) {
        ctr[k] = 0.5 * (mn[k] + mx[k]);
        half[k] = 0.5 * (mx[k] - mn[k]) * 1.000001 + 1e-12;
        if (!(half[k] > 0)) half[k] = 1e-6;
    }
    root_ = newNode(-1, ctr, half);
    for (int i = 0; i < (int)bodies_.size(); i++)
        if (bodies_[i].alive) insertBody(root_, i, 0);

    // Bottom-up: mass / centre of mass / max radius.
    // Children always have higher indices than their parent, so a reverse
    // sweep finalises children before their contributions reach the parent.
    // INVARIANT: during the sweep nd.com holds the RAW WEIGHTED SUM
    // sum(pos_j * m_j) over the subtree (not the average).  Each node pushes
    // that raw sum to its parent as-is; the final pass divides by mass.
    for (int n = (int)nodes_.size() - 1; n >= 0; n--) {
        TreeNode& nd = nodes_[n];
        if (nd.state == 1) {
            const Body& b = bodies_[nd.body];
            nd.mass = b.mass;
            for (int k = 0; k < dim_; k++) nd.com[k] = b.pos[k] * b.mass;
            nd.maxRad = b.radius;
        } else if (nd.state == 3) {
            double m = 0.0;
            double c[8] = {0};
            auto add = [&](int j) {
                const Body& b = bodies_[j];
                m += b.mass;
                for (int k = 0; k < dim_; k++) c[k] += b.pos[k] * b.mass;
                nd.maxRad = std::max(nd.maxRad, b.radius);
            };
            if (nd.body >= 0) add(nd.body);
            for (int j : nd.extra) add(j);
            nd.mass = m;
            for (int k = 0; k < dim_; k++) nd.com[k] = c[k];
        }
        if (n != root_ && nd.mass > 0) {
            TreeNode& p = nodes_[nd.parent];
            for (int k = 0; k < dim_; k++) p.com[k] += nd.com[k];
            p.mass += nd.mass;
            p.maxRad = std::max(p.maxRad, nd.maxRad);
        } else if (n != root_) {
            TreeNode& p = nodes_[nd.parent];
            p.maxRad = std::max(p.maxRad, nd.maxRad);
        }
    }
    // Finalise: divide accumulated com sums, compute metric radius bounds.
    const int nc = 1 << dim_;
    for (auto& nd : nodes_) {
        if (nd.mass > 0) {
            for (int k = 0; k < dim_; k++) nd.com[k] /= nd.mass;
        } else {
            for (int k = 0; k < dim_; k++) nd.com[k] = nd.ctr[k];
        }
        double mr = 0.0;
        for (int c = 0; c < nc; c++) {
            double dd[8] = {0};
            for (int k = 0; k < dim_; k++)
                dd[k] = (nd.ctr[k] + (((c >> k) & 1) ? nd.half[k] : -nd.half[k])) - nd.com[k];
            mr = std::max(mr, metricDistance(dd));
        }
        nd.mrad = mr;
    }
}

// Exact lower bound of the metric distance from xi to any point of the cell.
// Valid for the whole separable family d = F(sum f(|dx_k|)): the objective
// sum f(|dx_k|) is minimised by minimising each |dx_k| independently, i.e. by
// clamping to the axis intervals.
double Simulation::minMetricDistToBox(const double* xi, const TreeNode& nd) const {
    double dd[8] = {0};
    for (int k = 0; k < dim_; k++) {
        double lo = nd.ctr[k] - nd.half[k], hi = nd.ctr[k] + nd.half[k];
        double t = xi[k] < lo ? lo : (xi[k] > hi ? hi : xi[k]);
        dd[k] = t - xi[k];
    }
    return metricDistance(dd);
}

// Central pair force: F > 0 pulls body i TOWARDS the source at posj.
//   a_i += m_j * s * (posj - x_i) / |posj - x_i|
// with s = F(d) outside the core, F(R)*d/R inside (linear regularisation).
void Simulation::applyOneSided(int i, const double* posj, double mj, double radj) {
    double delta[8];
    double rho2 = 0.0;
    for (int k = 0; k < dim_; k++) {
        delta[k] = posj[k] - bodies_[i].pos[k];
        rho2 += delta[k] * delta[k];
    }
    if (!(rho2 > 0.0)) return;  // coincident: no direction, zero force
    if (bodies_[i].fixed) return;  // pinned: exerts force, receives none
    double d = metricDistance(delta);
    if (!(d > 0.0) || !std::isfinite(d)) return;
    const double Rij = bodies_[i].radius + radj;
    double s;
    if (Rij > 0.0 && d < Rij) {
        s = G_ * F_.value(Rij) * (d / Rij);
    } else {
        s = G_ * F_.value(d);
    }
    if (!std::isfinite(s)) return;
    if (std::fabs(s) > kForceClamp) s = s > 0 ? kForceClamp : -kForceClamp;
    if (s == 0.0) return;
    const double w = s / std::sqrt(rho2);  // s per unit Euclidean length
    for (int k = 0; k < dim_; k++)
        bodies_[i].acc[k] += mj * w * delta[k];
}

double Simulation::pairUOneSided(int i, const double* posj, double mj, double radj) const {
    double delta[8];
    for (int k = 0; k < dim_; k++) delta[k] = posj[k] - bodies_[i].pos[k];
    double d = metricDistance(delta);
    if (!std::isfinite(d)) return 0.0;
    const double Rij = bodies_[i].radius + radj;
    double val = pairV(d, Rij);
    if (!std::isfinite(val)) return 0.0;
    return bodies_[i].mass * mj * val;
}

void Simulation::treeForce() {
    for (auto& b : bodies_) std::fill(b.acc.begin(), b.acc.end(), 0.0);
    if (root_ < 0) return;
    int stack[512];
    const int nc = 1 << dim_;
    for (int i = 0; i < (int)bodies_.size(); i++) {
        if (!bodies_[i].alive) continue;
        if (bodies_[i].fixed) continue;  // pinned: no force of its own
        double xi[8];
        for (int k = 0; k < dim_; k++) xi[k] = bodies_[i].pos[k];
        int sp = 0;
        stack[sp++] = root_;
        while (sp > 0) {
            int ni = stack[--sp];
            const TreeNode& nd = nodes_[ni];
            if (nd.mass <= 0.0) continue;
            if (nd.state == 1) {
                if (nd.body != i)
                    applyOneSided(i, bodies_[nd.body].pos.data(), bodies_[nd.body].mass,
                                  bodies_[nd.body].radius);
            } else if (nd.state == 3) {
                if (nd.body >= 0 && nd.body != i)
                    applyOneSided(i, bodies_[nd.body].pos.data(), bodies_[nd.body].mass,
                                  bodies_[nd.body].radius);
                for (int j : nd.extra)
                    if (j != i)
                        applyOneSided(i, bodies_[j].pos.data(), bodies_[j].mass,
                                      bodies_[j].radius);
            } else {
                double delta[8];
                for (int k = 0; k < dim_; k++) delta[k] = nd.com[k] - xi[k];
                double d = metricDistance(delta);
                if (d > 0.0 && nd.mrad < theta_ * d) {
                    applyOneSided(i, nd.com, nd.mass, 0.0);
                } else {
                    for (int k = 0; k < nc; k++) {
                        int c = nd.child[k];
                        if (c >= 0 && nodes_[c].mass > 0 && sp < 511) stack[sp++] = c;
                    }
                }
            }
        }
    }
}

double Simulation::treePotential() const {
    if (root_ < 0) return 0.0;
    // Potential tolerates a far tighter opening criterion than force:
    // at theta_PE = 0.25 the estimator error is ~0.3% (vs -7% at 0.6),
    // and PE is computed once per frame, not once per substep.
    const double thPE = std::min(theta_, 0.25);
    double U = 0.0;
    int stack[512];
    const int nc = 1 << dim_;
    for (int i = 0; i < (int)bodies_.size(); i++) {
        if (!bodies_[i].alive) continue;
        double xi[8];
        for (int k = 0; k < dim_; k++) xi[k] = bodies_[i].pos[k];
        double Ui = 0.0;
        int sp = 0;
        stack[sp++] = root_;
        while (sp > 0) {
            int ni = stack[--sp];
            const TreeNode& nd = nodes_[ni];
            if (nd.mass <= 0.0) continue;
            if (nd.state == 1) {
                if (nd.body != i)
                    Ui += pairUOneSided(i, bodies_[nd.body].pos.data(), bodies_[nd.body].mass,
                                        bodies_[nd.body].radius);
            } else if (nd.state == 3) {
                if (nd.body >= 0 && nd.body != i)
                    Ui += pairUOneSided(i, bodies_[nd.body].pos.data(), bodies_[nd.body].mass,
                                        bodies_[nd.body].radius);
                for (int j : nd.extra)
                    if (j != i)
                        Ui += pairUOneSided(i, bodies_[j].pos.data(), bodies_[j].mass,
                                            bodies_[j].radius);
            } else {
                double delta[8];
                for (int k = 0; k < dim_; k++) delta[k] = nd.com[k] - xi[k];
                double d = metricDistance(delta);
                if (d > 0.0 && nd.mrad < thPE * d) {
                    Ui += pairUOneSided(i, nd.com, nd.mass, 0.0);
                } else {
                    for (int k = 0; k < nc; k++) {
                        int c = nd.child[k];
                        if (c >= 0 && nodes_[c].mass > 0 && sp < 511) stack[sp++] = c;
                    }
                }
            }
        }
        U += Ui;
    }
    return 0.5 * U;
}

void Simulation::collectContacts(std::vector<std::pair<int, int>>& out) const {
    if (root_ < 0) return;
    int stack[512];
    const int nc = 1 << dim_;
    for (int i = 0; i < (int)bodies_.size(); i++) {
        if (!bodies_[i].alive) continue;
        const double Ri = bodies_[i].radius;
        double xi[8];
        for (int k = 0; k < dim_; k++) xi[k] = bodies_[i].pos[k];
        int sp = 0;
        stack[sp++] = root_;
        while (sp > 0) {
            int ni = stack[--sp];
            const TreeNode& nd = nodes_[ni];
            if (nd.mass <= 0.0) continue;
            if (nd.state == 1 || nd.state == 3) {
                auto test = [&](int j) {
                    if (j <= i || !bodies_[j].alive) return;
                    double delta[8];
                    for (int k = 0; k < dim_; k++) delta[k] = bodies_[j].pos[k] - xi[k];
                    double d = metricDistance(delta);
                    double Rij = Ri + bodies_[j].radius;
                    if (std::isfinite(d) && d <= Rij) out.emplace_back(i, j);
                };
                if (nd.state == 1) test(nd.body);
                else {
                    if (nd.body >= 0) test(nd.body);
                    for (int j : nd.extra) test(j);
                }
            } else {
                if (minMetricDistToBox(xi, nd) > Ri + nd.maxRad) continue;  // prune
                for (int k = 0; k < nc; k++) {
                    int c = nd.child[k];
                    if (c >= 0 && nodes_[c].mass > 0 && sp < 511) stack[sp++] = c;
                }
            }
        }
    }
}

int Simulation::addBody(double mass, const double* pos, const double* vel, bool fixed) {
    if (!ok_ || mass <= 0 || !std::isfinite(mass)) return -1;
    Body b;
    b.mass = mass;
    b.radius = radiusFromMass(mass);
    b.pos.assign(pos, pos + dim_);
    b.vel.assign(vel, vel + dim_);
    b.acc.assign(dim_, 0.0);
    b.fixed = fixed;
    if (fixed) std::fill(b.vel.begin(), b.vel.end(), 0.0);  // pinned: no motion
    for (double x : b.pos) if (!std::isfinite(x)) return -1;
    for (double x : b.vel) if (!std::isfinite(x)) return -1;
    bodies_.push_back(std::move(b));
    accValid_ = false;
    return (int)bodies_.size() - 1;
}

void Simulation::clearBodies() {
    bodies_.clear();
    t_ = 0.0;
    accValid_ = false;
}

// ------------------------------------------------------------ metric kernel

// Solve f(x) = S for x >= 0. f is strictly increasing.
double Simulation::invertF(double S) const {
    if (!(S > 0.0)) return 0.0;
    if (fIsPower_) return std::pow(S, 1.0 / fExp_);

    double lo = 0.0, hi = 1.0;
    double fhi = f_.value(hi);
    int expand = 0;
    while (fhi < S) {
        lo = hi;
        hi *= 2.0;
        fhi = f_.value(hi);
        if (++expand > 200) return hi;  // bounded f: distance effectively unbounded
    }
    if (f_.value(lo) > S) return 0.0;

    // Safeguarded Newton (derivative from dual eval).
    double x = 0.5 * (lo + hi);
    for (int it = 0; it < 80; it++) {
        Dual fx = f_.eval(x);
        if (!std::isfinite(fx.v)) break;
        if (std::fabs(fx.v - S) <= 1e-14 * std::max(1.0, S)) return x;
        if (fx.v < S) lo = x; else hi = x;
        if (hi - lo <= 1e-15 * std::max(1.0, hi)) return 0.5 * (lo + hi);
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
        for (int k = 0; k < dim_; k++) S += std::pow(std::fabs(delta[k]), p);
        return S > 0.0 ? std::pow(S, 1.0 / p) : 0.0;
    }
    double S = 0.0;
    for (int k = 0; k < dim_; k++) S += f_.value(std::fabs(delta[k]));
    return invertF(S);
}

bool Simulation::pairKernel(const double* delta, double& dOut, double* grad) const {
    double d;
    if (fIsPower_) {
        const double p = fExp_;
        double S = 0.0;
        for (int k = 0; k < dim_; k++) {
            double e = std::fabs(delta[k]);
            S += std::pow(e, p);
        }
        d = S > 0.0 ? std::pow(S, 1.0 / p) : 0.0;
        if (!(d > 0.0) || !std::isfinite(d)) {
            dOut = d > 0 ? d : 0.0;
            for (int k = 0; k < dim_; k++) grad[k] = 0.0;
            return d > 0.0 && std::isfinite(d);
        }
        // grad[k] = (|dx_k| / d)^(p-1) * sgn(dx_k)
        for (int k = 0; k < dim_; k++) {
            double ratio = std::fabs(delta[k]) / d;  // in [0,1]
            double g = std::pow(ratio, p - 1.0);
            if (!std::isfinite(g)) g = (p < 1.0 && ratio > 0.0) ? kGradClamp : (ratio == 0.0 ? 0.0 : 1.0);
            g *= (delta[k] > 0.0) - (delta[k] < 0.0);
            grad[k] = std::clamp(g, -kGradClamp, kGradClamp);
        }
    } else {
        double S = 0.0;
        double fp[8];
        for (int k = 0; k < dim_; k++) {
            Dual fk = f_.eval(std::fabs(delta[k]));
            S += fk.v;
            fp[k] = fk.d;
        }
        d = invertF(S);
        if (!std::isfinite(d)) {
            dOut = 0.0;
            for (int k = 0; k < dim_; k++) grad[k] = 0.0;
            return false;
        }
        double fpd = f_.eval(d).d;
        if (std::fabs(fpd) < 1e-300) fpd = fpd >= 0 ? 1e-300 : -1e-300;
        for (int k = 0; k < dim_; k++) {
            double g = fp[k] / fpd;
            g *= (delta[k] > 0.0) - (delta[k] < 0.0);
            grad[k] = std::clamp(g, -kGradClamp, kGradClamp);
        }
    }
    dOut = d;
    return true;
}

// ------------------------------------------------------------------ physics

void Simulation::computeAccel() {
    if (treeActive()) {
        buildTree();
        treeForce();
        accValid_ = true;
        return;
    }
    const int n = (int)bodies_.size();
    for (auto& b : bodies_) std::fill(b.acc.begin(), b.acc.end(), 0.0);
    double delta[8];
    for (int i = 0; i < n; i++) {
        if (!bodies_[i].alive) continue;
        for (int j = i + 1; j < n; j++) {
            if (!bodies_[j].alive) continue;
            double rho2 = 0.0;
            for (int k = 0; k < dim_; k++) {
                delta[k] = bodies_[j].pos[k] - bodies_[i].pos[k];
                rho2 += delta[k] * delta[k];
            }
            if (!(rho2 > 0.0)) continue;  // coincident
            double d = metricDistance(delta);
            if (!(d > 0.0) || !std::isfinite(d)) continue;
            const double Rij = bodies_[i].radius + bodies_[j].radius;
            // Force scale: G*F(d) outside the core; inside the core the
            // linear regularisation G*F(R)*(d/R) (C^0 in the force at d = R).
            double s;
            if (Rij > 0.0 && d < Rij) {
                s = G_ * F_.value(Rij) * (d / Rij);
            } else {
                s = G_ * F_.value(d);
            }
            if (!std::isfinite(s)) continue;
            if (std::fabs(s) > kForceClamp) s = s > 0 ? kForceClamp : -kForceClamp;
            if (s == 0.0) continue;
            const double mi = bodies_[i].mass, mj = bodies_[j].mass;
            const double w = s / std::sqrt(rho2);
            // central force: equal and opposite, along the separating line
            // (pinned bodies exert force but receive none)
            for (int k = 0; k < dim_; k++) {
                double f = w * delta[k];
                if (!bodies_[i].fixed) bodies_[i].acc[k] += mj * f;
                if (!bodies_[j].fixed) bodies_[j].acc[k] -= mi * f;
            }
        }
    }
}

void Simulation::resolveContact(int i, int j) {
    Body& A = bodies_[i];
    Body& B = bodies_[j];
    if (!A.alive || !B.alive) return;
    double delta[8];
    for (int k = 0; k < dim_; k++) delta[k] = B.pos[k] - A.pos[k];
    double d = metricDistance(delta);
    if (!std::isfinite(d)) return;
    const double Rij = A.radius + B.radius;
    if (d > Rij) return;

    if (collisionMode_ == Collision::Merge) {
        // A pinned body never takes part in a merger: the halo must stay
        // intact, so the pair simply passes through (core-regularised).
        if (A.fixed || B.fixed) return;
        // Merge j into i (inelastic, momentum-conserving).
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
        // Fully elastic impulse along the contact normal.  Consistent with
        // the central-force model, the normal is the line joining the two
        // centres (identical to the metric-gradient normal for f = t^2).
        double gn2 = 0.0;
        for (int k = 0; k < dim_; k++) gn2 += delta[k] * delta[k];
        double gn = std::sqrt(gn2);
        if (gn < 1e-300) return;
        double nrm[8];
        double w = 0.0;  // closing speed along the normal
        for (int k = 0; k < dim_; k++) {
            nrm[k] = delta[k] / gn;  // A -> B
            w += (A.vel[k] - B.vel[k]) * nrm[k];
        }
        if (w <= 0.0) return;  // already separating
        // Impulse in inverse-mass form: a pinned body has zero inverse mass,
        // i.e. infinite inertia, and reflects the other one elastically.
        const double invA = A.fixed ? 0.0 : 1.0 / A.mass;
        const double invB = B.fixed ? 0.0 : 1.0 / B.mass;
        const double tot = invA + invB;
        if (tot <= 0.0) return;  // both pinned
        const double e = 1.0;  // perfectly elastic
        for (int k = 0; k < dim_; k++) {
            A.vel[k] -= (1.0 + e) * (invA / tot) * w * nrm[k];
            B.vel[k] += (1.0 + e) * (invB / tot) * w * nrm[k];
        }
    }
}

void Simulation::handleCollisions() {
    if (collisionMode_ == Collision::Pass) return;
    if (treeActive()) {
        // The tree from this substep's computeAccel() is still current
        // (kicks do not move positions); merges inside resolveContact()
        // re-execute computeAccel() themselves.
        std::vector<std::pair<int, int>> pairs;
        collectContacts(pairs);
        for (const auto& pr : pairs) resolveContact(pr.first, pr.second);
        return;
    }
    const int n = (int)bodies_.size();
    for (int i = 0; i < n; i++) {
        if (!bodies_[i].alive) continue;
        for (int j = i + 1; j < n; j++) {
            if (!bodies_[j].alive) continue;
            resolveContact(i, j);
        }
    }
}

bool Simulation::step(double dt, int substeps) {
    if (!ok_) return false;
    if (substeps < 1) substeps = 1;
    if (!(dt > 0.0) || !std::isfinite(dt)) return false;
    if (!accValid_) computeAccel();
    const double h = dt;
    for (int s = 0; s < substeps; s++) {
        // kick (half) + drift
        for (auto& b : bodies_) {
            if (!b.alive || b.fixed) continue;  // pinned bodies never move
            for (int k = 0; k < dim_; k++) {
                b.vel[k] += 0.5 * h * b.acc[k];
                b.pos[k] += b.vel[k] * h;
            }
        }
        // new accelerations
        computeAccel();
        // kick (half)
        for (auto& b : bodies_) {
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
    for (auto& b : bodies_)
        if (b.alive) c++;
    return c;
}

double Simulation::potentialEnergy() const {
    if (treeActive()) {
        buildTree();
        return treePotential();
    }
    double U = 0.0;
    const int n = (int)bodies_.size();
    double delta[8];
    for (int i = 0; i < n; i++) {
        if (!bodies_[i].alive) continue;
        for (int j = i + 1; j < n; j++) {
            if (!bodies_[j].alive) continue;
            for (int k = 0; k < dim_; k++) delta[k] = bodies_[j].pos[k] - bodies_[i].pos[k];
            double d = metricDistance(delta);
            if (!std::isfinite(d)) continue;
            const double Rij = bodies_[i].radius + bodies_[j].radius;
            double val = pairV(d, Rij);
            if (std::isfinite(val)) U += bodies_[i].mass * bodies_[j].mass * val;
        }
    }
    return U;
}

double Simulation::energy() const {
    double KE = 0.0;
    for (auto& b : bodies_) {
        if (!b.alive) continue;
        double v2 = 0.0;
        for (int k = 0; k < dim_; k++) v2 += b.vel[k] * b.vel[k];
        KE += 0.5 * b.mass * v2;
    }
    return KE + potentialEnergy();
}

int Simulation::getState(double* out, int maxBodies) const {
    int w = 0;
    for (auto& b : bodies_) {
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

void Simulation::samplePotential(int n, double rmin, double rmax, double* outR, double* outV) const {
    if (n < 1) return;
    if (!(rmin > 0.0)) rmin = 1e-6;
    if (rmax <= rmin) rmax = rmin * 1e6;
    for (int i = 0; i < n; i++) {
        double frac = (n == 1) ? 0.0 : (double)i / (n - 1);
        double r = std::exp(std::log(rmin) + frac * (std::log(rmax) - std::log(rmin)));
        double v = potentialAt(r);  // derived potential V = integral of F
        outR[i] = r;
        outV[i] = std::isfinite(v) ? v : NAN;
    }
}

int Simulation::unitBall2D(int n, double* outXY) const {
    if (n < 3) return 0;
    for (int i = 0; i < n; i++) {
        double th = 2.0 * kPi * i / n;
        double ux = std::cos(th), uy = std::sin(th);
        // Solve metricDistance(rho * u) == 1 for rho by bisection.
        double lo = 0.0, hi = 1.0;
        int expand = 0;
        auto dist = [&](double rho) {
            double d[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            d[0] = rho * ux;
            d[1] = rho * uy;
            return metricDistance(d);
        };
        while (dist(hi) < 1.0 && expand < 60) {
            hi *= 2.0;
            expand++;
        }
        if (dist(hi) < 1.0) {  // unreachable (bounded f)
            outXY[2 * i] = NAN;
            outXY[2 * i + 1] = NAN;
            continue;
        }
        for (int it = 0; it < 60; it++) {
            double mid = 0.5 * (lo + hi);
            if (dist(mid) < 1.0) lo = mid; else hi = mid;
        }
        double rho = 0.5 * (lo + hi);
        outXY[2 * i] = rho * ux;
        outXY[2 * i + 1] = rho * uy;
    }
    return n;
}

}  // namespace grav
