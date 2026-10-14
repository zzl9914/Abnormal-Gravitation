#pragma once
// Gravity simulation with fully user-defined metric and FORCE LAW.
//
// Metric:   d = F_m( sum_k f(|dx_k|) ),  F_m = f^-1,  f strictly increasing
//           on [0,inf), f(0)=0
// Force law (user expression, variable "r" = the metric distance d):
//   pair force on body i from body j:
//       a_i = m_j * F(d_ij) * (x_j - x_i) / |x_j - x_i|
//   F > 0 means ATTRACTION (the force points at the exerting body),
//   F < 0 means repulsion.  The force is always central (along the line
//   joining the two bodies), so total angular momentum is conserved
//   EXACTLY by velocity Verlet for every metric.
//
// Potential (derived, used for the energy monitor, the core regularisation
// and the Barnes-Hut PE estimator):  V(r) = integral of F dr  (V' = F).
// If F is a sum of power terms the antiderivative is exact and analytic;
// otherwise V is tabulated once (log-spaced, Gauss-Legendre quadrature,
// cubic-Hermite interpolation).  Energy is conserved exactly only for the
// isotropic metric f = t^2 (a central force whose magnitude depends on a
// non-radial distance is not curl-free); for anisotropic metrics the
// angular momentum stays exact while the reported energy drifts.
//
// Body radius from matter density rho:
//   R_i = ( m_i * Gamma(d/2+1) / (rho * pi^(d/2)) )^(1/d)
// (2D: disc area m = rho*pi*R^2; 3D: sphere volume m = rho*(4/3)*pi*R^3.)
//
// Core regularization: for pair distance d < R_ij = R_i + R_j the force
// magnitude is replaced by  F(R_ij) * d / R_ij  (linear from 0 at the
// centre, C^0 in the force and the gradient of the matched quadratic
//   W(x) = V(R) + F(R) * (x^2 - R^2) / (2R),
// so the r->0 singularity of F is never evaluated (only V(R), F(R) are).
//
// Collision behaviour when d <= R_ij:
//   Pass   : no collision handling (smooth regularized pass-through)
//   Merge  : inelastic merge, momentum-conserving
//   Bounce : fully elastic impulse along the contact normal (the line
//            joining the centres, consistent with the central force)
//
// Barnes-Hut acceleration: for large N (>= 1500 alive bodies, auto mode) the
// force / collision search uses a quadtree (2D) / octree (3D).  Far cells are
// collapsed into their centre of mass when  cellMetricRadius < theta * d.
// All distance bounds use the USER metric; the separable structure
// d = F_m(sum f(|dx_k|)) gives an exact lower bound of the distance from a
// point to an axis-aligned box (per-axis clamp), which prunes collision
// search exactly.  Tree interactions approximate: momentum is conserved only
// approximately (monopole), and the energy monitor uses the same tree
// approximation.  Direct O(N^2) summation is used below the threshold or in
// "off" mode.
//
// Integration: velocity Verlet (kick-drift-kick), symplectic (exactly so
// for f = t^2, where the force is the exact gradient of the derived V).

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "expr.hpp"

namespace grav {

enum class Collision : int { Pass = 0, Merge = 1, Bounce = 2 };

// Body flag bits (echoed in getState as element 7).
// Fixed bodies (dark-matter halo) never move, never merge and are skipped by
// the integrator; they still exert their force on everything else, and in
// Bounce mode they act as an immovable wall (elastic reflection).
enum BodyFlag : unsigned { kFlagFixed = 1u };

struct Body {
    std::vector<double> pos, vel, acc;
    double mass = 1.0;
    double radius = 0.0;  // physical radius from density
    bool alive = true;
    bool fixed = false;   // pinned (halo): no motion, no merging
};

// Quad-/octree node. state: 0 = empty, 1 = leaf (one body), 2 = internal,
// 3 = bucket (depth-capped leaf holding several coincident bodies).
struct TreeNode {
    double mass = 0.0;
    double com[8] = {};    // centre of mass
    double ctr[8] = {};    // geometric centre of the cell
    double half[8] = {};   // half extents of the cell
    double mrad = 0.0;     // metric radius bound: max d(com, corner)
    double maxRad = 0.0;   // max physical body radius in subtree
    int child[8];
    int body = -1;
    int parent = -1;
    uint8_t state = 0;
    std::vector<int> extra;  // overflow bodies at max depth
    TreeNode() {
        for (int i = 0; i < 8; i++) child[i] = -1;
    }
};

class Simulation {
public:
    // dims: any >= 1 (2 or 3 in practice). `forceExpr` is the force law F(r)
    // (F > 0 attractive). Returns with ok()==false + error() if the
    // expressions are invalid or f is not a valid metric generator.
    Simulation(int dims, const std::string& fExpr, const std::string& forceExpr);

    bool ok() const { return ok_; }
    const std::string& error() const { return err_; }

    void setCollision(Collision mode, double density);
    // Coupling constant of the force law: the pair force and the derived
    // potential are scaled by G (G > 0, default 1).  a_i = G * m_j * F(d) * u.
    // Used by the frontend to calibrate the evolution timescale: with
    // G * F(R/2) held fixed across configurations, the force at the
    // half-radius of a random cloud is the same in every setup.
    void setG(double G);
    // Raw force-law value F(r) at metric distance r (no G factor).
    double forceAt(double r) const;
    // treeMode: 0 = auto (use tree when alive >= 1500), 1 = always direct.
    // theta: Barnes-Hut opening angle (smaller = more accurate).
    void setTree(int treeMode, double theta);

    // `fixed` bodies are pinned: they exert force but never move, never merge
    // and reflect others elastically (Bounce) as an immovable wall.
    int  addBody(double mass, const double* pos, const double* vel, bool fixed = false);
    void clearBodies();  // also resets t to 0

    // `substeps` velocity-Verlet steps of size dt. Returns false on fatal error.
    bool step(double dt, int substeps);

    double time() const { return t_; }
    double energy() const;  // KE + PE, recomputed at the end of each step()
    int    aliveCount() const;
    int    dims() const { return dim_; }

    // Packed per alive body: [mass, x, y, z, vx, vy, vz, flags] (2D pads
    // zeros; flags = BodyFlag bits).  Returns number of bodies written.
    int getState(double* out, int maxBodies) const;

    // Analysis helpers (used by the web UI).
    // n log-spaced samples of the DERIVED potential V(r) on [rmin, rmax].
    void samplePotential(int n, double rmin, double rmax, double* outR, double* outV) const;
    // Unit ball of the metric, sliced on the first two axes: n points (x,y)
    // with metric distance 1 from the origin. Points with unreachable radius
    // are emitted as (NAN, NAN).
    int unitBall2D(int n, double* outXY) const;

    // Exposed for the test harness.
    double metricDistance(const double* delta) const;
    // d and grad(d) w.r.t. delta. grad[k] = d(d)/d(delta_k).
    bool pairKernel(const double* delta, double& dOut, double* grad) const;
    // Derived potential V(r) (test harness). r must be > 0.
    double potentialAt(double r) const;
    // Potential energy only (used by finite-difference force tests).
    double potentialEnergy() const;
    // Physical radius of a body of mass m under the current density.
    double radiusFromMass(double m) const;

private:
    void computeAccel();
    void handleCollisions();
    double invertF(double S) const;
    // pair potential incl. core regularisation at metric distance d, R=R_i+R_j
    double pairV(double d, double R) const;
    void buildPotentialTable();

    // ---- Barnes-Hut tree ----
    static const int kMaxDepth = 48;
    static const int kTreeThreshold = 1500;
    bool treeActive() const;
    int  newNode(int parent, const double* ctr, const double* half) const;
    void insertBody(int ni, int b, int depth) const;
    void buildTree() const;
    void treeForce();
    double treePotential() const;
    // one-sided interaction: add force on body i from a "source" at posj
    void applyOneSided(int i, const double* posj, double mj, double radj);
    double pairUOneSided(int i, const double* posj, double mj, double radj) const;
    // tree-pruned collection of overlapping (i<j) pairs for collisions
    void collectContacts(std::vector<std::pair<int, int>>& out) const;
    double minMetricDistToBox(const double* xi, const TreeNode& nd) const;
    // check one pair for contact and resolve it (merge or bounce)
    void resolveContact(int i, int j);

    int dim_ = 2;
    Expr f_, F_;   // metric generator f(t); force law F(r)
    bool ok_ = false;
    std::string err_;
    bool fIsPower_ = false;
    double fExp_ = 1.0;  // f(t) = t^fExp_
    // Derived potential V(r) = sum_i a_i r^{q_i} + vLn_*ln(r)  (analytic)
    // or a sampled table (numeric fallback).
    std::vector<std::pair<double, double>> vPow_;  // (q_i, a_i)
    double vLn_ = 0.0;
    bool vAnalytic_ = false;
    std::vector<double> vTabR_, vTabV_, vTabF_;  // log-spaced r, V, F=V'
    std::vector<Body> bodies_;
    double t_ = 0.0;
    Collision collisionMode_ = Collision::Pass;
    double density_ = 100.0;
    double G_ = 1.0;       // force-law coupling constant (scales F and V)
    int treeMode_ = 0;      // 0 auto, 1 off
    double theta_ = 0.6;
    mutable std::vector<TreeNode> nodes_;
    mutable int root_ = -1;
    bool accValid_ = false;  // acc up to date with positions/masses?
};

}  // namespace grav
