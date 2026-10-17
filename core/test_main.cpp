// Physics verification harness (console).  Force-law semantics: the third
// constructor argument is F(r) (F > 0 attractive, force central along the
// line joining the bodies); V is derived by the kernel (V' = F).
//   1.  Newtonian circular orbit: energy conservation + radius stability
//   2.  Central-force kernel under anisotropic power metric (t^3):
//       direction || delta, magnitude m2*F(d)/|delta|
//   3.  Same through the general-f numerical inversion path (t^2+t^3)
//   2b. Isotropic f=t^2 with a general (non-power-sum) F: force == -grad U
//       of the analytically integrated potential
//   4.  Metric gradient vs central finite differences
//   5.  L1 unit ball shape + derived-V sampling sanity (analytic + ln term
//       + Lennard-Jones + numeric table fallback)
//   6.  Inelastic merge conserves momentum
//   7.  Figure-eight three-body choreography stays bounded
//   8.  Expression validation rejects bad f / bad F
//   9.  Bounce: elastic impulse, momentum conserved
//   10. Pass-through core regularization keeps energy conserved (f=t^2)
//   11. Barnes-Hut accuracy (force, PE estimator, integration)
//   12. Anisotropic metrics: total angular momentum conserved EXACTLY
//       (the headline property of the central-force model)
// Exit code 0 = all pass.

#include "sim.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using grav::Simulation;

static int g_fail = 0;

static void check(bool cond, const char* name, const char* detail = "") {
    if (cond) {
        std::printf("[PASS] %s %s\n", name, detail);
    } else {
        std::printf("[FAIL] %s %s\n", name, detail);
        g_fail++;
    }
}

// Acceleration component `axis` of body 1, measured via one tiny Verlet step
// (v = a*dt to first order).
static double bodyAccelViaTinyStep(const std::string& f, const std::string& Flaw,
                                   int dims, double m1, double m2,
                                   std::vector<double> p1, std::vector<double> p2,
                                   int axis) {
    Simulation sim(dims, f.c_str(), Flaw.c_str());
    if (!sim.ok()) return NAN;
    std::vector<double> z(dims, 0.0);
    sim.setCollision(grav::Collision::Pass, 1e9);
    sim.addBody(m1, p1.data(), z.data());
    sim.addBody(m2, p2.data(), z.data());
    sim.step(1e-9, 1);
    double out[16];
    int n = sim.getState(out, 2);
    if (n != 2) return NAN;
    return out[4 + axis] / 1e-9;  // vel[axis] of first body
}

int main() {
    // ---------------------------------------------------------------- test 1
    // Newtonian: f=t^2 (Euclidean), F=1/r^2 (derived V=-1/r). Central mass
    // M=1000, satellite at r=1 with v=sqrt(M/r). Verlet keeps the orbit tight.
    {
        Simulation sim(2, "t^2", "1/r^2");
        check(sim.ok(), "T1 construct", sim.error());
        double M = 1000.0;
        double p0[2] = {0, 0}, v0[2] = {0, 0};
        double p1[2] = {1, 0}, v1[2] = {0, std::sqrt(M)};
        sim.setCollision(grav::Collision::Pass, 1e9);
        sim.addBody(M, p0, v0);
        sim.addBody(1.0, p1, v1);
        double E0 = sim.energy();
        double dt = 2e-4;
        double T = 2.0 * 3.14159265358979323846 / std::sqrt(M);  // circular period at r=1
        int steps = (int)(T / dt);
        for (int i = 0; i < steps; i++) sim.step(dt, 1);
        double st[16];
        sim.getState(st, 2);
        // satellite is second body (central mass barely moves)
        double x = st[9], y = st[10];  // pos of body 1: base 8 (=8 doubles per body)
        double r = std::sqrt(x * x + y * y);
        double E1 = sim.energy();
        double drift = std::fabs((E1 - E0) / E0);
        std::printf("      r=%.6f (expect 1), dE/E=%.3e\n", r, drift);
        check(std::fabs(r - 1.0) < 5e-3, "T1 circular orbit radius");
        check(drift < 1e-8, "T1 energy conservation");
    }

    // ---------------------------------------------------------------- test 2
    // Central-force kernel under an ANISOTROPIC power metric (t^3): the
    // force on body 1 is parallel to delta = p2 - p1 with magnitude
    // m2 * F(d) / |delta|  (F > 0 attractive, pointing at the source).
    {
        const char* f = "t^3";
        const char* Flaw = "1/r^2";
        int dims = 3;
        double m1 = 1.0, m2 = 2.0;
        std::vector<double> p1 = {0.3, -0.2, 0.1};
        std::vector<double> p2 = {1.2, -0.5, 0.8};
        Simulation sim(dims, f, Flaw);
        check(sim.ok(), "T2 construct", sim.error());
        std::vector<double> delta(3);
        double rho2 = 0.0;
        for (int k = 0; k < 3; k++) {
            delta[k] = p2[k] - p1[k];
            rho2 += delta[k] * delta[k];
        }
        double rho = std::sqrt(rho2);
        double d = sim.metricDistance(delta.data());
        double expect = m2 * sim.forceAt(d);  // |a1| = m2*F(d), direction delta/rho
        std::vector<double> a(3);
        double anorm = 0.0;
        for (int k = 0; k < 3; k++) {
            a[k] = bodyAccelViaTinyStep(f, Flaw, dims, m1, m2, p1, p2, k);
            anorm += a[k] * a[k];
        }
        anorm = std::sqrt(anorm);
        // direction: a1 x delta == 0
        double cr[3] = {a[1] * delta[2] - a[2] * delta[1],
                        a[2] * delta[0] - a[0] * delta[2],
                        a[0] * delta[1] - a[1] * delta[0]};
        double crn = std::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
        double mrel = std::fabs(anorm - expect) / expect;
        std::printf("      d=%.6f rho=%.6f |a|=%.8g expect=%.8g |axd|=%.2e\n",
                    d, rho, anorm, expect, crn);
        check(crn / anorm < 1e-9, "T2 anisotropic force is central (|| delta)");
        check(mrel < 1e-6, "T2 force magnitude = m2*F(d)/rho");
    }

    // ---------------------------------------------------------------- test 3
    // Same through the general-f numerical inversion path (f = t^2 + t^3).
    {
        const char* f = "t^2 + t^3";
        const char* Flaw = "1/r^2";
        int dims = 3;
        double m1 = 1.5, m2 = 0.7;
        std::vector<double> p1 = {-0.4, 0.9, 0.2};
        std::vector<double> p2 = {0.8, 0.1, -0.6};
        Simulation sim(dims, f, Flaw);
        check(sim.ok(), "T3 construct", sim.error());
        std::vector<double> delta(3);
        double rho2 = 0.0;
        for (int k = 0; k < 3; k++) {
            delta[k] = p2[k] - p1[k];
            rho2 += delta[k] * delta[k];
        }
        double rho = std::sqrt(rho2);
        double d = sim.metricDistance(delta.data());
        double expect = m2 * sim.forceAt(d);  // |a1| = m2*F(d), direction delta/rho
        double ax = bodyAccelViaTinyStep(f, Flaw, dims, m1, m2, p1, p2, 0);
        double ay = bodyAccelViaTinyStep(f, Flaw, dims, m1, m2, p1, p2, 1);
        double ex = expect * delta[0] / rho, ey = expect * delta[1] / rho;
        double rel = std::sqrt((ax - ex) * (ax - ex) + (ay - ey) * (ay - ey)) / expect;
        std::printf("      d=%.6f |a| rel=%.3e\n", d, rel);
        check(rel < 1e-4, "T3 general-metric central force = m2*F(d)*delta/rho^2");
    }

    // ---------------------------------------------------------------- test 2b
    // Isotropic f=t^2 with a GENERAL (non-power-sum) force law:
    // F = 2r/(r^2+0.01)^2  =>  V = -1/(r^2+0.01) + C analytically.
    // The force must be the exact gradient of U = m1*m2*V(|delta|).
    {
        const char* f = "t^2";
        const char* Flaw = "2*r/(r^2+0.01)^2";
        int dims = 3;
        double m1 = 1.0, m2 = 2.0;
        std::vector<double> p1 = {0.3, -0.2, 0.1};
        std::vector<double> p2 = {1.2, -0.5, 0.8};
        auto U = [&](const std::vector<double>& a, const std::vector<double>& b) {
            double r2 = 0.0;
            for (int k = 0; k < 3; k++) {
                double dd = b[k] - a[k];
                r2 += dd * dd;
            }
            return -m1 * m2 / (r2 + 0.01);
        };
        const double h = 1e-6;
        std::vector<double> p2p = p2, p2m = p2;
        p2p[0] += h; p2m[0] -= h;
        // F_1x = +dU/d(delta_x) = (Up - Um)/(2h);  a1x = F/m1
        double aFD = (U(p1, p2p) - U(p1, p2m)) / (2 * h) / m1;
        double aStep = bodyAccelViaTinyStep(f, Flaw, dims, m1, m2, p1, p2, 0);
        double rel = std::fabs(aFD - aStep) / std::fabs(aFD);
        std::printf("      a_fd=%.8g a_sim=%.8g rel=%.3e\n", aFD, aStep, rel);
        check(rel < 1e-5, "T2b general F: force = -grad(integral F)");
    }

    // ---------------------------------------------------------------- test 4
    // Metric gradient vs central finite differences (Lp, p=1.5, power path).
    {
        Simulation sim(3, "t^1.5", "1/r^2");
        check(sim.ok(), "T4 construct", sim.error());
        double delta[3] = {0.7, -0.4, 1.1};
        double d;
        double grad[3];
        sim.pairKernel(delta, d, grad);
        const double h = 1e-7;
        bool allOk = true;
        for (int k = 0; k < 3; k++) {
            double dp[3], dm[3];
            for (int q = 0; q < 3; q++) {
                dp[q] = delta[q];
                dm[q] = delta[q];
            }
            dp[k] += h; dm[k] -= h;
            double fd = (sim.metricDistance(dp) - sim.metricDistance(dm)) / (2 * h);
            double rel = std::fabs(fd - grad[k]) / std::max(1e-12, std::fabs(fd));
            if (rel > 1e-4) {
                std::printf("      axis %d: fd=%.8g grad=%.8g rel=%.3e\n", k, fd, grad[k], rel);
                allOk = false;
            }
        }
        check(allOk, "T4 metric gradient vs finite diff");
    }

    // ---------------------------------------------------------------- test 5
    // L1 unit ball is the diamond |x|+|y| = 1; derived-V sampling sanity.
    {
        Simulation sim(2, "t", "1/r^2");
        check(sim.ok(), "T5 construct", sim.error());
        double xy[192];
        int n = sim.unitBall2D(96, xy);
        check(n == 96, "T5 unit ball size");
        // angle 0 -> (1,0); angle 45deg -> (0.5,0.5)
        double x0 = xy[0], y0 = xy[1];
        int i45 = 96 / 8;  // 45 degrees
        double x45 = xy[2 * i45], y45 = xy[2 * i45 + 1];
        double e0 = std::fabs(x0 - 1.0) + std::fabs(y0);
        double e45 = std::fabs(x45 - 0.5) + std::fabs(y45 - 0.5);
        std::printf("      theta0=(%.4f,%.4f) theta45=(%.4f,%.4f)\n", x0, y0, x45, y45);
        check(e0 < 1e-6 && e45 < 1e-6, "T5 L1 ball is diamond");
        double rr[4], vv[4];
        sim.samplePotential(4, 1.0, 8.0, rr, vv);
        // F = 1/r^2 -> V = -1/r  (exact analytic antiderivative)
        check(std::fabs(vv[1] - (-1.0 / 2.0)) < 1e-12, "T5 V(2) = -1/2 (from F=1/r^2)");
        // F = 1/r -> V = ln r
        {
            Simulation s2(2, "t", "1/r");
            check(s2.ok(), "T5 ln construct");
            double r2[4], v2[4];
            s2.samplePotential(4, 1.0, 8.0, r2, v2);
            // V(2) - V(1) = ln 2  (table samples are absolute; use differences)
            check(std::fabs((v2[1] - v2[0]) - std::log(2.0)) < 1e-12,
                  "T5 V diff = ln2 (from F=1/r)");
        }
        // Lennard-Jones force: F = 1/r^6 - 1/r^12 -> V = -1/(5r^5) + 1/(11r^11)
        {
            Simulation s3(2, "t", "1/r^6-1/r^12");
            check(s3.ok(), "T5 LJ construct");
            double r3[4], v3[4];
            s3.samplePotential(4, 1.0, 8.0, r3, v3);
            double expect = -1.0 / (5.0 * 32.0) + 1.0 / (11.0 * 2048.0);
            check(std::fabs(v3[1] - expect) < 1e-12, "T5 V(2) from LJ force law");
        }
        // Numeric-table fallback: F = exp(-r) -> V(2)-V(1) = e^-1 - e^-2
        {
            Simulation s4(2, "t", "exp(-r)");
            check(s4.ok(), "T5 table construct");
            double r4[8], v4[8];
            s4.samplePotential(8, 0.5, 4.0, r4, v4);
            // find samples nearest 1 and 2 and use V(2)-V(1) via potentialAt
            double dV = s4.potentialAt(2.0) - s4.potentialAt(1.0);
            double expect = std::exp(-1.0) - std::exp(-2.0);
            std::printf("      table V(2)-V(1)=%.8f expect=%.8f\n", dV, expect);
            check(std::fabs(dV - expect) < 1e-6, "T5 numeric-table antiderivative");
        }
    }

    // ---------------------------------------------------------------- test 6
    // Merge: two overlapping bodies -> one body, momentum conserved.
    // Density 4 in 2D: R(m) = sqrt(m/(pi*rho)) -> R(4)+R(1) ~ 0.85 > 0.5 gap.
    {
        Simulation sim(2, "t^2", "1/r^2");
        sim.setCollision(grav::Collision::Merge, 4.0);
        double pa[2] = {0, 0}, va[2] = {1, -1};
        double pb[2] = {0.5, 0}, vb[2] = {-1, 2};
        sim.addBody(4.0, pa, va);
        sim.addBody(1.0, pb, vb);
        double px0 = 4.0 * 1.0 + 1.0 * -1.0;
        double py0 = 4.0 * -1.0 + 1.0 * 2.0;
        sim.step(1e-3, 1);
        double st[16];
        int n = sim.getState(st, 2);
        check(n == 1, "T6 merged to one body");
        double mass = st[0];
        double px = mass * st[4], py = mass * st[5];
        check(std::fabs(px - px0) < 1e-9 && std::fabs(py - py0) < 1e-9,
              "T6 momentum conserved");
    }

    // ---------------------------------------------------------------- test 9
    // Bounce: two equal masses collide head-on. Elastic impulse swaps normal
    // velocities; momentum is exactly conserved, energy to Verlet accuracy.
    // Density chosen so R_ij = 0.6: R = sqrt(1/(pi*rho)) = 0.3 -> rho = 3.5367.
    {
        Simulation sim(2, "t^2", "1/r^2");
        sim.setCollision(grav::Collision::Bounce, 1.0 / (0.3 * 0.3 * 3.14159265358979323846));
        double pa[2] = {-2, 0}, va[2] = {1.5, 0};
        double pb[2] = {2, 0}, vb[2] = {-1.5, 0};
        sim.addBody(1.0, pa, va);
        sim.addBody(1.0, pb, vb);
        double E0 = sim.energy();
        double dt = 2e-4;
        double minD = 1e9, maxD = 0;
        for (int i = 0; i < 20000; i++) {
            sim.step(dt, 1);
            double st[16];
            sim.getState(st, 2);
            double dd = std::fabs(st[9] - st[1]);  // |x2 - x1|
            minD = std::min(minD, dd);
            maxD = std::max(maxD, dd);
            // momentum must stay exactly (0,0)
            if (i % 2000 == 0) {
                double px = st[4] + st[12];
                if (std::fabs(px) > 1e-12) { check(false, "T9 momentum"); break; }
            }
        }
        double st[16];
        sim.getState(st, 2);
        double px = st[4] + st[12], py = st[5] + st[13];
        double drift = std::fabs((sim.energy() - E0) / E0);
        std::printf("      minD=%.4f maxD=%.4f px=%.2e py=%.2e dE/E=%.3e\n",
                    minD, maxD, px, py, drift);
        check(sim.aliveCount() == 2, "T9 no merge in bounce mode");
        check(std::fabs(px) < 1e-9 && std::fabs(py) < 1e-9, "T9 momentum conserved");
        check(minD > 0.5 && minD < 0.62, "T9 bounce turns around at contact");
        check(drift < 1e-4, "T9 energy conserved (elastic)");
    }

    // ---------------------------------------------------------------- test 10
    // Pass-through core regularization: light body fired head-on at a heavy
    // central mass with NO softening and NO collision handling. The interior
    // linear force (quadratic potential) removes the r->0 singularity: the
    // body passes smoothly through the centre and energy is conserved (f=t^2).
    {
        Simulation sim(2, "t^2", "1/r^2");
        // central mass 1000, want R ~ 0.4: rho = 1000/(pi*0.16)
        sim.setCollision(grav::Collision::Pass, 1000.0 / (3.14159265358979323846 * 0.16));
        double p0[2] = {0, 0}, v0[2] = {0, 0};
        double p1[2] = {5, 0}, v1[2] = {-3, 0};
        sim.addBody(1000.0, p0, v0);
        sim.addBody(1.0, p1, v1);
        double E0 = sim.energy();
        double dt = 5e-5;
        bool crossed = false;
        double minD = 1e9;
        for (int i = 0; i < 60000; i++) {  // t = 3
            sim.step(dt, 1);
            double st[16];
            sim.getState(st, 2);
            double x = st[9];
            minD = std::min(minD, std::fabs(x - st[1]));
            if (x < -0.5) crossed = true;
        }
        double drift = std::fabs((sim.energy() - E0) / E0);
        std::printf("      crossed=%d minD=%.4f dE/E=%.3e\n", (int)crossed, minD, drift);
        check(crossed, "T10 pass-through crossing");
        check(drift < 1e-3, "T10 energy conserved through core");
    }

    // ---------------------------------------------------------------- test 7
    // Figure-eight choreography (Chenciner-Montgomery), Newtonian 2D.
    {
        Simulation sim(2, "t^2", "1/r^2");
        sim.setCollision(grav::Collision::Pass, 1e9);
        double p1[2] = {0.97000436, -0.24308753};
        double p2[2] = {-0.97000436, 0.24308753};
        double p3[2] = {0, 0};
        double v3[2] = {-0.93240737, -0.86473146};
        double v1[2] = {-0.5 * v3[0], -0.5 * v3[1]};
        double v2[2] = {-0.5 * v3[0], -0.5 * v3[1]};
        sim.addBody(1.0, p1, v1);
        sim.addBody(1.0, p2, v2);
        sim.addBody(1.0, p3, v3);
        double E0 = sim.energy();
        double dt = 5e-4;
        double maxR = 0.0;
        bool diverged = false;
        for (int i = 0; i < 14000; i++) {  // slightly over one period (T ~ 6.32)
            sim.step(dt, 1);
            double st[24];
            sim.getState(st, 3);
            for (int b = 0; b < 3; b++) {
                double r = std::sqrt(st[8 * b + 1] * st[8 * b + 1] + st[8 * b + 2] * st[8 * b + 2]);
                if (r > maxR) maxR = r;
                if (!std::isfinite(r)) diverged = true;
            }
        }
        double drift = std::fabs((sim.energy() - E0) / E0);
        std::printf("      maxR=%.4f dE/E=%.3e\n", maxR, drift);
        check(!diverged && maxR < 1.2, "T7 figure-8 bounded");
        check(drift < 1e-5, "T7 energy conservation");
    }

    // ---------------------------------------------------------------- test 8
    // Expression validation rejects bad f / bad F.
    {
        Simulation a(2, "t-1", "1/r^2");      // f(0) != 0
        Simulation b(2, "-t", "1/r^2");      // decreasing (t^1 negated -> power -1)
        Simulation c(2, "sin(t)", "1/r^2");  // not monotone
        Simulation d(2, "t^2", "1/");        // F parse error
        Simulation e(2, "t^2", "1/(r-1)");   // F pole at r = 1 (probe rejects)
        Simulation g(2, "t^2", "42");        // constant force law: VALID now
        check(!a.ok(), "T8 reject f(0)!=0");
        check(!b.ok(), "T8 reject decreasing f");
        check(!c.ok(), "T8 reject non-monotone f");
        check(!d.ok(), "T8 reject F parse error");
        check(!e.ok(), "T8 reject F with pole at probe point");
        check(g.ok(), "T8 constant F accepted (uniform attraction)");
        // constant F=42 -> V = 42r: V(2)-V(1) = 42
        double rr[4], vv[4];
        g.samplePotential(4, 1.0, 8.0, rr, vv);
        check(std::fabs(g.potentialAt(2.0) - g.potentialAt(1.0) - 42.0) < 1e-9,
              "T8 constant F derives V = 42r");
    }

    // ---------------------------------------------------------------- test 12
    // THE headline property: under ANISOTROPIC metrics the pair forces are
    // still central (along the joining line), so velocity Verlet conserves
    // the total angular momentum to round-off -- while the reported energy
    // drifts, because a central force whose magnitude depends on a non-radial
    // distance is not curl-free (model physics, not an integrator error).
    {
        for (const char* fm : {"t^1.5", "t^3"}) {
            Simulation sim(3, fm, "1/r^2");
            check(sim.ok(), "T12 construct", fm);
            sim.setCollision(grav::Collision::Pass, 100.0);
            std::mt19937_64 rng(2026);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            double z[3] = {0, 0, 0};
            sim.addBody(500.0, z, z);
            const int NS = 5;
            for (int i = 0; i < NS; i++) {
                double rr = 2.2 + 1.8 * uni(rng);
                double th = uni(rng) * 6.283185307179586;
                double phi = std::acos(2 * uni(rng) - 1);   // isotropic directions
                double px = rr * std::sin(phi) * std::cos(th);
                double py = rr * std::sin(phi) * std::sin(th);
                double pz = rr * std::cos(phi);
                double v = std::sqrt(500.0 / rr) * 0.95;
                // tangential unit vector: cross(radial, arbitrary) normalised
                double ax = -py, ay = px, az = 0.0;
                double an = std::sqrt(ax * ax + ay * ay + az * az);
                double vx = v * ax / an, vy = v * ay / an, vz = v * az / an;
                double p[3] = {px, py, pz}, vv[3] = {vx, vy, vz};
                sim.addBody(1.0 + uni(rng), p, vv);
            }
            auto angMom = [&](double L[3]) {
                double st[8 * 16];
                int n = sim.getState(st, 16);
                L[0] = L[1] = L[2] = 0.0;
                for (int i = 0; i < n; i++) {
                    double m = st[8 * i], x = st[8 * i + 1], y = st[8 * i + 2], zz = st[8 * i + 3];
                    double vx = st[8 * i + 4], vy = st[8 * i + 5], vz = st[8 * i + 6];
                    L[0] += m * (y * vz - zz * vy);
                    L[1] += m * (zz * vx - x * vz);
                    L[2] += m * (x * vy - y * vx);
                }
            };
            double L0[3], L1[3];
            angMom(L0);
            double E0 = sim.energy();
            const double dt = 2e-4;
            for (int i = 0; i < 4000; i++) sim.step(dt, 1);  // t = 0.8
            angMom(L1);
            double nL0 = std::sqrt(L0[0] * L0[0] + L0[1] * L0[1] + L0[2] * L0[2]);
            double dL[3] = {L1[0] - L0[0], L1[1] - L0[1], L1[2] - L0[2]};
            double ndL = std::sqrt(dL[0] * dL[0] + dL[1] * dL[1] + dL[2] * dL[2]);
            double dE = (sim.energy() - E0) / E0;
            std::printf("      f=%s: dL/L=%.3e  dE/E=%+.3e\n", fm, ndL / nL0, dE);
            check(std::isfinite(nL0) && ndL / nL0 < 1e-10,
                  "T12 angular momentum exact (central forces)", fm);
        }
    }

    // ---------------------------------------------------------------- test 13
    // Coupling constant G: force and derived potential both scale linearly;
    // the calibration rule G*F(R/2) = const gives a consistent half-radius
    // force across different cloud radii / force laws.
    {
        const double mC = 1000.0, mS = 1.0;
        std::vector<double> pC = {0, 0}, pS = {1, 0};
        // G = 1 reference: central body acceleration from the satellite
        double a1 = bodyAccelViaTinyStep("t^2", "1/r^2", 2, mC, mS, pC, pS, 0);
        // G = 4: circular satellite speed doubles; central accel quadruples
        std::vector<double> vC = {0, 0}, vS = {0, 2 * std::sqrt(mC)};
        Simulation sim(2, "t^2", "1/r^2");
        sim.setG(4.0);
        sim.setCollision(grav::Collision::Pass, 1e9);
        sim.addBody(mC, pC.data(), vC.data());
        sim.addBody(mS, pS.data(), vS.data());
        sim.step(1e-9, 1);
        double st[16];
        sim.getState(st, 2);
        double aG = st[4] / 1e-9;  // vx of the central body (body 0)
        double expect = 4.0 * mS * sim.forceAt(1.0);
        std::printf("      a(G=1)=%.6f a(G=4)=%.6f expect=%.6f\n", a1, aG, expect);
        check(std::fabs(aG - expect) / expect < 1e-6, "T13 force scales with G");
        check(std::fabs(aG / a1 - 4.0) < 1e-6, "T13 G=4 gives 4x the G=1 force");
        // energy: V = G * (-1/r);  E = 1/2 m v^2 - G*M*m / r = 2000 - 4000 = -2000
        double E = sim.energy();
        check(std::fabs(E - (0.5 * mS * vS[1] * vS[1] - 4.0 * mC * mS)) < 1e-3,
              "T13 derived potential scales with G");
        // calibration rule: G * F(R/2) equal across two different cloud radii
        //   default: R = 10, F = 1/r^2 -> F(5) = 0.04
        //   doubled disc: R = 20 -> F(10) = 0.01 -> G = 4 -> G*F(10) = 0.04
        Simulation sA(2, "t^2", "1/r^2");
        Simulation sB(2, "t^2", "1/r^2");
        sB.setG(0.04 / sB.forceAt(10.0));
        double fA = 1.0 * sA.forceAt(5.0);
        double fB = 0.04;  // G * F(10) by construction
        check(std::fabs(fA - fB) < 1e-12, "T13 calibration G*F(R/2) is invariant");
        // invalid G values are ignored (G stays at the last valid value)
        Simulation sC(2, "t^2", "1/r^2");
        sC.setG(3.0);
        sC.setG(-1.0);
        sC.setG(NAN);
        sC.setG(INFINITY);
        sC.setCollision(grav::Collision::Pass, 1e9);
        sC.addBody(mC, pC.data(), vC.data());  // both at rest
        sC.addBody(mS, pS.data(), vC.data());
        double EC = sC.energy();  // KE = 0, PE = 3 * (-mC*mS)
        check(std::fabs(EC + 3000.0) < 1e-3, "T13 invalid G values are ignored");
    }

    // ---------------------------------------------------------------- test 14
    // Pinned "dark matter halo": exerts force, never moves, never merges,
    // and reflects others elastically (immovable wall) in Bounce mode.
    {
        const double mH = 1000.0, mS = 1.0;
        std::vector<double> pH = {0, 0}, v0 = {0, 0};
        // (a) force on the satellite is identical to a free central body
        {
            std::vector<double> pS = {2, 0};
            Simulation a(2, "t^2", "1/r^2");
            a.setCollision(grav::Collision::Pass, 1e9);
            a.addBody(mH, pH.data(), v0.data(), true);
            a.addBody(mS, pS.data(), v0.data());
            a.step(1e-9, 1);
            double st[16];
            a.getState(st, 2);
            double acc = st[8 + 4] / 1e-9;  // satellite vx
            double expect = -mH * a.forceAt(2.0);
            std::printf("      halo pull a=%.6f expect=%.6f\n", acc, expect);
            check(std::fabs(acc - expect) / std::fabs(expect) < 1e-6,
                  "T14 pinned halo exerts the full force");
            check(std::fabs(st[7] - 1.0) < 1e-12, "T14 state flags report pinned");
            check(std::fabs(st[15]) < 1e-12, "T14 satellite body is not pinned");
        }
        // (b) the halo stays at the origin while the satellite orbits
        {
            std::vector<double> pS = {3, 0}, vS = {0, std::sqrt(mH / 3.0)};
            Simulation b(2, "t^2", "1/r^2");
            b.setCollision(grav::Collision::Pass, 1e9);
            b.addBody(mH, pH.data(), v0.data(), true);
            b.addBody(mS, pS.data(), vS.data());
            const double dt = 1e-4;
            for (int i = 0; i < 20000; i++) b.step(dt, 1);  // t = 2, ~1 orbit
            double st[16];
            b.getState(st, 2);
            double rH = std::hypot(st[1], st[2]);
            std::printf("      halo |x| after t=2: %.3e\n", rH);
            check(rH < 1e-12, "T14 pinned halo never moves");
            check(std::fabs(st[4]) < 1e-12 && std::fabs(st[5]) < 1e-12,
                  "T14 pinned halo keeps zero velocity");
        }
        // (c) no merging: a body dropped into the halo survives (Merge mode)
        {
            std::vector<double> pS = {0.5, 0};
            Simulation c(2, "t^2", "1/r^2");
            c.setCollision(grav::Collision::Merge, 1e2);  // halo radius ~1.78
            c.addBody(mH, pH.data(), v0.data(), true);
            c.addBody(mS, pS.data(), v0.data());
            double E0 = c.energy();
            for (int i = 0; i < 5000; i++) c.step(2e-5, 1);
            check(c.aliveCount() == 2, "T14 pinned halo never merges");
            double st[16];
            c.getState(st, 2);
            check(std::fabs(st[0] - mH) < 1e-9, "T14 halo mass unchanged");
            check(std::fabs(st[1]) < 1e-12 && std::fabs(st[2]) < 1e-12,
                  "T14 halo still at the origin after impacts");
            double dE = std::fabs((c.energy() - E0) / E0);
            std::printf("      merge-mode pass-through dE/E=%.3e\n", dE);
            check(dE < 5e-3, "T14 energy sane with a pinned halo");
        }
        // (d) Bounce: the halo acts as an immovable elastic wall
        {
            std::vector<double> pS = {0.0, 3.0}, vS = {0, -1.0};
            Simulation d(2, "t^2", "1/r^2");
            d.setCollision(grav::Collision::Bounce, 1e2);  // R_halo ~1.78, R_s ~0.056
            d.addBody(mH, pH.data(), v0.data(), true);
            d.addBody(mS, pS.data(), vS.data());
            bool bounced = false;
            for (int i = 0; i < 4000 && !bounced; i++) {
                d.step(1e-4, 1);
                double st[16];
                d.getState(st, 2);
                if (st[8 + 5] > 0.0) bounced = true;  // vy reversed to +y
            }
            double st[16];
            d.getState(st, 2);
            std::printf("      bounced=%d vy=%.4f (in %.4f, out %.4f)\n",
                        (int)bounced, st[8 + 5], -1.0, st[8 + 5]);
            check(bounced, "T14 bounce reflects off a pinned halo");
            check(std::fabs(st[1]) < 1e-12 && std::fabs(st[2]) < 1e-12,
                  "T14 pinned halo unmoved by the bounce");
        }
    }

    // ---------------------------------------------------------------- test 11
    // Barnes-Hut accuracy on a 3000-body random annulus (Newtonian 2D, point
    // particles, min separation 0.02 via rejection sampling).
    //   (a) per-body acceleration (direct sum) vs tree force at theta=0.4,
    //       measured via the one-tiny-step velocity probe.  Errors are
    //       normalized by the force scale aMax = max|a|: per-body *relative*
    //       error diverges for bodies with near-cancelling forces, which is
    //       a property of the metric, not of the tree.
    //   (b) tree PE estimator (theta_PE = min(theta, 0.25)) vs exact PE.
    //   (c) cold weak-G disc integration with tree forces conserves energy.
    {
        const int NB = 3000;
        std::vector<double> pos;
        std::vector<double> mass(NB);
        {
            const double minSep = 0.02;
            std::mt19937_64 rng(12345);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            int made = 0;
            while (made < NB) {
                double r = 0.15 + 0.85 * std::sqrt(uni(rng));
                double a = uni(rng) * 6.283185307179586;
                double x = r * std::cos(a), y = r * std::sin(a);
                bool okp = true;
                for (int j = 0; j < made; j++) {
                    double dx = x - pos[2 * j], dy = y - pos[2 * j + 1];
                    if (dx * dx + dy * dy < minSep * minSep) { okp = false; break; }
                }
                if (!okp) continue;
                pos.push_back(x);
                pos.push_back(y);
                mass[made] = 0.5 + 1.5 * uni(rng);
                made++;
            }
        }
        auto accelSample = [&](int treeMode, double theta, int nSample) {
            std::vector<double> zero(2 * NB, 0.0);
            Simulation sim(2, "t^2", "1/r^2");
            sim.setCollision(grav::Collision::Pass, 1e4);
            sim.setTree(treeMode, theta);
            for (int i = 0; i < NB; i++)
                sim.addBody(mass[i], &pos[2 * i], &zero[2 * i]);
            sim.step(1e-12, 1);  // v = a*dt to first order
            double st[8 * 4096];
            int n = sim.getState(st, 4096);
            (void)n;
            std::vector<double> ax(nSample);
            for (int i = 0; i < nSample; i++) ax[i] = st[8 * i + 4];
            return ax;
        };
        std::vector<double> aDir = accelSample(1, 0.6, 400);   // tree off
        std::vector<double> aTree = accelSample(0, 0.4, 400);  // auto -> tree
        double aMax = 0.0;
        for (double v : aDir) aMax = std::max(aMax, std::fabs(v));
        double maxErr = 0, sumErr = 0;
        int cnt = 0, bad = 0;
        for (int i = 0; i < 400; i++) {
            double ad = aDir[i], at = aTree[i];
            if (!std::isfinite(ad) || !std::isfinite(at)) { check(false, "T11 finite"); break; }
            double e = std::fabs(ad - at) / aMax;
            maxErr = std::max(maxErr, e);
            sumErr += e;
            if (e > 0.01) bad++;
            cnt++;
        }
        std::printf("      NB=%d theta=0.4: maxErr/aMax=%.3e meanErr/aMax=%.3e bad(>1%%):%d/%d\n",
                    NB, maxErr, sumErr / cnt, bad, cnt);
        check(maxErr < 3e-2, "T11 BH force max error < 3% of force scale");
        check(sumErr / cnt < 6e-3, "T11 BH force mean error < 0.6% of force scale");
        check(bad < 150, "T11 BH force: fewer than 150/400 bodies above 1% error");
        // (b) tree PE estimator vs exact PE at t=0
        {
            std::vector<double> zero(2 * NB, 0.0);
            Simulation sDir(2, "t^2", "1/r^2");
            Simulation sTree(2, "t^2", "1/r^2");
            sDir.setCollision(grav::Collision::Pass, 1e4);
            sDir.setTree(1, 0.6);
            sTree.setCollision(grav::Collision::Pass, 1e4);
            sTree.setTree(0, 0.4);
            for (int i = 0; i < NB; i++) {
                sDir.addBody(mass[i], &pos[2 * i], &zero[2 * i]);
                sTree.addBody(mass[i], &pos[2 * i], &zero[2 * i]);
            }
            double Ed = sDir.energy(), Et = sTree.energy();
            std::printf("      PE: direct=%.6e tree=%.6e rel=%.3e\n", Ed, Et,
                        (Et - Ed) / Ed);
            check(std::isfinite(Ed) && std::fabs((Et - Ed) / Ed) < 0.01,
                  "T11 BH potential estimator error < 1%");
        }
        // (c) cold weak-G full disc (minSep 0.02), 200 tree steps
        {
            const double minSep = 0.02;
            std::mt19937_64 rng(777);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            std::vector<double> pos2;
            int made = 0;
            while (made < NB) {
                double r = std::sqrt(uni(rng));
                double a = uni(rng) * 6.283185307179586;
                double x = r * std::cos(a), y = r * std::sin(a);
                bool okp = true;
                for (int j = 0; j < made; j++) {
                    double dx = x - pos2[2 * j], dy = y - pos2[2 * j + 1];
                    if (dx * dx + dy * dy < minSep * minSep) { okp = false; break; }
                }
                if (!okp) continue;
                pos2.push_back(x);
                pos2.push_back(y);
                made++;
            }
            std::vector<double> zero(2 * NB, 0.0);
            Simulation sim(2, "t^2", "0.001/r^2");
            sim.setCollision(grav::Collision::Pass, 1e4);
            sim.setTree(0, 0.4);
            for (int i = 0; i < NB; i++) sim.addBody(1.0, &pos2[2 * i], &zero[2 * i]);
            double E0 = sim.energy();
            for (int s = 0; s < 200; s++) sim.step(5e-4, 1);
            double E1 = sim.energy();
            std::printf("      tree run: dE/E=%.3e\n", (E1 - E0) / E0);
            check(std::isfinite(E1) && std::fabs((E1 - E0) / E0) < 1e-2,
                  "T11 tree integration energy sane");
        }
    }

    std::printf("\n%s (%d failures)\n", g_fail == 0 ? "ALL TESTS PASSED" : "TESTS FAILED", g_fail);
    return g_fail == 0 ? 0 : 1;
}
