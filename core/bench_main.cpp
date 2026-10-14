// Performance benchmark: direct summation vs Barnes-Hut tree.
//   bin/gs_bench.exe [--3d] [--direct-max N] [--steps N]
// Prints ms/step and speedup tables for a uniform random disc/ball.
#include "sim.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using grav::Simulation;
using grav::Collision;

struct Result {
    double msPerStep;
    int steps;
};

static Result bench(int dims, int n, int treeMode, double theta, int wantSteps) {
    std::mt19937_64 rng(4242);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    Simulation sim(dims, "t^2", "1/r^2");
    if (!sim.ok()) {
        std::printf("  construct failed: %s\n", sim.error().c_str());
        return {0.0, 0};
    }
    sim.setCollision(Collision::Pass, 100.0);
    sim.setTree(treeMode, theta);
    std::vector<double> p(dims), v(dims, 0.0);
    for (int i = 0; i < n; i++) {
        double r = std::sqrt(uni(rng)) * (0.15 + 0.85 * uni(rng));
        if (dims == 2) {
            double a = uni(rng) * 6.283185307179586;
            p[0] = r * std::cos(a);
            p[1] = r * std::sin(a);
        } else {
            double z = 2.0 * uni(rng) - 1.0;
            double a = uni(rng) * 6.283185307179586;
            double rr = std::sqrt(std::max(0.0, 1.0 - z * z));
            p[0] = r * rr * std::cos(a);
            p[1] = r * rr * std::sin(a);
            p[2] = r * z;
        }
        sim.addBody(1.0, p.data(), v.data());
    }
    // Warm-up: one step forces tree build + one traversal + PE cache.
    sim.step(1e-4, 1);
    // Adapt step count so one measurement is ~1 s but never fewer than 2.
    int steps = wantSteps;
    auto t0 = std::chrono::steady_clock::now();
    for (int s = 0; s < steps; s++) sim.step(1e-4, 1);
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (sec < 0.2 && steps < 4096) {  // too fast to be reliable: rerun longer
        int more = (int)(0.8 / std::max(1e-9, sec / steps));
        if (more < steps) more = steps;
        if (more > 8192) more = 8192;
        t0 = std::chrono::steady_clock::now();
        for (int s = 0; s < more; s++) sim.step(1e-4, 1);
        sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        steps = more;
    }
    return {1000.0 * sec / steps, steps};
}

int main(int argc, char** argv) {
    int dims = 2;
    int directMax = 4000;
    int wantSteps = 4;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--3d") dims = 3;
        else if (a == "--direct-max" && i + 1 < argc) directMax = std::atoi(argv[++i]);
        else if (a == "--steps" && i + 1 < argc) wantSteps = std::atoi(argv[++i]);
    }
    std::printf("Gravity benchmark  dims=%d  f=t^2  F=1/r^2  collision=pass\n", dims);
    std::printf("%8s %10s %12s %14s\n", "N", "mode", "ms/step", "steps measured");
    const int Ns[] = {500, 1000, 2000, 4000, 8000, 16000, 32000};
    for (int n : Ns) {
        if (n <= directMax) {
            Result r = bench(dims, n, 1, 0.6, wantSteps);
            if (r.steps == 0) continue;
            std::printf("%8d %10s %12.2f %14d\n", n, "direct", r.msPerStep, r.steps);
        }
        for (double th : {0.6, 0.4}) {
            Result r = bench(dims, n, 0, th, wantSteps);
            if (r.steps == 0) continue;
            char mode[16];
            std::snprintf(mode, sizeof(mode), "tree %.2f", th);
            std::printf("%8d %10s %12.2f %14d\n", n, mode, r.msPerStep, r.steps);
        }
        std::fflush(stdout);
    }
    return 0;
}
