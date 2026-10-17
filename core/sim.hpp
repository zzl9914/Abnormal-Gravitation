#pragma once
// 自定义度量 + 自定义力律的多体引力模拟。
//
// 度量：   d = F_m( sum_k f(|dx_k|) ),  F_m = f^-1，f 在 [0,inf) 上严格单增、f(0)=0
// 力律（用户表达式，变量 r = 度量距离 d）：
//     a_i = m_j * F(d_ij) * (x_j - x_i) / |x_j - x_i|
//   F > 0 表示吸引（力指向施力天体），F < 0 表示排斥。力恒为中心力（沿连线），
//   因此**总角动量对任意度量都被速度 Verlet 精确守恒**。
//
// 势（派生，用于能量监视 / 核内正则 / 树的 PE 估计）：V(r) = ∫F dr，即 V' = F。
//   F 是幂和 → 解析精确；否则一次性建数值表（对数网格 4 点 Gauss-Legendre +
//   三次 Hermite 插值）。能量只在 f = t^2（各向同性）时严格守恒。
//
// 天体半径由物质密度导出：R = ( m·Γ(d/2+1) / (ρ·π^(d/2)) )^(1/d)
//
// 核内正则：d < R_ij 时力幅换成 F(R_ij)·d/R_ij，力在核内线性衰减到 0，
//   r→0 的奇点永不被求值（只需求 V(R)、F(R)）。
//
// 碰撞（d <= R_ij）：Pass 穿过（默认）/ Merge 完全非弹性并合 / Bounce 完全弹性脉冲。
//
// Barnes-Hut：auto 模式下存活体数 >= kTreeThreshold 启用；判据用度量距离
//   cellMetricRadius < theta * d。碰撞检测用「点到盒的度量距离下界」精确剪枝。
//
// ⚠ freestanding 约束（为 wasm 构建）：无 std::string / std::vector / new。
//   天体、树节点都是**固定容量**数组；默认构造必须平凡（静态存储靠 .bss 清零，
//   真正的默认值在 init() 里设置）。
//   ⚠ 树节点缓冲是**文件作用域共享静态数组**：同一时刻只有一个树在活动，
//     computeAccel / energy / collectContacts 每次都会重建。

#include <stdint.h>

#include "expr.hpp"

namespace grav {

enum class Collision : int { Pass = 0, Merge = 1, Bounce = 2 };

// Body flag bits（状态流第 8 个分量）。
// fixed（暗物质晕）永不移动、永不并合、被积分器跳过，但对所有天体施加引力；
// Bounce 模式下表现为不可推动的弹性壁。
enum BodyFlag : unsigned { kFlagFixed = 1u };

// 固定容量
enum {
    kMaxBodies = 4096,
    kMaxNodes = 16384,
    kMaxExtra = 2048,     // 深度上限处重合天体的溢出槽位
    kMaxContacts = 4096,  // 单次收集的接触对上限
    kTabN = 4096,         // 数值势表的点数
    kMaxPowTermsSim = 16
};

struct Body {
    double pos[8];
    double vel[8];
    double acc[8];
    double mass;
    double radius;
    bool alive;
    bool fixed;  // 钉住（晕）：不移动、不并合
};

// 四叉 / 八叉树节点。state: 0=空 1=叶 2=内部 3=bucket（深度封顶）
struct TreeNode {
    double mass;
    double com[8];   // 质心
    double ctr[8];   // 几何中心
    double half[8];  // 半边长
    double mrad;     // 度量半径界：max d(com, 角点)
    double maxRad;   // 子树内最大天体半径
    int child[8];
    int body;
    int parent;
    int extra;  // 溢出链表头（-1 表示无），替代原来的 vector<int>
    uint8_t state;
};

class Simulation {
public:
    // 平凡默认构造：配合 init() 使用，便于放进静态槽位数组（wasm 无 new）。
    Simulation() {}
    Simulation(int dims, const char* fExpr, const char* forceExpr) {
        init(dims, fExpr, forceExpr);
    }

    // 全量（重新）初始化。等价于原来的构造函数。
    bool init(int dims, const char* fExpr, const char* forceExpr);

    bool ok() const { return ok_; }
    const char* error() const { return err_; }

    void setCollision(Collision mode, double density);
    // 力律耦合常数：力与派生势同步乘 G（G > 0，默认 1）。
    // a_i = G·m_j·F(d)·u。前端用它标定时标：保持 G·F(R/2) 不变。
    void setG(double G);
    double g() const { return G_; }
    // 原始力律值 F(r)（不含 G）
    double forceAt(double r) const { return F_.value(r); }
    // treeMode: 0 = auto（alive >= kTreeThreshold 时用树），1 = 始终直接求和
    void setTree(int treeMode, double theta);

    // fixed 天体被钉住：施力但不动、不并合，反弹时按无穷大质量反射对方
    int addBody(double mass, const double* pos, const double* vel, bool fixed = false);
    void clearBodies();  // 同时把 t 归零

    bool step(double dt, int substeps);  // substeps 步速度 Verlet
    double time() const { return t_; }
    double energy() const;               // KE + PE
    int aliveCount() const;
    int dims() const { return dim_; }
    int bodyCount() const { return nbodies_; }

    // 打包每体 8 个 double：[mass, x, y, z, vx, vy, vz, flags]（2D 补零）
    int getState(double* out, int maxBodies) const;

    // 派生势 V(r) 在 [rmin, rmax] 上的 n 个对数采样点
    void samplePotential(int n, double rmin, double rmax, double* outR, double* outV) const;
    // 度量单位球（前两轴切片）：n 个点 (x, y)；不可达半径输出 (NAN, NAN)
    int unitBall2D(int n, double* outXY) const;

    // 测试专用
    double metricDistance(const double* delta) const;
    bool pairKernel(const double* delta, double& dOut, double* grad) const;
    double potentialAt(double r) const;
    double potentialEnergy() const;
    double radiusFromMass(double m) const;

private:
    void computeAccel();
    void handleCollisions();
    double invertF(double S) const;
    double pairV(double d, double R) const;
    void buildPotentialTable();

    // ---- Barnes-Hut ----
    static const int kMaxDepth = 48;
    static const int kTreeThreshold = 1500;
    bool treeActive() const;
    int newNode(int parent, const double* ctr, const double* half) const;
    void insertBody(int ni, int b, int depth) const;
    void buildTree() const;
    void treeForce();
    double treePotential() const;
    void applyOneSided(int i, const double* posj, double mj, double radj);
    double pairUOneSided(int i, const double* posj, double mj, double radj) const;
    void collectContacts() const;
    double minMetricDistToBox(const double* xi, const TreeNode& nd) const;
    void resolveContact(int i, int j);

    // ⚠ 所有成员都**没有**默认初始化器：静态实例靠 .bss 清零，
    // 真正的默认值（density=100, G=1, theta=0.6...）在 init() 中设置。
    int dim_;
    Expr f_, F_;  // 度量生成元 f(t)；力律 F(r)
    bool ok_;
    char err_[kMaxErr];
    bool fIsPower_;
    double fExp_;
    // 派生势：V = sum a_i r^{q_i} + vLn_·ln r（解析）或数值表
    double vPowP_[kMaxPowTermsSim];
    double vPowA_[kMaxPowTermsSim];
    int nPow_;
    double vLn_;
    bool vAnalytic_;
    double vTabR_[kTabN];
    double vTabV_[kTabN];
    double vTabF_[kTabN];
    Body bodies_[kMaxBodies];
    int nbodies_;
    double t_;
    Collision collisionMode_;
    double density_;
    double G_;
    int treeMode_;
    double theta_;
    bool accValid_;
};

}  // namespace grav
