#include "lqr_inverted_pendulum.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <numeric>

// ──────────────────────────────────────────────
//  构造一阶倒立摆状态空间矩阵
// ──────────────────────────────────────────────
void buildPendulumSystem(const PendulumParams& params,
                         Eigen::MatrixXd& A,
                         Eigen::MatrixXd& B) {
    const double M = params.M;
    const double m = params.m;
    const double l = params.l;
    const double g = params.g;
    const double d = params.d;

    // 总质量
    double Mt = M + m;

    // A = [[        0,     1 ],
    //      [ Mt*g/(M*l),  -d ]]
    A.resize(2, 2);
    A << 0.0,          1.0,
         Mt * g / (M * l), -d;

    // B = [[    0    ],
    //      [ -1/(M*l) ]]
    B.resize(2, 1);
    B << 0.0,
         -1.0 / (M * l);
}

// ──────────────────────────────────────────────
//  CARE:  AᵀP + PA − PBR⁻¹BᵀP + Q = 0
//  使用 Schur 分解法 (Hamiltonian 矩阵方法)
// ──────────────────────────────────────────────
LqrGain solveLQR(const Eigen::MatrixXd& A,
                 const Eigen::MatrixXd& B,
                 const Eigen::MatrixXd& Q,
                 const Eigen::MatrixXd& R) {
    const int n = A.rows();               // 状态维度
    const int m = B.cols();               // 输入维度

    LqrGain result;
    result.K.resize(m, n);
    result.P.resize(n, n);
    result.converged = false;

    // --- 输入检查 ---
    // Q 半正定, R 正定
    Eigen::LLT<Eigen::MatrixXd> lltR(R);
    if (lltR.info() != Eigen::Success) {
        std::cerr << "[LQR] R 不是正定矩阵！\n";
        return result;
    }

    // R^{-1}
    Eigen::MatrixXd Rinv = R.inverse();

    // ── 构造 Hamiltonian 矩阵 H (2n × 2n) ──
    //   H = [  A,   -B R^{-1} Bᵀ ]
    //       [ -Q,       -Aᵀ      ]
    //
    // 性质: H 的特征值关于虚轴对称。
    //   H 有 n 个稳定特征值 (Re < 0) 和 n 个不稳定特征值 (Re > 0).
    //   稳定特征值对应的不变子空间 span 了 CARE 的解。
    //
    // 注意: 对表达式模板显式 .eval() 避免 -O2 下的悬垂引用
    Eigen::MatrixXd S = B * Rinv * B.transpose();  // S = B R⁻¹ Bᵀ

    Eigen::MatrixXd H(2 * n, 2 * n);
    H << A,                -S,
         -Q,               -A.transpose();

    // ── 特征分解 ──
    // H * V = V * Λ,  其中 V = [V11  V12;  V21  V22]
    // 取稳定特征向量（Re(λ) < 0），这些向量组成 [V11; V21]
    // CARE 的解:  P = V21 * V11^{-1}
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> eigSolver(H);

    if (eigSolver.info() != Eigen::Success) {
        std::cerr << "[LQR] Hamiltonian 特征分解失败！\n";
        return result;
    }

    Eigen::VectorXcd eigvals  = eigSolver.eigenvalues();
    Eigen::MatrixXcd eigvecs  = eigSolver.eigenvectors();

    // 按特征值实部升序排列 (最 stable 的在前)
    std::vector<int> order(2 * n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return eigvals(a).real() < eigvals(b).real();
    });

    // 取前 n 个特征向量（属于稳定子空间）
    Eigen::MatrixXcd U11(n, n);
    Eigen::MatrixXcd U21(n, n);
 
    for (int k = 0; k < n; ++k) {
        int col = order[k];
        for (int i = 0; i < n; ++i) {
            U11(i, k) = eigvecs(i, col);       // 上半块 n×n
            U21(i, k) = eigvecs(n + i, col);   // 下半块 n×n
        }
    }

    // P = U21 * U11^{-1}   (理论上 P 为实对称矩阵)
    Eigen::MatrixXcd P_complex = U21 * U11.inverse();

    result.P = P_complex.real();                     // 取实部
    result.P = 0.5 * (result.P + result.P.transpose()); // 对称化

    // --- 验证 Riccati 方程残差 ---
    {
        Eigen::MatrixXd residual = A.transpose() * result.P
                                 + result.P * A
                                 - result.P * B * Rinv * B.transpose() * result.P
                                 + Q;
        double resNorm = residual.norm();
        if (resNorm > 1e-6) {
            std::cout << "[LQR] CARE 残差范数 = " << resNorm
                      << " (偏大，解可能不精确)\n";
        }
    }

    // --- 反馈增益 K = R^{-1} Bᵀ P ---
    result.K = Rinv * B.transpose() * result.P;
    result.converged = true;

    return result;
}

// ──────────────────────────────────────────────
//  Runge-Kutta 4 阶仿真
// ──────────────────────────────────────────────
Eigen::MatrixXd simulate(const Eigen::VectorXd& x0,
                         const PendulumParams& params,
                         const Eigen::MatrixXd& K,
                         double dt,
                         int steps) {
    int n = x0.size();
    Eigen::MatrixXd trajectory(steps + 1, n);
    Eigen::VectorXd x = x0;
    trajectory.row(0) = x.transpose();

    // 构造开环矩阵
    Eigen::MatrixXd A, B;
    buildPendulumSystem(params, A, B);

    // 闭环动态: dx/dt = (A - B K) x
    Eigen::MatrixXd Acl = A - B * K;

    for (int k = 0; k < steps; ++k) {
        // RK4
        Eigen::VectorXd k1 = Acl * x;
        Eigen::VectorXd k2 = Acl * (x + 0.5 * dt * k1);
        Eigen::VectorXd k3 = Acl * (x + 0.5 * dt * k2);
        Eigen::VectorXd k4 = Acl * (x + dt * k3);

        x = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        trajectory.row(k + 1) = x.transpose();
    }

    return trajectory;
}

// ──────────────────────────────────────────────
//  导出 CSV 文件
// ──────────────────────────────────────────────
void exportCSV(const Eigen::MatrixXd& trajectory,
               const Eigen::VectorXd& time,
               const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << "\n";
        return;
    }

    file << std::fixed << std::setprecision(8);
    // 纯数值 CSV (无表头), MATLAB readmatrix 直接可用

    int rows = static_cast<int>(trajectory.rows());
    for (int i = 0; i < rows; ++i) {
        file << time(i) << ","
             << trajectory(i, 0) << ","
             << trajectory(i, 1) << "\n";
    }

    file.close();
    std::cout << "仿真数据已保存至: " << filename << "\n";
}
