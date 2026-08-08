#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

/**
 * 一阶线性倒立摆 LQR 控制器
 *
 * 物理模型（小车-摆杆系统，线性化小角度近似 sinθ≈θ, cosθ≈1）:
 *
 *   状态变量:  x = [θ, θ̇]ᵀ         (θ: 摆杆偏角 [rad], θ̇: 角速度 [rad/s])
 *   控制输入:  u = F               (作用在小车上的力 [N])
 *
 *   连续时间状态空间:  dx/dt = A x + B u
 *
 *     A = [[       0,        1 ],
 *          [(M+m)g/(M*l),  -d ]]
 *
 *     B = [[      0    ],
 *          [ -1/(M*l)  ]]
 *
 *   参数:
 *     M — 小车质量 [kg]
 *     m — 摆杆质量 [kg]
 *     l — 摆杆长度（质心到转轴）[m]
 *     g — 重力加速度 [m/s²]
 *     d — 阻尼系数 [(N·m·s)/rad]
 *
 * LQR 代价函数:
 *     J = ∫₀∞ (xᵀQx + uᵀRu) dt
 *
 *   通过求解连续时间代数 Riccati 方程 (CARE):
 *     AᵀP + PA − PBR⁻¹BᵀP + Q = 0
 *
 *   得到最优状态反馈增益:
 *     K = R⁻¹BᵀP
 *     u = −Kx
 */

struct PendulumParams {
    double M = 1.0;    // 小车质量 [kg]
    double m = 0.1;    // 摆杆质量 [kg]
    double l = 0.5;    // 摆杆长度 [m]
    double g = 9.81;   // 重力加速度 [m/s²]
    double d = 0.05;   // 阻尼系数
};

struct LqrGain {
    Eigen::MatrixXd K;   // 反馈增益矩阵 (1×n)
    Eigen::MatrixXd P;   // Riccati 解 (n×n)
    bool converged;      // CARE 是否收敛
};

// 构造一阶倒立摆的状态矩阵 A (2×2) 和输入矩阵 B (2×1)
void buildPendulumSystem(const PendulumParams& params,
                         Eigen::MatrixXd& A,
                         Eigen::MatrixXd& B);

// 求解连续时间代数 Riccati 方程 (CARE) 使用 Schur 分解法
//   输入: A (n×n), B (n×m), Q (n×n 半正定), R (m×m 正定)
//   输出: LqrGain 包含增益 K 与 Riccati 解 P
LqrGain solveLQR(const Eigen::MatrixXd& A,
                 const Eigen::MatrixXd& B,
                 const Eigen::MatrixXd& Q,
                 const Eigen::MatrixXd& R);

// 使用四阶 Runge-Kutta 仿真闭环系统
//   x0:     初始状态
//   params: 摆杆参数
//   K:      反馈增益
//   dt:     仿真步长 [s]
//   steps:  仿真步数
//   返回:   每个时刻的状态 (steps×n), 每行一个时刻
Eigen::MatrixXd simulate(const Eigen::VectorXd& x0,
                         const PendulumParams& params,
                         const Eigen::MatrixXd& K,
                         double dt,
                         int steps);

// 将仿真结果输出到 CSV 文件，方便 MATLAB 对照
void exportCSV(const Eigen::MatrixXd& trajectory,
               const Eigen::VectorXd& time,
               const std::string& filename);
