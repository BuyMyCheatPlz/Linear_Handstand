#include "lqr_inverted_pendulum.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "═══════════════════════════════════════════════\n";
    std::cout << "  一阶线性倒立摆 LQR 控制器\n";
    std::cout << "═══════════════════════════════════════════════\n\n";

    // ── 1. 定义物理参数 ──────────────────────────
    PendulumParams params;
    params.M = 1.0;     // 小车质量 1 kg
    params.m = 0.1;     // 摆杆质量 0.1 kg
    params.l = 0.5;     // 摆杆长度 0.5 m
    params.g = 9.81;    // 重力加速度
    params.d = 0.05;    // 阻尼系数

    std::cout << "── 物理参数 ──\n";
    std::cout << "  小车质量 M    = " << params.M << " kg\n";
    std::cout << "  摆杆质量 m    = " << params.m << " kg\n";
    std::cout << "  摆杆长度 l    = " << params.l << " m\n";
    std::cout << "  重力加速度 g  = " << params.g << " m/s²\n";
    std::cout << "  阻尼系数 d    = " << params.d << "\n\n";

    // ── 2. 构造状态空间矩阵 ──────────────────────
    Eigen::MatrixXd A, B;
    buildPendulumSystem(params, A, B);

    std::cout << "── 开环系统矩阵 ──\n";
    std::cout << "A = \n" << A << "\n\n";
    std::cout << "B = \n" << B << "\n\n";

    // 开环特征值（应有正实部，不稳定）
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> openEig(A);
    std::cout << "开环特征值: ";
    for (int i = 0; i < openEig.eigenvalues().size(); ++i) {
        std::cout << openEig.eigenvalues()(i).real()
                  << (openEig.eigenvalues()(i).imag() >= 0 ? " + " : " - ")
                  << std::abs(openEig.eigenvalues()(i).imag()) << "i  ";
    }
    std::cout << "\n\n";

    // ── 3. 定义 LQR 权重矩阵 ─────────────────────
    // Q = diag(q_theta, q_thetadot): 惩罚状态偏差
    //    增大 q_theta  → 更快回到竖直位置
    //    增大 q_thetadot → 抑制摆动速度
    // R: 惩罚控制输入（力），增大 R → 更保守的控制

    Eigen::MatrixXd Q(2, 2);
    Q << 100.0, 0.0,       // 角度偏差权重大 → 快速稳定
           0.0, 1.0;       // 角速度权重较小

    Eigen::MatrixXd R(1, 1);
    R << 0.01;             // 控制代价小 → 允许较大的力

    std::cout << "── LQR 权重矩阵 ──\n";
    std::cout << "Q = \n" << Q << "\n\n";
    std::cout << "R = \n" << R << "\n\n";

    // ── 4. 求解 LQR ─────────────────────────────
    LqrGain gain = solveLQR(A, B, Q, R);

    if (!gain.converged) {
        std::cerr << "LQR 求解失败！\n";
        return 1;
    }

    std::cout << "── LQR 结果 ──\n";
    std::cout << "Riccati 解 P = \n" << gain.P << "\n\n";
    std::cout << "最优反馈增益 K = [ " << gain.K(0, 0) << "  "
              << gain.K(0, 1) << " ]\n\n";

    // 闭环特征值
    Eigen::MatrixXd Acl = A - B * gain.K;
    Eigen::ComplexEigenSolver<Eigen::MatrixXd> clEig(Acl);
    std::cout << "闭环特征值: ";
    for (int i = 0; i < clEig.eigenvalues().size(); ++i) {
        std::cout << clEig.eigenvalues()(i).real()
                  << (clEig.eigenvalues()(i).imag() >= 0 ? " + " : " - ")
                  << std::abs(clEig.eigenvalues()(i).imag()) << "i  ";
    }
    std::cout << "\n(特征值实部为负 → 系统稳定)\n\n";

    // ── 5. 仿真 ─────────────────────────────────
    double dt = 0.001;      // 1 ms 步长
    int steps = 5000;       // 仿真 5 秒

    // 初始状态: θ = 15°, θ̇ = 0
    Eigen::VectorXd x0(2);
    x0 << 15.0 * M_PI / 180.0,  // 15° → rad
           0.0;

    std::cout << "── 仿真设置 ──\n";
    std::cout << "  初始角度 θ₀ = 15°\n";
    std::cout << "  初始角速度 θ̇₀ = 0\n";
    std::cout << "  步长 dt = " << dt << " s\n";
    std::cout << "  步数 = " << steps << " (共 " << dt * steps << " s)\n\n";

    Eigen::MatrixXd traj = simulate(x0, params, gain.K, dt, steps);

    // 构造时间向量
    Eigen::VectorXd time(steps + 1);
    for (int i = 0; i <= steps; ++i) time(i) = i * dt;

    // ── 6. 输出关键指标 ─────────────────────────
    double theta_steady = traj(steps, 0);
    double thetadot_steady = traj(steps, 1);

    // 计算调节时间（角度进入 ±2% 初始偏差范围内）
    double target_band = 0.02 * std::abs(x0(0));
    double settling_time = 0.0;
    for (int i = 0; i <= steps; ++i) {
        bool settled = true;
        for (int j = i; j <= steps; ++j) {
            if (std::abs(traj(j, 0)) > target_band) {
                settled = false;
                break;
            }
        }
        if (settled) {
            settling_time = i * dt;
            break;
        }
    }

    // 最大控制力
    double max_force = 0.0;
    for (int i = 0; i < steps; ++i) {
        double u = -(gain.K * traj.row(i).transpose())(0);
        if (std::abs(u) > max_force) max_force = std::abs(u);
    }

    std::cout << "── 仿真结果 ──\n";
    std::cout << "  稳态角度 θ(∞)  = " << theta_steady * 180.0 / M_PI << "°\n";
    std::cout << "  稳态角速度      = " << thetadot_steady << " rad/s\n";
    std::cout << "  调节时间 (±2%)  = " << settling_time << " s\n";
    std::cout << "  最大控制力      = " << max_force << " N\n";

    // ── 7. 打印部分轨迹 ─────────────────────────
    std::cout << "\n── 时间历程 (前 20 步 + 末尾 5 步) ──\n";
    std::cout << "  t(s)        θ(°)        θ̇(rad/s)    u(N)\n";
    std::cout << "  ──────────────────────────────────────────\n";
    // 前 10 步 (每 5 步打一次)
    for (int i = 0; i <= 50; i += 5) {
        double u = -(gain.K * traj.row(i).transpose())(0);
        std::cout << "  " << std::setw(10) << time(i)
                  << "  " << std::setw(10) << traj(i, 0) * 180.0 / M_PI
                  << "  " << std::setw(10) << traj(i, 1)
                  << "  " << std::setw(10) << u << "\n";
    }
    std::cout << "  ...\n";
    for (int i = steps - 5; i <= steps; i++) {
        double u = -(gain.K * traj.row(i).transpose())(0);
        std::cout << "  " << std::setw(10) << time(i)
                  << "  " << std::setw(10) << traj(i, 0) * 180.0 / M_PI
                  << "  " << std::setw(10) << traj(i, 1)
                  << "  " << std::setw(10) << u << "\n";
    }

    // ── 8. 导出数据给 MATLAB ────────────────────
    exportCSV(traj, time, "pendulum_sim.csv");

    // 同时导出 LQR 增益供 MATLAB 对照（纯数值，无表头，readmatrix 直接可用）
    {
        std::ofstream fk("lqr_gain_cpp.csv");
        fk << std::fixed << std::setprecision(12);
        fk << gain.K(0, 0) << "," << gain.K(0, 1) << "\n";
        fk.close();
        std::cout << "反馈增益已保存至: lqr_gain_cpp.csv\n";
    }

    {
        std::ofstream fp("riccati_P_cpp.csv");
        fp << std::fixed << std::setprecision(12);
        fp << gain.P(0, 0) << "," << gain.P(0, 1) << "\n";
        fp << gain.P(1, 0) << "," << gain.P(1, 1) << "\n";
        fp.close();
        std::cout << "Riccati 矩阵已保存至: riccati_P_cpp.csv\n";
    }

    std::cout << "\n✅ LQR 求解与仿真完成！\n";
    return 0;
}
