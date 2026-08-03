%% verify_lqr.m
% 一阶线性倒立摆 LQR 控制器 — MATLAB 验证脚本
%
% 功能:
%   1. 独立实现 LQR 求解与仿真（使用 MATLAB 内置 care / lqr 函数）
%   2. 读取 C++ 输出结果，进行数值比对验证
%   3. 绘制状态轨迹、相平面图、控制输入曲线
%
% 运行方式: 在 MATLAB 中直接运行此脚本

clear; clc; close all;

fprintf('═══════════════════════════════════════════════════\n');
fprintf('  一阶线性倒立摆 LQR — MATLAB 验证\n');
fprintf('═══════════════════════════════════════════════════\n\n');

%% ── 1. 物理参数 ────────────────────────────────────
M = 1.0;      % 小车质量 [kg]
m = 0.1;      % 摆杆质量 [kg]
l = 0.5;      % 摆杆长度 [m]
g = 9.81;     % 重力加速度 [m/s^2]
d = 0.05;     % 阻尼系数

Mt = M + m;   % 总质量

fprintf('── 物理参数 ──\n');
fprintf('  M = %.2f kg,  m = %.2f kg,  l = %.2f m,  g = %.2f,  d = %.3f\n\n', M, m, l, g, d);

%% ── 2. 状态空间矩阵 ────────────────────────────────
% 状态 x = [θ; θ̇]
% dx/dt = A x + B u

A = [0,         1;
     Mt*g/(M*l), -d];

B = [0;
     -1/(M*l)];

fprintf('── 开环系统矩阵 ──\n');
fprintf('A = \n'); disp(A);
fprintf('B = \n'); disp(B);

fprintf('开环特征值: ');
disp(eig(A)');

%% ── 3. LQR 设计 ────────────────────────────────────
% 权重矩阵（与 C++ 一致）
Q = diag([100, 1]);   % 角度偏差权重 >> 角速度权重
R = 0.01;             % 控制代价

fprintf('\n── LQR 权重 ──\n');
fprintf('Q = \n'); disp(Q);
fprintf('R = \n'); disp(R);

% 方法 1: 使用 MATLAB 内置 lqr 函数（推荐）
[K_lqr, P_lqr, cl_eig] = lqr(A, B, Q, R);

% 方法 2: 使用 care 求解器手动计算（验证用）
[P_care, ~, L_care] = care(A, B, Q, R);
K_care = R \ (B' * P_care);

fprintf('\n── LQR 结果 (MATLAB lqr) ──\n');
fprintf('Riccati 解 P = \n'); disp(P_lqr);
fprintf('反馈增益 K = [ %.8f  %.8f ]\n', K_lqr(1), K_lqr(2));
fprintf('闭环特征值: '); disp(cl_eig');

% 验证 care 和 lqr 一致性
fprintf('‖K_lqr - K_care‖ = %.2e\n', norm(K_lqr - K_care));
fprintf('‖P_lqr - P_care‖ = %.2e\n', norm(P_lqr - P_care));

%% ── 4. 仿真 ────────────────────────────────────────
dt = 0.001;           % 1 ms 步长
T_final = 5.0;        % 仿真 5 秒
t = 0:dt:T_final;     % 时间向量
N = length(t);

theta0 = 15 * pi / 180;   % 初始角度 15°
x0 = [theta0; 0];         % 初始状态

% 闭环系统: dx/dt = (A - B*K) x
Acl = A - B * K_lqr;

% 使用 ode45 高精度求解（MATLAB 的 Runge-Kutta (4,5) 变步长）
[t_ode, x_ode] = ode45(@(t, x) Acl * x, [0 T_final], x0);

% 同时用固定步长 RK4（与 C++ 一致），方便对照
x_rk4 = zeros(N, 2);
x_rk4(1, :) = x0';
x = x0;

for k = 1:N-1
    % RK4
    k1 = Acl * x;
    k2 = Acl * (x + 0.5*dt * k1);
    k3 = Acl * (x + 0.5*dt * k2);
    k4 = Acl * (x + dt * k3);
    x = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    x_rk4(k+1, :) = x';
end

% 控制力: u = -K * x
u_rk4 = -K_lqr * x_rk4';

fprintf('\n── 仿真结果 ──\n');
fprintf('  RK4  稳态 θ(5s) = %.6f°\n', x_rk4(end,1) * 180/pi);
fprintf('  ode45 稳态 θ(5s) = %.6f°\n', x_ode(end,1) * 180/pi);
fprintf('  RK4 vs ode45 ‖Δθ‖_∞ = %.2e°\n', ...
        max(abs(x_rk4(:,1) - interp1(t_ode, x_ode(:,1), t'))) * 180/pi);

%% ── 5. 读取 C++ 输出进行验证 ────────────────────────
fprintf('\n── C++ 与 MATLAB 交叉验证 ──\n');

cpp_gain_found = false;
cpp_traj_found = false;
cpp_P_found = false;

% 5a. 读取 C++ 的 LQR 增益 (纯数值 CSV, 无表头, 1×2)
if isfile('lqr_gain_cpp.csv')
    cpp_gain = readmatrix('lqr_gain_cpp.csv');
    K_cpp = cpp_gain(1, :);  % 唯一一行: [K1, K2]
    fprintf('  C++ 反馈增益 K = [ %.8f  %.8f ]\n', K_cpp(1), K_cpp(2));
    fprintf('  MATLAB 反馈增益  = [ %.8f  %.8f ]\n', K_lqr(1), K_lqr(2));
    fprintf('  ‖K_cpp - K_matlab‖ = %.2e\n', norm(K_cpp - K_lqr));
    if norm(K_cpp - K_lqr) < 1e-6
        fprintf('  ✅ K 增益一致！\n');
    else
        fprintf('  ⚠️  K 增益有偏差（检查特征向量排序）\n');
    end
    cpp_gain_found = true;
else
    fprintf('  ⚠️  未找到 lqr_gain_cpp.csv (请先运行 C++ 程序)\n');
end

% 5b. 读取 C++ 的 Riccati 矩阵
if isfile('riccati_P_cpp.csv')
    cpp_P = readmatrix('riccati_P_cpp.csv');
    fprintf('\n  C++ Riccati 解 P = \n'); disp(cpp_P);
    fprintf('  MATLAB Riccati 解 P = \n'); disp(P_lqr);
    fprintf('  ‖P_cpp - P_matlab‖ = %.2e\n', norm(cpp_P - P_lqr));
    if norm(cpp_P - P_lqr) < 5e-6
        fprintf('  ✅ P 矩阵一致！\n');
    else
        fprintf('  ⚠️  P 矩阵有偏差\n');
    end
    cpp_P_found = true;
else
    fprintf('  ⚠️  未找到 riccati_P_cpp.csv\n');
end

% 5c. 读取 C++ 仿真轨迹
if isfile('pendulum_sim.csv')
    cpp_data = readmatrix('pendulum_sim.csv');
    t_cpp = cpp_data(:, 1);
    theta_cpp = cpp_data(:, 2);      % rad
    thetadot_cpp = cpp_data(:, 3);   % rad/s

    % 插值到相同时间点进行逐点比较
    theta_ml = interp1(t', x_rk4(:,1), t_cpp);
    thetadot_ml = interp1(t', x_rk4(:,2), t_cpp);

    fprintf('\n  C++ 稳态 θ(5s)     = %.6f°\n', theta_cpp(end) * 180/pi);
    fprintf('  MATLAB RK4 稳态     = %.6f°\n', x_rk4(end,1) * 180/pi);
    fprintf('  ‖θ_cpp - θ_matlab‖_∞ = %.2e°\n', ...
            max(abs(theta_cpp - theta_ml)) * 180/pi);
    fprintf('  ‖θ̇_cpp - θ̇_matlab‖_∞ = %.2e rad/s\n', ...
            max(abs(thetadot_cpp - thetadot_ml)));

    if max(abs(theta_cpp - theta_ml)) < 1e-8
        fprintf('  ✅ 仿真轨迹一致！\n');
    else
        fprintf('  ⚠️  仿真轨迹有微小偏差（浮点舍入误差可接受）\n');
    end
    cpp_traj_found = true;
else
    fprintf('  ⚠️  未找到 pendulum_sim.csv\n');
end

%% ── 6. 绘图 ────────────────────────────────────────
figure('Name', '一阶倒立摆 LQR 控制', 'Position', [100 100 1200 800]);

% ── 子图 1: 角度 vs 时间 ──
subplot(2, 3, 1);
plot(t, x_rk4(:,1) * 180/pi, 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--');
xlabel('时间 (s)'); ylabel('角度 \theta (°)');
title('摆杆角度'); grid on;
legend('LQR 控制', 'Location', 'northeast');

% ── 子图 2: 角速度 vs 时间 ──
subplot(2, 3, 2);
plot(t, x_rk4(:,2), 'r-', 'LineWidth', 1.5); hold on;
yline(0, 'k--');
xlabel('时间 (s)'); ylabel('角速度 \theta'' (rad/s)');
title('摆杆角速度'); grid on;

% ── 子图 3: 控制力 vs 时间 ──
subplot(2, 3, 3);
plot(t, u_rk4, 'Color', [0 0.5 0], 'LineWidth', 1.5); hold on;
yline(0, 'k--');
xlabel('时间 (s)'); ylabel('控制力 F (N)');
title('控制输入'); grid on;

% ── 子图 4: 相平面 (θ vs θ̇) ──
subplot(2, 3, 4);
plot(x_rk4(:,1) * 180/pi, x_rk4(:,2), 'm-', 'LineWidth', 1.5); hold on;
plot(x_rk4(1,1) * 180/pi, x_rk4(1,2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot(x_rk4(end,1) * 180/pi, x_rk4(end,2), 'g*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('\theta (°)'); ylabel('\theta'' (rad/s)');
title('相平面轨迹'); grid on;
legend('轨迹', '起点', '终点', 'Location', 'best');

% ── 子图 5: C++ vs MATLAB 对比 ──
subplot(2, 3, 5);
if cpp_traj_found
    plot(t_cpp, (theta_cpp - theta_ml) * 180/pi, 'k-', 'LineWidth', 1);
    xlabel('时间 (s)'); ylabel('\Delta\theta (°)');
    title('C++ vs MATLAB 角度差异'); grid on;
else
    text(0.5, 0.5, '无 C++ 数据\n(请运行 C++ 程序)', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    title('C++ vs MATLAB (无数据)');
end

% ── 子图 6: 闭环极点 ──
subplot(2, 3, 6);
plot(real(cl_eig), imag(cl_eig), 'rx', 'MarkerSize', 12, 'LineWidth', 2); hold on;
plot(real(eig(A)), imag(eig(A)), 'bo', 'MarkerSize', 10, 'LineWidth', 1.5);
xline(0, 'k--'); yline(0, 'k--');
xlabel('Real'); ylabel('Imag');
title('开环 vs 闭环极点');
legend('闭环极点 (LQR)', '开环极点', 'Location', 'best');
grid on; axis equal;

sgtitle('一阶线性倒立摆 LQR 控制 — MATLAB 验证');

%% ── 7. 权重参数扫描（可选分析） ────────────────────
fprintf('\n── 权重参数分析 ──\n');

% 固定 R, 变化 Q(1,1) / Q(2,2) 比值
q_ratios = [1, 10, 50, 100, 500, 1000];
fprintf('  R=%.3f, 变化 Q(1,1)/Q(2,2) 比值:\n', R);
fprintf('  %-10s %-20s %-20s\n', 'Q_ratio', 'K(1)', 'K(2)');
for qi = 1:length(q_ratios)
    Q_test = diag([q_ratios(qi), 1]);
    [K_test, ~, ~] = lqr(A, B, Q_test, R);
    fprintf('  %-10d %-20.8f %-20.8f\n', q_ratios(qi), K_test(1), K_test(2));
end

% 固定 Q, 变化 R
r_values = [0.001, 0.01, 0.1, 1.0, 10.0];
fprintf('\n  Q=diag([100,1]), 变化 R:\n');
fprintf('  %-10s %-20s %-20s\n', 'R', 'K(1)', 'K(2)');
for ri = 1:length(r_values)
    [K_test, ~, ~] = lqr(A, B, Q, r_values(ri));
    fprintf('  %-10.3f %-20.8f %-20.8f\n', r_values(ri), K_test(1), K_test(2));
end

fprintf('\n✅ MATLAB 验证完成！\n');
fprintf('  解读: Q(1,1)↑ → 角度偏差惩罚↑ → K(1)↑ → 更快回到平衡位置\n');
fprintf('        R↑ → 控制代价↑ → K 幅值↓ → 更保守的控制，调节更慢\n');
