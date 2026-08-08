function run_lqr_sim()
    % 如果尚未生成 MEX，就先编译包装器。
    if ~exist('lqr_mex', 'file')
        fprintf('Compiling lqr_mex.cpp...\n');
        mex -R2018a lqr_mex.cpp
    end

    % 定义一个简单的离散系统，演示 LQR 求解流程。
    A = [1 1; 0 1];
    B = [0.5; 1.0];
    Q = diag([8, 1]);
    R = 10;

    % 通过 MEX 调用得到最优反馈增益和代价矩阵。
    [K, P] = lqr_mex(A, B, Q, R);

    fprintf('K =\n');
    disp(K);
    fprintf('P =\n');
    disp(P);

    % 使用 u = -Kx 做闭环仿真。
    x0 = [1; -1];
    N = 20;
    x = zeros(2, N + 1);
    u = zeros(1, N);
    x(:, 1) = x0;

    for k = 1:N
        u(k) = -K * x(:, k);
        x(:, k + 1) = A * x(:, k) + B * u(k);
    end

    % 画出状态和控制输入曲线。
    t = 0:N;
    figure('Name', 'LQR simulation');
    subplot(2, 1, 1);
    plot(t, x(1, :), '-o', t, x(2, :), '-s');
    xlabel('Step');
    ylabel('State');
    legend('x_1', 'x_2');
    grid on;

    subplot(2, 1, 2);
    stairs(t(1:end-1), u, 'r');
    xlabel('Step');
    ylabel('Control');
    title('u = -Kx');
    grid on;
end
