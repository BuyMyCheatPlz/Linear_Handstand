function [K, P] = lqr_matlab_interface(A, B, Q, R)
%LQR_MATLAB_INTERFACE 通过 MEX 接口调用离散时间 LQR 求解器
%
%   [K, P] = lqr_matlab_interface(A, B, Q, R)
%
%   输入:
%       A - 状态转移矩阵，n x n
%       B - 输入矩阵，n x m
%       Q - 状态权重矩阵，n x n
%       R - 输入权重矩阵，m x m
%
%   输出:
%       K - 反馈增益矩阵，m x n
%       P - 代价矩阵，n x n
%
%   示例:
%       A = [1 1; 0 1];
%       B = [0.5; 1.0];
%       Q = diag([10, 1]);
%       R = 0.1;
%       [K, P] = lqr_matlab_interface(A, B, Q, R);

    % 如果 MEX 文件尚未编译，提前报错，避免后续调用失败。
    if ~exist('lqr_mex', 'file')
        error('MEX 文件 lqr_mex 尚未编译，请先运行 run_lqr_sim.m');
    end

    % 直接把求解工作交给 C++ 实现。
    [K, P] = lqr_mex(A, B, Q, R);
end
