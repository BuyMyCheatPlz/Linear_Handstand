/**
 * lqr_eigen.cpp — Standalone Eigen-based discrete-time LQR solver (no MATLAB)
 *
 * Compile (MSVC):
 *   cl /EHsc /std:c++17 /I %CPP_EIGEN%\include lqr_eigen.cpp /Fe:lqr_eigen.exe
 *
 * Compile (GCC/MinGW):
 *   g++ -std=c++17 -O2 -I $CPP_EIGEN/include lqr_eigen.cpp -o lqr_eigen.exe
 *
 * Run:
 *   ./lqr_eigen.exe
 */

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using Eigen::MatrixXd;

// --- Pretty-print helper ---

void printMatrix(const MatrixXd& m, const std::string& name, int precision = 8) {
    std::cout << name << " =\n";
    std::cout << std::setprecision(precision) << std::fixed;
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            std::cout << std::setw(12) << m(i, j) << ' ';
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}

// --- LQR solver (same algorithm as lqr_mex_eigen.cpp) ---

bool solveDiscreteLQR(const MatrixXd& A,
                      const MatrixXd& B,
                      const MatrixXd& Q,
                      const MatrixXd& R,
                      MatrixXd& K,
                      MatrixXd& P,
                      int    maxIterations = 1000,
                      double tolerance     = 1e-10) {
    const int n = static_cast<int>(A.rows());
    const int m = static_cast<int>(B.cols());

    if (A.rows() != A.cols())
        throw std::invalid_argument("A must be square");
    if (B.rows() != n)
        throw std::invalid_argument("B row count must match A size");
    if (Q.rows() != n || Q.cols() != n)
        throw std::invalid_argument("Q must match A size");
    if (R.rows() != m || R.cols() != m)
        throw std::invalid_argument("R must be square with input dimension");

    MatrixXd AT = A.transpose();
    MatrixXd BT = B.transpose();

    P = Q;

    for (int iter = 0; iter < maxIterations; ++iter) {
        MatrixXd S      = R + BT * P * B;
        MatrixXd S_inv  = S.inverse();
        MatrixXd P_next = Q + AT * P * A
                          - (AT * P * B) * S_inv * (BT * P * A);

        if ((P_next - P).cwiseAbs().maxCoeff() < tolerance) {
            P = P_next;
            K = S_inv * (BT * P * A);
            std::cout << "[DARE] Converged in " << (iter + 1)
                      << " iterations.\n\n";
            return true;
        }

        P = P_next;
    }

    // Max iterations — finalize with current P
    MatrixXd S     = R + BT * P * B;
    MatrixXd S_inv = S.inverse();
    K = S_inv * (BT * P * A);
    std::cout << "[DARE] Reached max iterations (" << maxIterations
              << ").\n\n";
    return true;
}

// --- Demo ---

int main() {
    try {
        // ---- Example 1: Simple 2-state, 1-input system ----
        std::cout << "===== Example 1: 2-state integrator =====\n\n";

        MatrixXd A(2, 2);
        A << 1.0, 1.0,
             0.0, 1.0;

        MatrixXd B(2, 1);
        B << 0.5,
             1.0;

        MatrixXd Q(2, 2);
        Q << 10.0, 0.0,
              0.0, 1.0;

        MatrixXd R(1, 1);
        R << 0.1;

        printMatrix(A, "A");
        printMatrix(B, "B");
        printMatrix(Q, "Q");
        printMatrix(R, "R");

        // Solve
        MatrixXd K, P;
        solveDiscreteLQR(A, B, Q, R, K, P);

        printMatrix(K, "K (optimal gain)");
        printMatrix(P, "P (DARE solution)");

        // ---- Verify: closed-loop eigenvalues of (A - B K) ----
        Eigen::MatrixXd A_cl = A - B * K;
        Eigen::VectorXcd eig = A_cl.eigenvalues();
        std::cout << "Closed-loop eigenvalues (|λ| should be < 1):\n";
        for (int i = 0; i < eig.size(); ++i) {
            std::cout << "  λ" << i << " = " << eig(i)
                      << "  (|λ| = " << std::abs(eig(i)) << ")\n";
        }
        std::cout << std::endl;

        // ---- Verify: DARE residual ----
        MatrixXd S_check = R + B.transpose() * P * B;
        MatrixXd dare_residual =
            A.transpose() * P * A - P
            - A.transpose() * P * B * S_check.inverse() * B.transpose() * P * A
            + Q;
        double residual_norm = dare_residual.cwiseAbs().maxCoeff();
        std::cout << "||DARE residual||_∞ = " << residual_norm << "\n\n";

        // ---- Example 2: 3-state, 2-input system ----
        std::cout << "===== Example 2: 3-state, 2-input system =====\n\n";

        MatrixXd A2(3, 3);
        A2 << 1.0,  0.1,  0.0,
              0.0,  1.0,  0.1,
              0.0,  0.0,  1.0;

        MatrixXd B2(3, 2);
        B2 << 0.0,  0.0,
              0.1,  0.0,
              0.0,  0.1;

        MatrixXd Q2 = MatrixXd::Identity(3, 3) * 10.0;
        MatrixXd R2 = MatrixXd::Identity(2, 2) * 1.0;

        MatrixXd K2, P2;
        solveDiscreteLQR(A2, B2, Q2, R2, K2, P2);

        printMatrix(K2, "K (optimal gain)");
        printMatrix(P2, "P (DARE solution)");

        // Verify closed-loop stability
        Eigen::MatrixXd A_cl2 = A2 - B2 * K2;
        Eigen::VectorXcd eig2 = A_cl2.eigenvalues();
        std::cout << "Closed-loop eigenvalues:\n";
        for (int i = 0; i < eig2.size(); ++i) {
            std::cout << "  λ" << i << " = " << eig2(i)
                      << "  (|λ| = " << std::abs(eig2(i)) << ")\n";
        }

        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
