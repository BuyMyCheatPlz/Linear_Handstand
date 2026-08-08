/**
 * lqr_mex_eigen.cpp — Eigen-based MEX wrapper for discrete-time LQR solver
 *
 * Usage in MATLAB:
 *   mex -R2018a lqr_mex_eigen.cpp -I%CPP_EIGEN%\include
 *   [K, P] = lqr_mex_eigen(A, B, Q, R);
 *
 * Dependencies: Eigen 3.x (header-only, set EIGEN_PATH or adjust -I flag)
 */

#include "mex.h"
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using Eigen::MatrixXd;

// --- MATLAB <-> Eigen conversion helpers ---

/** Convert a real double MATLAB array to an Eigen dynamic matrix.
 *  Eigen stores in column-major order by default, which matches MATLAB's
 *  Fortran-order layout — we can copy column by column efficiently.
 */
static MatrixXd mxArrayToEigen(const mxArray* arr) {
    if (!mxIsDouble(arr) || mxIsComplex(arr)) {
        mexErrMsgIdAndTxt("lqr_mex_eigen:invalidInput",
                          "Inputs must be real double arrays.");
    }

    const mwSize rows = mxGetM(arr);
    const mwSize cols = mxGetN(arr);
    const double* data = static_cast<const double*>(mxGetPr(arr));

    MatrixXd m(static_cast<int>(rows), static_cast<int>(cols));
    // MATLAB column-major: element (i,j) at offset i + j*rows
    for (mwSize j = 0; j < cols; ++j) {
        for (mwSize i = 0; i < rows; ++i) {
            m(static_cast<int>(i), static_cast<int>(j)) = data[i + j * rows];
        }
    }
    return m;
}

/** Convert an Eigen dynamic matrix back to a MATLAB double array. */
static mxArray* eigenToMxArray(const MatrixXd& m) {
    mxArray* out = mxCreateDoubleMatrix(
        static_cast<mwSize>(m.rows()),
        static_cast<mwSize>(m.cols()),
        mxREAL);
    double* data = static_cast<double*>(mxGetPr(out));
    for (int j = 0; j < m.cols(); ++j) {
        for (int i = 0; i < m.rows(); ++i) {
            data[i + j * m.rows()] = m(i, j);
        }
    }
    return out;
}

// --- Discrete-time LQR solver ---

/**
 * Solve the discrete-time LQR problem via Riccati iteration.
 *
 * System:  x_{k+1} = A x_k + B u_k
 * Cost:    J = Σ (x_k' Q x_k + u_k' R u_k)
 *
 * Discrete Algebraic Riccati Equation (DARE):
 *   P = Q + A' P A - A' P B (R + B' P B)^{-1} B' P A
 *
 * Optimal feedback gain:
 *   K = (R + B' P B)^{-1} B' P A
 *
 * @param[in]  A             State transition matrix (n × n)
 * @param[in]  B             Input matrix (n × m)
 * @param[in]  Q             State cost matrix (n × n), positive semi-definite
 * @param[in]  R             Input cost matrix (m × m), positive definite
 * @param[out] K             Optimal feedback gain (m × n)
 * @param[out] P             Solution to the DARE (n × n)
 * @param[in]  maxIterations Maximum Riccati iterations (default 1000)
 * @param[in]  tolerance     Convergence threshold on ||P_next - P||_∞
 * @return    true on success (always returns true; throws on invalid input)
 */
bool solveDiscreteLQR(const MatrixXd& A,
                      const MatrixXd& B,
                      const MatrixXd& Q,
                      const MatrixXd& R,
                      MatrixXd& K,
                      MatrixXd& P,
                      int    maxIterations = 1000,
                      double tolerance     = 1e-10) {
    // --- Dimension validation ---
    const int n = static_cast<int>(A.rows());   // state dimension
    const int m = static_cast<int>(B.cols());   // input dimension

    if (A.rows() != A.cols()) {
        throw std::invalid_argument("A must be square");
    }
    if (B.rows() != n) {
        throw std::invalid_argument("B row count must match A size");
    }
    if (Q.rows() != n || Q.cols() != n) {
        throw std::invalid_argument("Q must be n × n (match A)");
    }
    if (R.rows() != m || R.cols() != m) {
        throw std::invalid_argument("R must be m × m (match input dimension)");
    }

    // Pre-compute transposes (used every iteration)
    MatrixXd AT = A.transpose();
    MatrixXd BT = B.transpose();

    P = Q;

    for (int iter = 0; iter < maxIterations; ++iter) {
        // S = R + B' P B   (m × m symmetric positive definite)
        MatrixXd S = R + BT * P * B;

        // Riccati update:
        //   P_next = Q + A' P A - A' P B * S^{-1} * B' P A
        MatrixXd S_inv  = S.inverse();
        MatrixXd P_next = Q + AT * P * A
                          - (AT * P * B) * S_inv * (BT * P * A);

        // Check convergence in infinity-norm
        if ((P_next - P).cwiseAbs().maxCoeff() < tolerance) {
            P = P_next;
            // Compute gain from the converged P:
            //   K = S^{-1} B' P A
            K = S_inv * (BT * P * A);
            return true;
        }

        P = P_next;
    }

    // Max iterations reached — compute K with the final P estimate
    MatrixXd S     = R + BT * P * B;
    MatrixXd S_inv = S.inverse();
    K = S_inv * (BT * P * A);
    return true;
}

// --- MEX entry point ---------------------------------------------------------

/**
 * MEX gateway function.
 *
 * Inputs  (4 required):
 *   prhs[0] — A  (n × n)
 *   prhs[1] — B  (n × m)
 *   prhs[2] — Q  (n × n)
 *   prhs[3] — R  (m × m)
 *
 * Outputs (up to 2):
 *   plhs[0] — K  (m × n)  optimal gain
 *   plhs[1] — P  (n × n)  DARE solution (optional)
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("lqr_mex_eigen:nargin",
                          "Usage: [K, P] = lqr_mex_eigen(A, B, Q, R)");
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("lqr_mex_eigen:nargout",
                          "At most two outputs are supported.");
    }

    try {
        MatrixXd A = mxArrayToEigen(prhs[0]);
        MatrixXd B = mxArrayToEigen(prhs[1]);
        MatrixXd Q = mxArrayToEigen(prhs[2]);
        MatrixXd R = mxArrayToEigen(prhs[3]);

        MatrixXd K, P;
        if (!solveDiscreteLQR(A, B, Q, R, K, P)) {
            mexErrMsgTxt("LQR solve failed to converge.");
        }

        plhs[0] = eigenToMxArray(K);
        if (nlhs >= 2) {
            plhs[1] = eigenToMxArray(P);
        }
    } catch (const std::exception& ex) {
        mexErrMsgTxt(ex.what());
    }
}
