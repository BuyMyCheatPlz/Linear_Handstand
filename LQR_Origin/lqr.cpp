#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "mex.h"

class Matrix {
public:
    // 构造一个空矩阵，默认尺寸为 0x0。
    Matrix() : rows_(0), cols_(0) {}
    // 构造一个 rows x cols 的矩阵，并用 value 初始化所有元素。
    Matrix(int rows, int cols, double value = 0.0)
        : rows_(rows), cols_(cols), data_(static_cast<size_t>(rows) * cols, value) {}

    // 返回矩阵行数。
    int rows() const { return rows_; }
    // 返回矩阵列数。
    int cols() const { return cols_; }

    // 以可写方式访问指定位置的元素。
    double& operator()(int row, int col) { return data_[index(row, col)]; }
    // 以只读方式访问指定位置的元素。
    double operator()(int row, int col) const { return data_[index(row, col)]; }

    // 生成 n x n 单位矩阵。
    static Matrix identity(int n) {
        Matrix result(n, n, 0.0);
        for (int i = 0; i < n; ++i) {
            result(i, i) = 1.0;
        }
        return result;
    }

private:
    // 将二维下标转换成一维下标，并检查越界。
    size_t index(int row, int col) const {
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::out_of_range("Matrix index out of range");
        }
        return static_cast<size_t>(row) * cols_ + static_cast<size_t>(col);
    }

    int rows_;
    int cols_;
    std::vector<double> data_;
};

// 返回矩阵的转置。
Matrix transpose(const Matrix& m) {
    Matrix result(m.cols(), m.rows(), 0.0);
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            result(j, i) = m(i, j);
        }
    }
    return result;
}

// 对两个同型矩阵逐元素相加。
Matrix add(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("add: dimension mismatch");
    }
    Matrix result(a.rows(), a.cols(), 0.0);
    for (int i = 0; i < a.rows(); ++i) {
        for (int j = 0; j < a.cols(); ++j) {
            result(i, j) = a(i, j) + b(i, j);
        }
    }
    return result;
}

// 对两个同型矩阵逐元素相减。
Matrix subtract(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("subtract: dimension mismatch");
    }
    Matrix result(a.rows(), a.cols(), 0.0);
    for (int i = 0; i < a.rows(); ++i) {
        for (int j = 0; j < a.cols(); ++j) {
            result(i, j) = a(i, j) - b(i, j);
        }
    }
    return result;
}

// 计算矩阵乘法 a * b。
Matrix multiply(const Matrix& a, const Matrix& b) {
    if (a.cols() != b.rows()) {
        throw std::invalid_argument("multiply: dimension mismatch");
    }
    Matrix result(a.rows(), b.cols(), 0.0);
    for (int i = 0; i < a.rows(); ++i) {
        for (int k = 0; k < a.cols(); ++k) {
            const double aik = a(i, k);
            for (int j = 0; j < b.cols(); ++j) {
                result(i, j) += aik * b(k, j);
            }
        }
    }
    return result;
}

// 计算两个矩阵对应元素的最大绝对差，用于判断迭代是否收敛。
double maxAbsDiff(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("maxAbsDiff: dimension mismatch");
    }
    double diff = 0.0;
    for (int i = 0; i < a.rows(); ++i) {
        for (int j = 0; j < a.cols(); ++j) {
            diff = std::max(diff, std::fabs(a(i, j) - b(i, j)));
        }
    }
    return diff;
}

// 使用高斯-约旦消元求矩阵逆。
Matrix inverse(Matrix m) {
    if (m.rows() != m.cols()) {
        throw std::invalid_argument("inverse: matrix must be square");
    }

    const int n = m.rows();
    Matrix inv = Matrix::identity(n);

    for (int col = 0; col < n; ++col) {
        int pivotRow = col;
        double pivotAbs = std::fabs(m(col, col));
        for (int row =
             col + 1; row < n; ++row) {
            const double candidate = std::fabs(m(row, col));
            if (candidate > pivotAbs) {
                pivotAbs = candidate;
                pivotRow = row;
            }
        }

        if (pivotAbs < 1e-12) {
            throw std::runtime_error("inverse: singular matrix");
        }

        if (pivotRow != col) {
            for (int j = 0; j < n; ++j) {
                std::swap(m(col, j), m(pivotRow, j));
                std::swap(inv(col, j), inv(pivotRow, j));
            }
        }

        const double pivot = m(col, col);
        for (int j = 0; j < n; ++j) {
            m(col, j) /= pivot;
            inv(col, j) /= pivot;
        }

        for (int row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = m(row, col);
            if (std::fabs(factor) < 1e-15) {
                continue;
            }
            for (int j = 0; j < n; ++j) {
                m(row, j) -= factor * m(col, j);
                inv(row, j) -= factor * inv(col, j);
            }
        }
    }

    return inv;
}

// 用离散时间 Riccati 迭代求解 LQR 增益 K 和代价矩阵 P。
bool solveDiscreteLQR(const Matrix& A,
                      const Matrix& B,
                      const Matrix& Q,
                      const Matrix& R,
                      Matrix& K,
                      Matrix& P,
                      int maxIterations = 1000,
                      double tolerance = 1e-10) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("A must be square");
    }
    if (B.rows() != A.rows()) {
        throw std::invalid_argument("B row count must match A size");
    }
    if (Q.rows() != A.rows() || Q.cols() != A.cols()) {
        throw std::invalid_argument("Q must match A size");
    }
    if (R.rows() != B.cols() || R.cols() != B.cols()) {
        throw std::invalid_argument("R must be square with input dimension");
    }

    P = Q;
    Matrix A_T = transpose(A);
    Matrix B_T = transpose(B);

    for (int iter = 0; iter < maxIterations; ++iter) {
        Matrix BT_P = multiply(B_T, P);
        Matrix BT_P_B = multiply(BT_P, B);
        Matrix S = add(R, BT_P_B);
        Matrix S_inv = inverse(S);

        Matrix BT_P_A = multiply(BT_P, A);
        Matrix A_T_P = multiply(A_T, P);
        Matrix A_T_P_A = multiply(A_T_P, A);
        Matrix A_T_P_B = multiply(A_T_P, B);
        Matrix gainTerm = multiply(A_T_P_B, multiply(S_inv, BT_P_A));

        Matrix P_next = add(Q, subtract(A_T_P_A, gainTerm));

        if (maxAbsDiff(P_next, P) < tolerance) {
            P = P_next;
            K = multiply(S_inv, BT_P_A);
            return true;
        }

        P = P_next;
    }

    Matrix BT_P = multiply(B_T, P);
    Matrix S = add(R, multiply(BT_P, B));
    Matrix S_inv = inverse(S);
    K = multiply(S_inv, multiply(BT_P, A));
    return true;
}

// 将矩阵按整齐格式打印到标准输出。
void printMatrix(const Matrix& m, const std::string& name) {
    std::cout << name << " =\n";
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            std::cout << std::setw(12) << std::setprecision(8) << std::fixed << m(i, j) << ' ';
        }
        std::cout << '\n';
    }
}

// 演示入口：构造一个二阶离散系统并求解 LQR。
int main() {
    try {
        Matrix A(2, 2, 0.0);
        A(0, 0) = 1.0;
        A(0, 1) = 1.0;
        A(1, 0) = 0.0;
        A(1, 1) = 1.0;

        Matrix B(2, 1, 0.0);
        B(0, 0) = 0.5;
        B(1, 0) = 1.0;

        Matrix Q(2, 2, 0.0);
        Q(0, 0) = 10.0;
        Q(1, 1) = 1.0;

        Matrix R(1, 1, 0.0);
        R(0, 0) = 0.1;

        Matrix K;
        Matrix P;
        if (!solveDiscreteLQR(A, B, Q, R, K, P)) {
            std::cerr << "LQR 求解失败\n";
            return 1;
        }

        printMatrix(K, "K");
        printMatrix(P, "P");
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}