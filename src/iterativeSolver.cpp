#include "mat.hpp"

template <typename T>
T norm_2(const Mat<T>& x) {
    T result = 0;
    for (int64_t i = 0; i < x.getRows(); i++) {
        result += pow(x[i][0], 2);
    }
    return sqrt(result);
}

template <typename T>
Mat<T> jacobi(const Mat<T>& A,
              const Mat<T>& b,
              const int64_t& maxIter,
              const double& conv = 1e-6) {
    if (A.getRows() != b.getRows()) {
        std::cerr << "Illegal format of A and b"
                  << "\n";
    }
    int64_t rows = A.getRows();
    int64_t cols = A.getCols();
    Mat<T> xnew(rows, 1);
    Mat<T> xold(rows, 1);
    int iter = 0;
    T tmp = 0;
    T conv_ = 1;
    while (iter <= maxIter and conv_ >= conv) {
        for (int64_t i = 0; i < rows; i++) {
            for (int64_t j = 0; j < cols; j++) {
                tmp += A[i][j] * xold[j][0];
            }
            xnew[i][0] = b[i][0] / A[i][i] - tmp / A[i][i] + xold[i][0];
            tmp = 0;
        }
        conv_ = norm_2(xnew - xold);
        xold = xnew;
        iter++;
        if (iter == maxIter) {
            std::cout << "iter number reach the maximum iter number"
                      << std::endl;
        }
    }
    return xold;
}

template <typename T>
Mat<T> gauss_seidel(const Mat<T>& A,
                    const Mat<T>& b,
                    const int64_t& maxIter,
                    const double& conv = 1e-6) {
    if (A.getRows() != b.getRows()) {
        std::cerr << "Illegal format of A and b"
                  << "\n";
    }
    int64_t rows = A.getRows();
    int64_t cols = A.getCols();
    Mat<T> xnew(rows, 1);
    Mat<T> xold(rows, 1);
    int iter = 0;
    T tmp = 0;
    T conv_ = 1;
    while (iter <= maxIter and conv_ >= conv) {
        for (int64_t i = 0; i < rows; i++) {
            for (int64_t j = 0; j < cols; j++) {
                // tmp += A[i][j] * xold[j][0];
                if (j < i) {
                    tmp += A[i][j] * xnew[j][0];
                } else if (j > i) {
                    tmp += A[i][j] * xold[j][0];
                } else {
                    tmp += 0;
                }
            }
            xnew[i][0] = 1. / A[i][i] * (b[i][0] - tmp);
            tmp = 0;
        }
        conv_ = norm_2(xnew - xold);
        xold = xnew;
        iter++;
        if (iter == maxIter) {
            std::cout << "iter number reach the maximum iter number"
                      << std::endl;
        }
    }
    return xold;
}