//
// Created by songx on 2024/1/26.
//
#include "mat.hpp"

template <typename T>
Mat<T>::Mat(const Mat& other) {
    m_rows = other.m_rows;
    m_cols = other.m_cols;
    m_data = new T[m_rows * m_cols];
    for (int64_t i = 0; i < m_cols * m_rows; ++i) {
        m_data[i] = other.m_data[i];
    }
}

template <typename T>
Mat<T>::Mat(Mat<T>&& other) noexcept {
    m_rows = other.m_rows;
    m_cols = other.m_cols;
    m_data = other.m_data;
    other.m_data = nullptr;
}

template <typename T>
Mat<T>& Mat<T>::operator=(const Mat<T>& other) {
    if (this != &other) {
        delete[] m_data;
        m_rows = other.m_rows;
        m_cols = other.m_cols;
        m_data = new T[m_cols * m_rows];
        for (int i = 0; i < m_rows * m_cols; ++i) {
            m_data[i] = other.m_data[i];
        }
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator=(Mat<T>&& other) noexcept {
    if (this != &other) {
        delete[] m_data;
        m_rows = other.m_rows;
        m_cols = other.m_cols;
        m_data = other.m_data;
        other.m_data = nullptr;
    }
    return *this;
}

template <typename T>
Mat<T> operator+(const Mat<T>& A, const Mat<T>& B) {
    if (A.m_rows == B.m_rows and A.m_cols == B.m_cols) {
        Mat<T> sum(A.m_rows, A.m_cols);
        for (int64_t i = 0; i < A.m_rows * A.m_cols; ++i) {
            sum.m_data[i] = A.m_data[i] + B.m_data[i];
        }
        return sum;
    } else {
        std::cerr << "The Dimensions of Two Matrix Should be Equal."
                  << std::endl;
    }
}

template <typename T, typename T2>
Mat<T> operator+(const Mat<T>& A, const T2& B) {
    Mat<T> sum(A.m_rows, A.m_cols);
    for (int i = 0; i < A.m_rows * A.m_cols; ++i) {
        sum.m_data[i] = A.m_data[i] + B;
    }
    return sum;
}

template <typename T, typename T2>
Mat<T> operator+(const T2& B, const Mat<T>& A) {
    Mat<T> sum(A.m_rows, A.m_cols);
    for (int i = 0; i < A.m_rows * A.m_cols; ++i) {
        sum.m_data[i] = A.m_data[i] + B;
    }
    return sum;
}

template <typename T>
Mat<T> operator-(const Mat<T>& A, const Mat<T>& B) {
    Mat<T> minus(A.m_rows, A.m_cols);
    if (A.m_rows == B.m_rows and A.m_cols == B.m_cols) {
        for (int i = 0; i < A.m_rows * A.m_cols; ++i) {
            minus.m_data[i] = A.m_data[i] - B.m_data[i];
        }
    } else {
        std::cerr << "The Dimensions of Two Matrix Should be Equal."
                  << std::endl;
    }
    return minus;
}

template <typename T, typename T2>
Mat<T> operator-(const Mat<T>& A, const T2& B) {
    Mat<T> sum(A.m_rows, A.m_cols);
    for (int i = 0; i < A.m_rows * A.m_cols; ++i) {
        sum.m_data[i] = A.m_data[i] - B;
    }
    return sum;
}

template <typename T, typename T2>
Mat<T> operator-(const T2& B, const Mat<T>& A) {
    Mat<T> sum(A.m_rows, A.m_cols);
    for (int i = 0; i < A.m_rows * A.m_cols; ++i) {
        sum.m_data[i] = B - A.m_data[i];
    }
    return sum;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Mat<T>& mat) {
    std::stringstream ss;
    ss.precision(2);
    ss.setf(std::ios::fixed);
    for (int64_t i = 0; i < mat.m_rows; i++) {
        for (int64_t j = 0; j < mat.m_cols; j++) {
            ss << mat.m_data[i * mat.m_cols + j] << " ";
        }
        ss << "\n";
    }
    os << ss.str() + "\n";
    return os;
}

template <typename T>
Mat<T> operator*(const Mat<T>& A, const Mat<T>& B) {
    int64_t rows = A.m_rows;
    int64_t cols = B.m_cols;
    if (A.m_cols != B.m_rows) {
        std::cerr << "Illegal product"
                  << "\n";
    }
    Mat<T> C(rows, cols);
#pragma omp parallel for collapse(2) shared(A, B, C)
    for (int64_t i = 0; i < rows; i++) {
        for (int64_t j = 0; j < cols; j++) {
            for (int64_t k = 0; k < A.m_cols; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

template <typename T>
T* Mat<T>::operator[](int row) const {
    return m_data + (row)*m_cols;
}

template <typename T>
void Mat<T>::switch_row(int64_t row1, int64_t row2) {
    double tmp;
    for (int64_t i = 0; i < m_cols; i++) {
        tmp = m_data[row1 * m_cols + i];
        m_data[row1 * m_cols + i] = m_data[row2 * m_cols + i];
        m_data[row2 * m_cols + i] = tmp;
    }
}

template <typename T>
void Mat<T>::switch_col(int64_t col1, int64_t col2) {
    double tmp;
    for (int64_t i = 0; i < m_rows; i++) {
        tmp = m_data[i * m_cols + col1];
        m_data[i * m_cols + col1] = m_data[i * m_cols + col2];
        m_data[i * m_cols + col2] = tmp;
    }
}

template <typename T>
void Mat<T>::transpose() {
    T* tmp = new T[m_cols * m_rows];

    for (int64_t i = 0; i < m_rows; i++) {
        for (int64_t j = 0; j < m_cols; j++) {
            std::cout << m_data[i * m_cols + j] << "\n";
            tmp[j * m_rows + i] = m_data[i * m_cols + j];
        }
    }
    for (int64_t i = 0; i < m_cols * m_rows; i++) {
        m_data[i] = tmp[i];
    }
    delete[] tmp;
    int64_t new_rows = m_cols;
    m_cols = m_rows;
    m_rows = new_rows;
}

template <typename T>
Mat<T> Mat<T>::trans() {
    Mat<T> transp(m_cols, m_rows);
    for (int64_t i = 0; i < m_rows; i++) {
        for (int64_t j = 0; j < m_cols; j++) {
            transp.m_data[j * m_rows + i] = this->m_data[i * m_cols + j];
        }
    }
    return transp;
}

template <typename T>
std::tuple<Mat<T>, Mat<T>, Mat<T>> Mat<T>::decomposeLU() {
    if (m_cols != m_rows) {
        std::cerr << "invalide decompose" << std::endl;
    }

    int64_t n = this->m_rows;
    Mat<T> L(n, n, "unity");
    Mat<T> U(n, n);
    Mat<T> A(n, n, this->m_data);
    Mat<T> P(n, n, "unity");
    for (int64_t k = 0; k < n - 1; k++) {
        int64_t l = k;
        double max = abs(A[k][k]);
        for (int64_t i = k; i < n; i++) {
            if (abs(A[i][k]) > max) {
                l = i;
                max = A[i][k];
            }
        }
        if (l != k) {
            A.switch_row(k, l);
            P.switch_row(k, l);
        }
    }

    for (int64_t k = 0; k < n; k++) {
        for (int64_t j = k; j < n; j++) {
            U[k][j] = A[k][j];
        }
        for (int64_t i = k + 1; i < n; i++) {
            L[i][k] = A[i][k] / A[k][k];
            for (int64_t j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - L[i][k] * U[k][j];
            }
        }
    }
    return {P, L, U};
}

template <typename T>
std::pair<Mat<T>, Mat<T>> Mat<T>::decomposeCholesky() {
    Mat<T> L(m_rows, m_cols);
    for (int64_t j = 0; j < m_rows; j++) {
        for (int64_t i = j; i < m_cols; i++) {
            if (i == j) {
                double sum = 0.0;
                for (int64_t k = 0; k < j; k++) {
                    // std::cout << k << std::endl;
                    sum += L[j][k] * L[j][k];
                }
                L[i][j] = sqrt(this->m_data[j * m_cols + j] - sum);
            } else {
                double sum = 0.0;
                for (int64_t k = 0; k < j; k++) {
                    std::cout << k << std::endl;
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = 1. / L[j][j] * (this->m_data[i * m_cols + j] - sum);
            }
        }
    }
    return {L, L.trans()};
}

template <typename T>
int64_t Mat<T>::getRows() const {
    return this->m_rows;
}

template <typename T>
int64_t Mat<T>::getCols() const {
    return this->m_cols;
}

template <typename T>
void Mat<T>::printShape() {
    std::cout << "(" << this->m_rows << ',' << this->m_cols << ")" << std::endl;
}

template <typename T>
Mat<T>::~Mat() {
    delete[] m_data;
    m_rows = 0;
    m_cols = 0;
}