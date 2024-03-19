//
// Created by RetrO on 2024/1/23.
//

#ifndef MATRIX_MAT_HPP
#define MATRIX_MAT_HPP

#include <omp.h>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>

template <typename T>
class Mat {
   private:
    int64_t m_rows{}, m_cols{};
    T* m_data;

   public:
    Mat() : m_rows(0), m_cols(0), m_data(nullptr) {}

    Mat(int64_t rows, int64_t cols) : m_rows(rows), m_cols(cols) {
        m_data = new T[cols * rows];
        memset(m_data, 0, cols * rows);
    }

    Mat(int64_t rows, int64_t cols, T value) : m_rows(rows), m_cols(cols) {
        m_data = new T[rows * cols];
        for (int64_t i = 0; i < rows * cols; ++i) {
            m_data[i] = value;
        }
    }

    Mat(int64_t rows, int64_t cols, T* value)
        : m_rows(rows), m_cols(cols), m_data(value) {}

    Mat(int64_t rows, int64_t cols, std::string str)
        : m_rows(rows), m_cols(cols) {
        if (str == "unity") {
            m_data = new T[cols * rows];
            for (int64_t i = 0; i < rows; i++) {
                for (int64_t j = 0; j < cols; j++) {
                    if (i == j) {
                        m_data[i * cols + j] = 1.0;
                    } else {
                        m_data[i * cols + j] = 0.0;
                    }
                }
            }
        } else {
            std::cerr << "invalid type" << std::endl;
        }
    }

    Mat(const Mat& other);

    Mat(Mat<T>&& other) noexcept;

    Mat& operator=(const Mat& other);

    Mat& operator=(Mat&& other) noexcept;

    template <typename Y>
    friend Mat<Y> operator+(const Mat<Y>& A, const Mat<Y>& B);

    template <typename Y, typename W>
    friend Mat<Y> operator+(const Mat<Y>& A, const W& B);

    template <typename Y, typename W>
    friend Mat<Y> operator+(const W& B, const Mat<Y>& A);

    template <typename Y>
    friend Mat<Y> operator-(const Mat<Y>& A, const Mat<Y>& B);

    template <typename Y, typename W>
    friend Mat<Y> operator-(const Mat<Y>& A, const W& B);

    template <typename Y, typename W>
    friend Mat<Y> operator-(const W& B, const Mat<Y>& A);

    template <typename Y>
    friend std::ostream& operator<<(std::ostream& os, const Mat<Y>& mat);

    template <typename Y>
    friend Mat<Y> operator*(const Mat<Y>& A, const Mat<Y>& B);

    T* operator[](int row) const;

    void switch_row(int64_t row1, int64_t row2);

    void switch_col(int64_t col1, int64_t col2);

    void transpose();

    Mat<T> trans();

    std::tuple<Mat<T>, Mat<T>, Mat<T>> decomposeLU();

    std::pair<Mat<T>, Mat<T>> decomposeCholesky();

    int64_t getRows() const;

    int64_t getCols() const;

    void printShape();

    ~Mat();
};

#endif  // MATRIX_MAT_HPP
