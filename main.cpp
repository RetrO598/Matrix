#include <iostream>
#include "src/iterativeSolver.cpp"
#include "src/mat.cpp"

int main() {
    // double valueA[6] = {1., 2., 3.2, 4., 2.1, 3.4};
    // double valueB[2] = {2.2, 3.1};
    // Mat<double> A(3, 2, valueA);
    // Mat<double> B(2, 1, valueB);
    // Mat<double> C;
    // C = A * B;
    // Mat<double> D(1, 2, 1);
    // std::cout << A << std::endl;
    // std::cout << B << std::endl;
    // std::cout << C << std::endl;
    // std::cout << D << std::endl;
    // double value[9] = {10, -1, 0, -1, 10, -2, 0, -2, 5};
    // double b[3] = {9, -5, 12};
    // Mat<double> A(3, 3, value);
    // Mat<double> B(3, 1, b);
    // std::cout << A << std::endl;
    // std::cout << B << std::endl;
    // std::cout << norm_2(B) << std::endl;
    // std::cout << gauss_seidel(A, B, 20) << std::endl;
    // Mat<double> x = gauss_seidel(A, B, 20);
    // std::cout << A * x << std::endl;

    double value[16] = {4, 2, 8, 0, 2, 10, 10, 9, 8, 10, 21, 6, 0, 9, 6, 34};
    Mat<double> A(4, 4, value);
    auto [L, Lt] = A.decomposeCholesky();
    std::cout << L << std::endl;
    std::cout << Lt << std::endl;
    std::cout << L * Lt << std::endl;
    return 0;
}
