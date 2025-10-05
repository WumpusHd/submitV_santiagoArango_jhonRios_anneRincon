#include <iostream>
#include <vector>
#include <iomanip>

int main() {
    // Matriz A (2x3) y vector x (3)
    std::vector<std::vector<double>> A = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0}
    };
    std::vector<double> x = {7.0, 8.0, 9.0};

    const size_t m = A.size();
    const size_t n = A[0].size();

    std::vector<double> y(m, 0.0);

    // y = A * x
    for (size_t i = 0; i < m; ++i) {
        double acc = 0.0;
        for (size_t j = 0; j < n; ++j) {
            acc += A[i][j] * x[j];
        }
        y[i] = acc;
    }

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);
    for (size_t i = 0; i < m; ++i) {
        std::cout << y[i] << (i + 1 == m ? '\n' : ' ');
    }
    return 0;
}
// g++ -std=c++17 -o matvec_ser matvec_ser.cpp
// ./matvec_ser