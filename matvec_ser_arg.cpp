#include <iostream>
#include <vector>
#include <iomanip>
#include <random>
#include <string>

int main(int argc, char** argv) {
    // Suponemos entrada correcta: argc==3, argv[1]=m, argv[2]=n
    const size_t m = std::stoull(std::string(argv[1]));
    const size_t n = std::stoull(std::string(argv[2]));

    // Generador aleatorio [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Matriz A (m x n) y vector x (n)
    std::vector<std::vector<double>> A(m, std::vector<double>(n));
    std::vector<double> x(n), y(m, 0.0);

    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            A[i][j] = dist(gen);

    for (size_t j = 0; j < n; ++j)
        x[j] = dist(gen);

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

// g++ -O3 -std=c++17 matvec_ser_arg.cpp -o matvec_ser_arg
// ./matvec_ser_arg 2 3