#include <iostream>
#include <vector>
#include <iomanip>
#include <random>
#include <string>
#include <omp.h>

int main(int argc, char** argv) {
    const size_t m = std::stoull(std::string(argv[1]));
    const size_t n = std::stoull(std::string(argv[2]));

    // Generador aleatorio [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Matriz A (m x n), vector x (n) y resultado y (m)
    std::vector<std::vector<double>> A(m, std::vector<double>(n));
    std::vector<double> x(n), y(m, 0.0);

    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            A[i][j] = dist(gen);

    for (size_t j = 0; j < n; ++j)
        x[j] = dist(gen);

    // Paralelismo por filas (row-wise): cada hilo procesa bloques de filas
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(m); ++i) {
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
// g++ -fopenmp -O3 -std=c++14 matvec_omp.cpp -o matvec_omp
// OMP_NUM_THREADS=8 ./matvec_omp 2 3
// OMP_NUM_THREADS=8 ./matvec_omp 1000 2000