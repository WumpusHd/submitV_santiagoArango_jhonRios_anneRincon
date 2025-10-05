#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <omp.h>

int main(int argc, char** argv) {
    const int N      = std::stoi(std::string(argv[1]));
    const int STEPS  = std::stoi(std::string(argv[2]));
    const double TH  = 0.1;        // Umbral "hot"
    const double SRC = 100.0;      // Fuente central
    const int BS     = 32;         // Tamaño de bloque (subdominio)

    using Grid = std::vector<std::vector<double>>;
    Grid A(N, std::vector<double>(N, 0.0));   // estado actual
    Grid B(N, std::vector<double>(N, 0.0));   // estado siguiente

    // RNG [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Interior aleatorio, bordes 0
    for (int i = 1; i < N - 1; ++i)
        for (int j = 1; j < N - 1; ++j)
            A[i][j] = dist(gen);

    // Fuente fija en el centro
    const int ci = N / 2, cj = N / 2;
    A[ci][cj] = SRC;
    B[ci][cj] = SRC;

    // Iteraciones Jacobi con descomposición por dominios 2D (tiles)
    #pragma omp parallel
    {
        for (int step = 0; step < STEPS; ++step) {

            // Reparto de bloques entre hilos
            #pragma omp for schedule(static)
            for (int bi = 1; bi < N - 1; bi += BS) {
                for (int bj = 1; bj < N - 1; bj += BS) {

                    const int i_end = std::min(bi + BS, N - 1);
                    const int j_end = std::min(bj + BS, N - 1);

                    // Actualizar celdas interiores del bloque usando vecinos de A
                    for (int i = bi; i < i_end; ++i) {
                        for (int j = bj; j < j_end; ++j) {
                            if (i == ci && j == cj) { B[i][j] = SRC; continue; }
                            double up    = A[i-1][j];
                            double down  = A[i+1][j];
                            double left  = A[i][j-1];
                            double right = A[i][j+1];
                            double center= A[i][j];
                            B[i][j] = (up + down + left + right + 4.0 * center) / 8.0;
                        }
                    }
                }
            }

            // Sincronización: asegurar que todos terminaron su subdominio
            #pragma omp barrier

            // Mantener fuente fija y hacer swap A<->B (una vez por iteración)
            #pragma omp single
            {
                std::swap(A, B);
                A[ci][cj] = SRC;
                B[ci][cj] = SRC;
            }

            // Barrera para que todos vean el nuevo A antes de la siguiente iteración
            #pragma omp barrier
        }
    }

    // Conteo de "hot cells" en el estado final A:
    // 1 iff |A(i,j) - (vecinos_prom)/4| > TH
    int hot_count = 0;
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            double avg_nb = (A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1]) / 4.0;
            if (std::fabs(A[i][j] - avg_nb) > TH) ++hot_count;
        }
    }

    std::cout << hot_count << "\n";
    return 0;
}
// g++ -fopenmp -O3 -std=c++14 hot_plate_domain_omp.cpp -o hot_plate_domain_omp
// Ejemplo: N=256, 100 iteraciones, 8 hilos
// OMP_NUM_THREADS=8 ./hot_plate_domain_omp 256 100