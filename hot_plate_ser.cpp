#include <iostream>
#include <vector>
#include <cmath>

int main() {
    const int N = 256;                 // Tamaño de la placa: 256x256
    const double TH = 0.1;             // Umbral de “caliente”
    const double SOURCE = 100.0;       // Temperatura de la fuente central
    using Grid = std::vector<std::vector<double>>;

    Grid cur(N, std::vector<double>(N, 0.0));  // Estado inicial (0)
    Grid nxt(N, std::vector<double>(N, 0.0));  // Nuevo estado (una pasada)

    // Bordes fijos en 0.0 por construcción
    // Fuente de calor fija en el centro
    const int ci = N / 2;
    const int cj = N / 2;
    cur[ci][cj] = SOURCE;
    nxt[ci][cj] = SOURCE;  // Mantener fija la fuente al actualizar

    // Actualización (una pasada) usando la fórmula dada (solo celdas interiores)
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            if (i == ci && j == cj) continue; // la fuente permanece fija
            double up    = cur[i-1][j];
            double down  = cur[i+1][j];
            double left  = cur[i][j-1];
            double right = cur[i][j+1];
            double center= cur[i][j];
            nxt[i][j] = (up + down + left + right + 4.0 * center) / 8.0;
        }
    }

    // Contar “celdas calientes” con la condición indicada, sobre el estado actualizado
    int hot_count = 0;
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            double up    = nxt[i-1][j];
            double down  = nxt[i+1][j];
            double left  = nxt[i][j-1];
            double right = nxt[i][j+1];
            double avg_nb = (up + down + left + right) / 4.0;
            double diff = std::fabs(nxt[i][j] - avg_nb);
            if (diff > TH) ++hot_count;
        }
    }

    // Imprimir el número de celdas “calientes”
    std::cout << hot_count << "\n";
    return 0;
}

// g++ -O3 -std=c++14 hot_plate_ser.cpp -o hot_plate_ser
// ./hot_plate_ser