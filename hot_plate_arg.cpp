#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>

int main(int argc, char** argv) {
    const int N = std::stoi(std::string(argv[1]));
    const double TH = 0.1;           // Umbral "caliente"
    const double SOURCE = 100.0;     // Temperatura de la fuente central
    using Grid = std::vector<std::vector<double>>;

    Grid cur(N, std::vector<double>(N, 0.0));
    Grid nxt(N, std::vector<double>(N, 0.0));

    // Generador aleatorio [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Inicialización aleatoria del interior; bordes fijos en 0
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            cur[i][j] = dist(gen);
        }
    }

    // Fuente de calor fija en el centro
    const int ci = N / 2;
    const int cj = N / 2;
    cur[ci][cj] = SOURCE;
    nxt[ci][cj] = SOURCE;  // mantenerla fija tras la actualización

    // Actualización (una pasada) con la fórmula dada (solo interior)
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

    // Conteo de "hot cells" en el estado actualizado nxt:
    // 1 iff |(x,y) - ((x+1,y)+(x-1,y)+(x,y+1)+(x,y-1))/4| > 0.1
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

    std::cout << hot_count << "\n";
    return 0;
}

// g++ -O3 -std=c++14 hot_plate_arg.cpp -o hot_plate_arg
// ./hot_plate_arg 256