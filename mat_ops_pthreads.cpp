#include <pthread.h>
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>

using Matrix = std::vector<std::vector<double>>;

// ----------- Utilidades de matrices -----------
Matrix rand_matrix(int rows, int cols, std::mt19937& gen, std::uniform_real_distribution<double>& dist) {
    Matrix M(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M[i][j] = dist(gen);
    return M;
}

void print_matrix(const Matrix& M) {
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);
    int r = (int)M.size(), c = (int)M[0].size();
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j)
            std::cout << M[i][j] << (j + 1 == c ? '\n' : ' ');
    }
}

// ----------- Argumentos de hilos -----------
struct AddArgs {
    const Matrix* A;
    const Matrix* Bsum;
    Matrix* Csum;
    int m, n;
    pthread_barrier_t* barrier;
};

struct MulArgs {
    const Matrix* A;
    const Matrix* Bmul;
    Matrix* Cmul;
    int m, n, p;
    pthread_barrier_t* barrier;
};

struct TransArgs {
    const Matrix* A;
    Matrix* AT;
    int m, n;
    pthread_barrier_t* barrier;
};

// ----------- Tareas (hilos) -----------
void* task_add(void* arg) {
    auto* args = static_cast<AddArgs*>(arg);
    pthread_barrier_wait(args->barrier); // sincronizar inicio
    for (int i = 0; i < args->m; ++i)
        for (int j = 0; j < args->n; ++j)
            (*(args->Csum))[i][j] = (*(args->A))[i][j] + (*(args->Bsum))[i][j];
    return nullptr;
}

void* task_mul(void* arg) {
    auto* args = static_cast<MulArgs*>(arg);
    pthread_barrier_wait(args->barrier); // sincronizar inicio
    for (int i = 0; i < args->m; ++i) {
        for (int k = 0; k < args->n; ++k) {
            double aik = (*(args->A))[i][k];
            for (int j = 0; j < args->p; ++j) {
                (*(args->Cmul))[i][j] += aik * (*(args->Bmul))[k][j];
            }
        }
    }
    return nullptr;
}

void* task_transpose(void* arg) {
    auto* args = static_cast<TransArgs*>(arg);
    pthread_barrier_wait(args->barrier); // sincronizar inicio
    for (int i = 0; i < args->m; ++i)
        for (int j = 0; j < args->n; ++j)
            (*(args->AT))[j][i] = (*(args->A))[i][j];
    return nullptr;
}

int main() {
    // --- Entradas embebidas ---
    const int m = 3;   // filas de A
    const int n = 4;   // columnas de A (= filas de B_mul)
    const int p = 5;   // columnas de B_mul

    // RNG compartido para inicialización (en main, no en hilos)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Matrices de entrada
    Matrix A     = rand_matrix(m, n, gen, dist);
    Matrix Bsum  = rand_matrix(m, n, gen, dist); // para suma
    Matrix Bmul  = rand_matrix(n, p, gen, dist); // para multiplicación

    // Resultados
    Matrix Csum(m, std::vector<double>(n, 0.0)); // A + Bsum
    Matrix Cmul(m, std::vector<double>(p, 0.0)); // A * Bmul
    Matrix AT(n, std::vector<double>(m, 0.0));   // A^T

    // Barrera para arrancar las 3 tareas a la vez
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, nullptr, 3);

    // Crear hilos
    pthread_t t_add, t_mul, t_tr;
    AddArgs   aargs{&A, &Bsum, &Csum, m, n, &barrier};
    MulArgs   margs{&A, &Bmul, &Cmul, m, n, p, &barrier};
    TransArgs trargs{&A, &AT, m, n, &barrier};

    pthread_create(&t_add, nullptr, task_add, &aargs);
    pthread_create(&t_mul, nullptr, task_mul, &margs);
    pthread_create(&t_tr,  nullptr, task_transpose, &trargs);

    // Esperar finalización (barrera final implícita con join)
    pthread_join(t_add, nullptr);
    pthread_join(t_mul, nullptr);
    pthread_join(t_tr,  nullptr);

    pthread_barrier_destroy(&barrier);

    // --- Impresión (en main, para no mezclar salidas de hilos) ---
    std::cout << "A (" << m << "x" << n << "):\n";      print_matrix(A);
    std::cout << "B_sum (" << m << "x" << n << "):\n";  print_matrix(Bsum);
    std::cout << "C_sum = A + B_sum (" << m << "x" << n << "):\n"; print_matrix(Csum);

    std::cout << "B_mul (" << n << "x" << p << "):\n";  print_matrix(Bmul);
    std::cout << "C_mul = A * B_mul (" << m << "x" << p << "):\n"; print_matrix(Cmul);

    std::cout << "A^T (" << n << "x" << m << "):\n";    print_matrix(AT);

    return 0;
}

// g++ -O3 -std=c++17 -pthread mat_ops_pthreads.cpp -o mat_ops_pthreads
// ./mat_ops_pthreads