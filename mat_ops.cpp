#include <iostream>
#include <vector>
#include <iomanip>
#include <random>

using Matrix = std::vector<std::vector<double>>;

Matrix rand_matrix(int rows, int cols, std::mt19937& gen, std::uniform_real_distribution<double>& dist) {
    Matrix M(rows, std::vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M[i][j] = dist(gen);
    return M;
}

Matrix add(const Matrix& A, const Matrix& B) {
    int m = (int)A.size();
    int n = (int)A[0].size();
    Matrix S(m, std::vector<double>(n, 0.0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            S[i][j] = A[i][j] + B[i][j];
    return S;
}

Matrix mul(const Matrix& A, const Matrix& B) {
    int m = (int)A.size();
    int n = (int)A[0].size();   // cols of A == rows of B
    int p = (int)B[0].size();
    Matrix C(m, std::vector<double>(p, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < n; ++k) {
            double aik = A[i][k];
            for (int j = 0; j < p; ++j) {
                C[i][j] += aik * B[k][j];
            }
        }
    }
    return C;
}

Matrix transpose(const Matrix& A) {
    int m = (int)A.size();
    int n = (int)A[0].size();
    Matrix AT(n, std::vector<double>(m, 0.0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            AT[j][i] = A[i][j];
    return AT;
}

void print_matrix(const Matrix& M) {
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);
    int m = (int)M.size();
    int n = (int)M[0].size();
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << M[i][j] << (j + 1 == n ? '\n' : ' ');
        }
    }
}

int main() {
    const int m = 3;   // filas de A
    const int n = 4;   // columnas de A (= filas de B_mul)
    const int p = 5;   // columnas de B_mul

    // RNG [0,1)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Matrices
    Matrix A      = rand_matrix(m, n, gen, dist);
    Matrix B_sum  = rand_matrix(m, n, gen, dist); // para suma
    Matrix B_mul  = rand_matrix(n, p, gen, dist); // para multiplicación

    Matrix S  = add(A, B_sum);   // A + B_sum  (m x n)
    Matrix M  = mul(A, B_mul);   // A * B_mul  (m x p)
    Matrix AT = transpose(A);    // A^T        (n x m)

    // Impresión
    std::cout << "A (" << m << "x" << n << "):\n";
    print_matrix(A);
    std::cout << "B_sum (" << m << "x" << n << "):\n";
    print_matrix(B_sum);
    std::cout << "S = A + B_sum (" << m << "x" << n << "):\n";
    print_matrix(S);

    std::cout << "B_mul (" << n << "x" << p << "):\n";
    print_matrix(B_mul);
    std::cout << "M = A * B_mul (" << m << "x" << p << "):\n";
    print_matrix(M);

    std::cout << "A^T (" << n << "x" << m << "):\n";
    print_matrix(AT);

    return 0;
}

// g++ -O3 -std=c++17 mat_ops.cpp -o mat_ops
// ./mat_ops