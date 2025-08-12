#ifndef AMG_SOLVER_H
#define AMG_SOLVER_H

#include <vector>
#include <algorithm>
#include <cmath>

// Utility vector types
struct Vec2i {
    int x, y;
    Vec2i() : x(0), y(0) {}
    Vec2i(int x_, int y_) : x(x_), y(y_) {}
};

// Fixed sparse matrix structure for AMG solver
template<class T>
class FixedSparseMatrix {
public:
    struct Entry {
        int row, col;
        T value;
        Entry(int r, int c, T v) : row(r), col(c), value(v) {}
    };
    
    std::vector<Entry> entries;
    int rows, cols;
    
    FixedSparseMatrix(int r, int c) : rows(r), cols(c) {}
    
    void addEntry(int row, int col, T value) {
        entries.push_back(Entry(row, col, value));
    }
    
    void clear() {
        entries.clear();
    }
};

// BLAS-like operations namespace
namespace BLAS {
    template<class T>
    void zero(std::vector<T>& x) {
        std::fill(x.begin(), x.end(), T(0));
    }
    
    template<class T>
    T dot(const std::vector<T>& x, const std::vector<T>& y) {
        T result = T(0);
        for (size_t i = 0; i < x.size(); ++i) {
            result += x[i] * y[i];
        }
        return result;
    }
    
    template<class T>
    T abs_max(const std::vector<T>& x) {
        T result = T(0);
        for (size_t i = 0; i < x.size(); ++i) {
            result = std::max(result, std::abs(x[i]));
        }
        return result;
    }
    
    template<class T>
    T mean(const std::vector<T>& x) {
        T sum = T(0);
        for (size_t i = 0; i < x.size(); ++i) {
            sum += x[i];
        }
        return sum / T(x.size());
    }
    
    template<class T>
    void add_scaled(T alpha, const std::vector<T>& x, std::vector<T>& y) {
        for (size_t i = 0; i < x.size() && i < y.size(); ++i) {
            y[i] += alpha * x[i];
        }
    }
    
    template<class T>
    void subtractConst(std::vector<T>& x, T c) {
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] -= c;
        }
    }
}

// Zero vector utility
template<class T>
void zero(std::vector<T>& x) {
    BLAS::zero(x);
}

// Matrix-vector multiplication
template<class T>
void multiply(const FixedSparseMatrix<T>& matrix, const std::vector<T>& x, std::vector<T>& y) {
    BLAS::zero(y);
    for (const auto& entry : matrix.entries) {
        y[entry.row] += entry.value * x[entry.col];
    }
}

// Simple AMG preconditioner (placeholder implementation)
template<class T>
void amgPrecond2D(std::vector<FixedSparseMatrix<T>*>& A_L,
                  std::vector<FixedSparseMatrix<T>*>& R_L,
                  std::vector<FixedSparseMatrix<T>*>& P_L,
                  std::vector<Vec2i>& S_L,
                  std::vector<T>& z,
                  const std::vector<T>& r) {
    // Simple Jacobi preconditioner as placeholder
    z = r;
    for (size_t i = 0; i < z.size(); ++i) {
        z[i] *= T(0.6);  // Simple damping factor
    }
}

// Poisson matrix setup for 2D grid
template<class T>
void setupPoissonMatrix2D(FixedSparseMatrix<T>& matrix, int ni, int nj) {
    matrix.clear();
    
    for (int j = 0; j < nj; ++j) {
        for (int i = 0; i < ni; ++i) {
            int idx = j * ni + i;
            
            // Diagonal entry
            T diag_val = T(4);
            matrix.addEntry(idx, idx, diag_val);
            
            // Off-diagonal entries
            if (i > 0) {
                matrix.addEntry(idx, idx - 1, T(-1));  // Left
            }
            if (i < ni - 1) {
                matrix.addEntry(idx, idx + 1, T(-1));  // Right
            }
            if (j > 0) {
                matrix.addEntry(idx, idx - ni, T(-1)); // Bottom
            }
            if (j < nj - 1) {
                matrix.addEntry(idx, idx + ni, T(-1)); // Top
            }
        }
    }
}

// AMG PCG Solver (main function provided by user)
template<class T>
bool AMGPCGSolvePrebuilt2D(const FixedSparseMatrix<T> &fixed_matrix,
                           const std::vector<T> &rhs,
                           std::vector<T> &result,
                           std::vector<FixedSparseMatrix<T> *> &A_L,
                           std::vector<FixedSparseMatrix<T> *> &R_L,
                           std::vector<FixedSparseMatrix<T> *> &P_L,
                           std::vector<Vec2i> &S_L,
                           const int total_level,
                           T tolerance_factor,
                           int max_iterations,
                           T &residual_out,
                           int &iterations_out,
                           int ni, int nj,
                           bool PURE_NEUMANN) {
    std::vector<T> m, z, s, r, v;
    unsigned int n = ni * nj;
    if (m.size() != n) {
        m.resize(n);
        s.resize(n);
        z.resize(n);
        r.resize(n);
        if (PURE_NEUMANN) { v.resize(n); }
    }
    zero(result);
    r = rhs;
    if (PURE_NEUMANN) {
        T mean = BLAS::mean(r);
        v = r;
        BLAS::subtractConst(v, mean);
    }
    residual_out = PURE_NEUMANN ? BLAS::abs_max(v) : BLAS::abs_max(r);
    if (residual_out == 0) {
        iterations_out = 0;
        return true;
    }
    T tol = std::max(tolerance_factor * residual_out, tolerance_factor);
    if (PURE_NEUMANN) r = v;
    amgPrecond2D(A_L, R_L, P_L, S_L, z, r);
    T rho = BLAS::dot(z, r);
    if (rho == 0 || rho != rho) {
        iterations_out = 0;
        return false;
    }

    s = z;

    int iteration;
    for (iteration = 0; iteration < max_iterations; ++iteration) {
        multiply(fixed_matrix, s, z);
        T alpha = rho / BLAS::dot(s, z);
        BLAS::add_scaled(alpha, s, result);
        BLAS::add_scaled(-alpha, z, r);
        if (PURE_NEUMANN) {
            T mean = BLAS::mean(r);
            v = r;
            BLAS::subtractConst(v, mean);
        }
        residual_out = PURE_NEUMANN ? BLAS::abs_max(v) : BLAS::abs_max(r);
        if (residual_out <= tol) {
            iterations_out = iteration + 1;
            return true;
        }
        if (PURE_NEUMANN) r = v;
        amgPrecond2D(A_L, R_L, P_L, S_L, z, r);
        T rho_new = BLAS::dot(z, r);
        T beta = rho_new / rho;
        BLAS::add_scaled(beta, s, z);
        s.swap(z); // s=beta*s+z
        rho = rho_new;
    }
    iterations_out = iteration;
    return false;
}

#endif // AMG_SOLVER_H
