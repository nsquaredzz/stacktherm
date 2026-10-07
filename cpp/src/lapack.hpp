// Thin wrappers over the Fortran BLAS / LAPACK symbols (Accelerate on macOS,
// OpenBLAS or the reference library elsewhere).  Matrices are column-major.
#pragma once

#include <stdexcept>
#include <vector>

extern "C" {
void dsyevd_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double* w, double* work,
             const int* lwork, int* iwork, const int* liwork, int* info);
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha,
            const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c,
            const int* ldc);
}

namespace st {

// eigenvalues ascending in w; a is overwritten with the eigenvectors (columns)
inline void symmetric_eigen(int n, std::vector<double>& a, std::vector<double>& w) {
    w.assign(n, 0.0);
    int info = 0, lwork = -1, liwork = -1, iquery = 0;
    double query = 0;
    dsyevd_("V", "U", &n, a.data(), &n, w.data(), &query, &lwork, &iquery, &liwork, &info);
    lwork = int(query);
    liwork = iquery;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    dsyevd_("V", "U", &n, a.data(), &n, w.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
    if (info != 0) throw std::runtime_error("eigenvalue decomposition failed");
}

// c = op(a) * op(b) for square n x n matrices; ta / tb are 'N' or 'T'
inline void matmul(int n, char ta, const std::vector<double>& a, char tb, const std::vector<double>& b,
                   std::vector<double>& c) {
    const double one = 1.0, zero = 0.0;
    c.resize(size_t(n) * n);
    dgemm_(&ta, &tb, &n, &n, &n, &one, a.data(), &n, b.data(), &n, &zero, c.data(), &n);
}

}  // namespace st
