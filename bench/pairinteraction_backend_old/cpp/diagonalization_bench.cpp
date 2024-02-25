#include <benchmark/benchmark.h>
#include <eigen3/Eigen/Dense>
#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

#include <cassert>

constexpr int const N = 333;

static void BM_zheev(benchmark::State &state) {
    Eigen::MatrixXcd hamiltonian = Eigen::MatrixXcd::Random(N, N);
    hamiltonian = (hamiltonian + hamiltonian.adjoint()).eval();
    assert(hamiltonian == hamiltonian.adjoint());
    for (auto _ : state) {
        Eigen::MatrixXcd m_eivec = hamiltonian;
        lapack_int n = m_eivec.cols();
        Eigen::VectorXd m_eivalues(n);
        lapack_int lda = m_eivec.outerStride();
        char jobz = 'V';
        char uplo = 'L';
        lapack_complex_double *A = (lapack_complex_double *)m_eivec.data();
        double *W = m_eivalues.data();
        lapack_int info = LAPACKE_zheev(LAPACK_COL_MAJOR, jobz, uplo, n, A, lda, W);
        assert(info == 0);
        benchmark::DoNotOptimize(A);
        benchmark::DoNotOptimize(W);
    }
}

static void BM_zheevd(benchmark::State &state) {
    Eigen::MatrixXcd hamiltonian = Eigen::MatrixXcd::Random(N, N);
    hamiltonian = (hamiltonian + hamiltonian.adjoint()).eval();
    assert(hamiltonian == hamiltonian.adjoint());
    for (auto _ : state) {
        Eigen::MatrixXcd m_eivec = hamiltonian;
        lapack_int n = m_eivec.cols();
        Eigen::VectorXd m_eivalues(n);
        lapack_int lda = m_eivec.outerStride();
        char jobz = 'V';
        char uplo = 'L';
        lapack_complex_double *A = (lapack_complex_double *)m_eivec.data();
        double *W = m_eivalues.data();
        lapack_int info = LAPACKE_zheevd(LAPACK_COL_MAJOR, jobz, uplo, n, A, lda, W);
        assert(info == 0);
        benchmark::DoNotOptimize(A);
        benchmark::DoNotOptimize(W);
    }
}

static void BM_Eigen(benchmark::State &state) {
    Eigen::MatrixXcd hamiltonian = Eigen::MatrixXcd::Random(N, N);
    hamiltonian = (hamiltonian + hamiltonian.adjoint()).eval();
    assert(hamiltonian == hamiltonian.adjoint());
    for (auto _ : state) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(hamiltonian);
        benchmark::DoNotOptimize(eigensolver.eigenvalues());
        benchmark::DoNotOptimize(eigensolver.eigenvectors());
    }
}

BENCHMARK(BM_zheev);
BENCHMARK(BM_zheevd);
BENCHMARK(BM_Eigen);
BENCHMARK_MAIN();
