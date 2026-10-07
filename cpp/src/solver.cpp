// Conjugate gradients preconditioned by classical (Ruge-Stuben) algebraic
// multigrid, via AMGCL.  The operator is an M-matrix; strength-based coarsening
// follows the strong direction of the nanometre-to-micrometre cell anisotropy
// on its own, and the cost stays linear in the number of cells.
#include "stacktherm/stacktherm.hpp"

#include <chrono>
#include <stdexcept>
#include <tuple>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/solver/cg.hpp>

namespace st {

namespace {
using Backend = amgcl::backend::builtin<double>;
using Amg = amgcl::amg<Backend, amgcl::coarsening::ruge_stuben, amgcl::relaxation::gauss_seidel>;
using Solver = amgcl::make_solver<Amg, amgcl::solver::cg<Backend>>;

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
}  // namespace

struct LinearSolver::Impl {
    std::unique_ptr<Solver> solver;
    int n = 0;
    double setup_s = 0;
};

LinearSolver::LinearSolver(const Csr& A, double tol, int maxiter) : impl_(new Impl) {
    auto t0 = std::chrono::steady_clock::now();
    Solver::params prm;
    prm.solver.tol = tol;
    prm.solver.maxiter = maxiter;
    prm.precond.coarsening.eps_strong = 0.25f;
    impl_->n = A.n;
    impl_->solver.reset(new Solver(std::tie(A.n, A.ptr, A.col, A.val), prm));
    impl_->setup_s = seconds_since(t0);
}

LinearSolver::~LinearSolver() = default;

Vec LinearSolver::solve(const Vec& b, SolveInfo* info) const {
    auto t0 = std::chrono::steady_clock::now();
    Vec x(impl_->n, 0.0);
    bool zero = true;
    for (double v : b) if (v != 0.0) { zero = false; break; }
    int iterations = 0;
    double error = 0.0;
    if (!zero) {
        std::tie(iterations, error) = (*impl_->solver)(b, x);
        if (!(error < 1e-6)) throw std::runtime_error("the linear solver did not converge");
    }
    if (info) *info = SolveInfo{iterations, error, impl_->setup_s, seconds_since(t0)};
    return x;
}

}  // namespace st
