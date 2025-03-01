#include "SolverHeff.hpp"

#include <numbers>
#include <cmath>
#include <stdexcept>


double EquationSolver::equation(double H) const {
    return K_z * std::sin(std::numbers::pi * H0 / (2 * H)) - (std::numbers::pi * H0 / (2 * H));
}


EquationSolver::EquationSolver(double K_z, double H0, double tol, int maxIter)
    : K_z(K_z), H0(H0), tolerance(tol), maxIterations(maxIter) {}


double EquationSolver::solve(double H_left, double H_right) {
    if (equation(H_left) * equation(H_right) > 0) {
        throw std::invalid_argument("Root is not bracK_zeted within the given range.");
    }

    double H_mid;
    for (int i = 0; i < maxIterations; ++i) {
        H_mid = (H_left + H_right) / 2.0;
        double f_mid = equation(H_mid);

        if (std::abs(f_mid) < tolerance) {
            return H_mid;
        }

        if (equation(H_left) * f_mid < 0) {
            H_right = H_mid;
        }
        else {
            H_left = H_mid;
        }
    }
    return H_mid;
}
