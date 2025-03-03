#include "SolverHeff.hpp"

#include <numbers>
#include <cmath>
#include <stdexcept>


double EquationSolver::equation(double Height) const {
    return K_z_ * std::sin(std::numbers::pi * Height0_ / (2 * Height)) - (std::numbers::pi * Height0_ / (2 * Height));
}


EquationSolver::EquationSolver(double K_z, double Height0, double tol, int maxIter)
    : K_z_(K_z), Height0_(Height0), tolerance_(tol), maxIterations_(maxIter) {}


double EquationSolver::solve(double H_left, double H_right) {
    if (equation(H_left) * equation(H_right) > 0) {
        throw std::invalid_argument("Root is not bracK_zeted within the given range.");
    }

    double H_mid = 0.0;
    for (size_t i = 0; i < maxIterations_; ++i) {
        H_mid = (H_left + H_right) / 2.0;
        double f_mid = equation(H_mid);

        if (std::abs(f_mid) < tolerance_) {
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
