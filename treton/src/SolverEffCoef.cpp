#include "SolverEffCoef.hpp"

#include <numbers>
#include <cmath>
#include <stdexcept>

EquationSolver::EquationSolver(double param1, double param2, EquationType type, double tol, int maxIter)
    : K_eff_(param1)
    , geom_size_(param2)
    , tolerance_(tol)
    , maxIterations_(maxIter)
    , type_(type)
{}

double EquationSolver::equation(double x) const {
    if (type_ == EquationType::Height) {
        return K_eff_ * std::sin(std::numbers::pi * geom_size_ / (2 * x)) - (std::numbers::pi * geom_size_ / (2 * x));
    }
    else {
        return 2 * K_eff_ * std::cyl_bessel_j(1, 2.41 * geom_size_ / x) - 2.41 * (geom_size_ / x);
    }
}

double EquationSolver::solve(double left, double right) {
    if (equation(left) * equation(right) > 0) {
        throw std::invalid_argument("Root is not bracketted within the given range.");
    }

    double mid = 0.0;
    for (size_t i = 0; i < maxIterations_; ++i) {
        mid = (left + right) / 2.0;
        double f_mid = equation(mid);

        if (std::abs(f_mid) < tolerance_) {
            return mid;
        }

        if (equation(left) * f_mid < 0) {
            right = mid;
        }
        else {
            left = mid;
        }
    }
    return mid;
}
