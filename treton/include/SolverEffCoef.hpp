#ifndef SOLVER_HEFF_HPP
#define SOLVER_HEFF_HPP

#include <numbers>
#include <cmath>
#include <stdexcept>


class EquationSolver {
public:
    enum class EquationType { Height, Radius };

    EquationSolver(double param1, double param2, EquationType type, double tol = 1e-6, int maxIter = 100);

    double solve(double left = 0.1, double right = 10);

private:
    double equation(double x) const;

    double K_eff_;
    double geom_size_;

    double tolerance_;
    int maxIterations_;

    EquationType type_;
};

#endif
