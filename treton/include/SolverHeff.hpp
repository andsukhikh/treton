#ifndef SOLVER_HEFF_HPP
#define SOLVER_HEFF_HPP

#include <numbers>
#include <cmath>
#include <stdexcept>


class EquationSolver {
private:
    double K_z_;
    double Height0_;
    double tolerance_;
    int maxIterations_;

    double equation(double H) const;

public:
    EquationSolver(double K_z, double Height0, double tol = 1e-6, int maxIter = 100);

    double solve(double H_left = 0.1, double H_right = 10);
};

#endif
