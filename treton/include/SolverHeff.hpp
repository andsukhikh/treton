#ifndef SOLVER_HEFF_HPP
#define SOLVER_HEFF_HPP

#include <numbers>
#include <cmath>
#include <stdexcept>


class EquationSolver {
private:
    double K_z;
    double H0;
    double tolerance;
    int maxIterations;

    double equation(double H) const;

public:
    EquationSolver(double K_z, double H0, double tol = 1e-6, int maxIter = 100);

    double solve(double H_left = 0.1, double H_right = 10);
};

#endif
