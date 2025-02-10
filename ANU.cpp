#include <cmath>
#include "headers/ANU.hpp"

double ANU(double RE, double PR, double& A) {
    // Calculation of Nusselt number for water
    if (RE <= 2300.0) {
        A = 3.0;
        return 3.66;
    } else if (RE > 10000.0) {
        A = 1.0;
        auto g = std::pow(PR, 0.33);
        auto gg = std::pow(RE, 0.87);
        return 0.012 * std::pow(PR, 0.33) * std::pow(RE, 0.87);
    } else {
        A = 2.0;
        return 4.82 + 0.0185 * std::pow(PR, 0.33) * std::pow(RE, 0.8);
    }
}
