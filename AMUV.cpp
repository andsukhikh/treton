#include <cmath>
#include "headers/AMUV.hpp"

double AMUV(double temperature) {
    // Calculation of dynamic viscosity coefficient for water
    double density = 49 - 0.223 * temperature - 1.75 * 1e-5 * temperature * temperature;
    return 1.2162 * 1e-5 * pow(density, 1 / 3) * exp(0.6976 * density / temperature);
}
