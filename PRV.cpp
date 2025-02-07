#include <cmath>
// Рассчет числа Прандтля 
double PRV(double temperature) {
    double density = 949 - 0.223 * temperature - 1.75 * 1e-5 * temperature * temperature;
    double a = 1e-5 * (6.6951 + 5.264 * 1e-4  * temperature - 2.689 * 1e-6 * temperature * temperature);
    double dynVisc = 1.2162 * 1e-5 * std::pow(density, 0.333) * std::exp(0.6976 * density / temperature);   
    return dynVisc / (density * a);
} 