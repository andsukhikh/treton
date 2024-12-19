#include <cmath>

double CpSodium(double temperature) {
    return 1436 - 0.5805 * temperature + 4.62 * 1e-4 * temperature * temperature;
}
