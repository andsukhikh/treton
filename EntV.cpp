#include "headers/EntV.hpp"

// Рассчет энтальпия натрия
double EntV(double temperature) {
    return -140954 + 1437 * temperature - 0.29 * temperature * temperature + 1.54 * 1E-04 * temperature * temperature * temperature;
}

