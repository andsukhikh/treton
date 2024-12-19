#include <cmath>
#include "headers/EntV.hpp"

void SodiumEV(const double& pressure, double& enthalpy, const double& temperature, double& volume) {
        enthalpy = -140954 + 1437 * temperature - 0.29 * temperature * temperature + 1.54 * 1E-04 * temperature * temperature * temperature;
        volume = 1 / (949 - 0.223 * temperature - 1.75 * 1E-5 * temperature * temperature);
}


void SodiumTV(const double& pressure, const double& enthalpy, double& temperature, double& volume) {       
        temperature = 0.00079 * enthalpy + 123.35;
        volume = 1 / (949 - 0.223 * temperature - 1.75 * 1E-5 * temperature * temperature);

}
