#include <exception>
#include <stdexcept>
#include <cmath>

#include "FuelMaterials.hpp"

double UO2::Density(double temperature)
{
    params_.density_ = 0.0;
    return params_.density_;
}

double UO2::HeatCapacity(double temperature)
{
    params_.heat_cap_ = 0.0;
    return params_.heat_cap_;
}

double UO2::ThermalConductivity(double temperature)
{
    if (temperature < 598)
    {
        params_.therm_cond_ = 0.0;
        return params_.therm_cond_;
    }
    if (temperature < 1327)
    {
        params_.therm_cond_ = 0.0;
        return params_.therm_cond_;
    }
    throw std::invalid_argument("temperature of coolant exceeded the permissible temperature of the approximation formula");
}
