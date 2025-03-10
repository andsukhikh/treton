#include <exception>
#include <stdexcept>
#include <cmath>

#include "CladdingMaterials.hpp"

double Zirconium::Density(double temperature)
{
    params_.density_ = 6340 - 0.15 * (temperature + 273 - 2128);
    return params_.density_;
}

double Zirconium::HeatCapacity(double temperature)
{
    params_.heat_cap_ = 0.245 + 0.5645 * 1e-4 * (temperature + 273);
    return params_.heat_cap_;
}

double Zirconium::ThermalConductivity(double temperature)
{
    if (temperature < 598)
    {
        params_.therm_cond_ = 16.49 + 8.7 * 1e-4 * (temperature + 273);
        return params_.therm_cond_;
    }
    if (temperature < 1327)
    {
        params_.therm_cond_ = -9.29 + 3.559 * 1e-2 * (temperature + 273) - 6.01 * 1e-6 * std::pow((temperature + 273), 2);
        return params_.therm_cond_;
    }
    throw std::invalid_argument("temperature of coolant exceeded the permissible temperature of the approximation formula");
}
