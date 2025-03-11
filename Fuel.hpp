#ifndef FUEL_HPP
#define FUEL_HPP

#include <vector>

#include "FuelMaterials.hpp"


template<typename Material>
class Fuel : public Material
{
private:

    double outer_radius_;
    double inner_radius_;

    Fuel(double orad, double irad)
        : outer_radius_(orad)
        , inner_radius_(irad)
    {}
public:
    double get_outer_radius() const
    {
        return outer_radius_;
    }

    double get_inner_radius() const
    {
        return inner_radius_;
    }

    double get_volume_heat_capacity(double temperature)
    {
        return Density(temperature) * ThermalConductivity(temperature);
    }

    double get_thermal_conductivity(double temperature)
    {
        return ThermalConductivity(temperature);
    }

    double area()
    {
        return std::numbers::pi * (std::pow(outer_radius_, 2) - std::pow(inner_radius_, 2));
    }
};

#endif 
