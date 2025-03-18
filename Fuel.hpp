#ifndef FUEL_HPP
#define FUEL_HPP

#include <vector>

#include "FuelMaterials.hpp"
#include "Mesh.hpp"


template<typename Material>
class Fuel : public Material, public Mesh
{
private:
    double inner_radius_;
    double outer_radius_;

public:
    Fuel(double irad, double orad)
        , inner_radius_(irad)
        , outer_radius_(orad)
    {}

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

    double area(double radius_first_boundary, double radius_second_boundary)
    {
        return std::numbers::pi * (std::pow(radius_second_boundary, 2) - std::pow(radius_first_boundary, 2));
    }

    void construct_mesh_from(double partitionNumber)
    {
        init_mesh<partitionNumber>(inner_radius_, outer_radius_);
    }
};

#endif 
