#ifndef ClADDING_HPP
#define ClADDING_HPP

#include <vector>
#include <cmath>
#include <numbers>

#include "CladdingMaterials.hpp"
#include "Mesh.hpp"


template<typename Material>
class Cladding : public Material
{
private:
	double outer_radius_;
	double inner_radius_;

public:
    Mesh coord;
    Mesh temperature;

    Cladding(double irad, double orad)
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
        mesh_.init_mesh<partitionNumber>(inner_radius_, outer_radius_);
    }
};

#endif

