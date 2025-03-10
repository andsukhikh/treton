#ifndef CLADDING_MATERIALS_HPP
#define CLADDING_MATERIALS_HPP

#include "IMaterial.hpp"

class Zirconium : protected IMaterial
{
public:
    double Density(double temperature)                  override;
    double HeatCapacity(double temperature)             override;
    double ThermalConductivity(double temperature)      override;
};

#endif 


