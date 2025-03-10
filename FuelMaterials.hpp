#ifndef FUEL_MATERIALS
#define FUEL_MATERIALS

#include "IMaterial.hpp"

class UO2 : protected IMaterial
{
public:
    double Density(double temperature)                  override;
    double HeatCapacity(double temperature)             override;
    double ThermalConductivity(double temperature)      override;
};


class UC : protected IMaterial
{
public:
    double Density(double temperature)                  override;
    double HeatCapacity(double temperature)             override;
    double ThermalConductivity(double temperature)      override;
};


class UN : protected IMaterial
{
public:
    double Density(double temperature)                  override;
    double HeatCapacity(double temperature)             override;
    double ThermalConductivity(double temperature)      override;
};


#endif 

