#ifndef IMATERIAL_HPP
#define IMATERIAL_HPP

#include <utility>

struct MaterialParams
{
    double density_ = 0.0;
    double therm_cond_ = 0.0;
    double heat_cap_ = 0.0;
};


class IMaterial
{
protected:
    virtual double Density(const double temperature) = 0;
    virtual double HeatCapacity(const double temperature) = 0;
    virtual double ThermalConductivity(const double temperature) = 0;

    virtual ~IMaterial() = default;

    IMaterial(MaterialParams&& params);

protected:
    MaterialParams params_;
};

#endif

