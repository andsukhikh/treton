#ifndef ROD_HPP
#define ROD_HPP

#include "Gap.hpp"
#include "Cladding.hpp"
#include "Fuel.hpp"

struct RodGeomParams
{
    double fuel_outer_radius_ = 0.0;
    double fuel_inner_radius_ = 0.0;

    double contact_resistanñe_ = 0.0;

    double cladding_outer_radius_ = 0.0;
    double cladding_inner_radius_ = 0.0;
};


template<typename CladdingMaterial, typename FuelMaterial>
class Rod
{
protected:
	Cladding<CladdingMaterial> cladding;
	Gap gap;
	Fuel<FuelMaterial> fuel;

    Rod(RodGeomParams&& params_) 
        : cladding(params_.cladding_outer_radius_, params_.cladding_inner_radius_)
        , gap(params_.contact_resistanñe_)
        , fuel(params_.fuel_outer_radius_, params_.fuel_inner_radius_)
    {}
public:
    double perimeter()
    {
        return 2 * std::numbers::pi * cladding.get_outer_radius();
    }
};

#endif
