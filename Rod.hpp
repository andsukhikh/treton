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
        : cladding      (params_.cladding_outer_radius_, params_.cladding_inner_radius_)
        , gap           (params_.contact_resistanñe_)
        , fuel          (params_.fuel_outer_radius_, params_.fuel_inner_radius_)
    {}
public:
    double perimeter()
    {
        return 2 * std::numbers::pi * cladding.get_outer_radius;
    }
};


//struct RodMeshParams
//{
//    double mash_step_;
//    double fuel_layers_number_;
//    double cladding_layers_number_;
//};


//template<typename CladdingMaterial, typename FuelMaterial>
//class RodMesh : protected Rod<CladdingMaterial, FuelMaterial>
//{
//protected:
//    double mash_step_;
//    double fuel_layers_number_;
//    double cladding_layers_number_;
//
//    std::vector<double> nodes_;
//    std::vector<double> fuel_layers_radius_;
//    std::vector<double> cladding_layers_radius_;
//
//    RodMesh(RodGeomParams&& geom_params, RodMeshParams&& mesh_params)
//        : Rod                       (std::move(geom_params))
//        , mash_step_                (mesh_params.mash_step_)
//        , fuel_layers_number_       (mesh_params.fuel_layers_number_)
//        , cladding_layers_number_   (mesh_params.cladding_layers_number_)
//
//        , nodes_                    (mesh_params.cladding_layers_number_ +  mesh_params.fuel_layers_number_ + 2)
//        , fuel_layers_radius_       (mesh_params.fuel_layers_number_ + 1)
//        , cladding_layers_radius_   (mesh_params.cladding_layers_number_ + 1)
//    {
//        decl_mesh_nodes();
//    }
//public:
//    void decl_fuel_mesh()
//    {
//        *fuel_layers_radius_.begin() = fuel.get_inner_radius();
//        *fuel_layers_radius_.rbegin() = fuel.get_outer_radius();
//
//        double first_rad = fuel.get_inner_radius();
//        double last_rad = fuel.get_outer_radius();
//
//        double eq_area = (std::pow(first_rad, 2) - std::pow(last_rad, 2)) / (2 * fuel_layers_number_);
//
//        for (size_t i = 1; i < fuel_layers_radius_.size(); ++i)
//        {
//            fuel_layers_radius_[i] = std::sqrt(std::pow(first_rad, 2) + (2 * i - 1) * eq_area);
//        }
//    }
//
//    void decl_cladding_mesh()
//    {
//        *cladding_layers_radius_.begin() = cladding.get_inner_radius();
//        *cladding_layers_radius_.rbegin() = cladding.get_outer_radius();
//
//        double first_rad = cladding.get_inner_radius();
//        double last_rad = cladding.get_outer_radius();
//
//        double eq_area = (std::pow(first_rad, 2) - std::pow(last_rad, 2)) / (2 * cladding_layers_number_);
//
//        for (size_t i = 1; i < cladding_layers_radius_.size(); ++i)
//        {
//            cladding_layers_radius_[i] = std::sqrt(std::pow(first_rad, 2) + (2 * i - 1) * eq_area);
//        }
//    }
//
//    void decl_mesh_nodes()
//    {
//        decl_fuel_mesh();
//        decl_cladding_mesh();
//
//        *nodes_.begin() = fuel_layers_radius_.begin();
//        for (size_t i = 1; i < fuel_layers_number_; ++i)
//        {
//            nodes_[i] = ( fuel_layers_radius_[i] + fuel_layers_radius_[i + 1] ) / 2;
//        }
//        nodes_[fuel_layers_number_] = fuel_layers_radius_.back();
//
//
//        for (size_t i = 1; i < cladding_layers_number_; ++i)
//        {
//            nodes_[fuel_layers_number_ + i] = cladding_layers_radius_[i];
//        }
//        nodes_[fuel_layers_number_ + cladding_layers_number_] = cladding_layers_radius_.back;
//    }
//};

#endif
