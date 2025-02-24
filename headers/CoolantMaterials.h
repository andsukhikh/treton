#ifndef COOLANT_MATERIALS_HPP
#define COOLANT_MATERIALS_HPP

#include <memory>
#include <string_view>

#include "headers/Coolant.h"

class Sodium : public Coolant {
public:
	double Entalpy (const double temperature)			const override;
	double Volume (const double temperature)			const override;
	double Temperature (const double enthalpy)			const override;
	double Pr (const double temperature)				const override;
	double density (const double temperature)			const override;
	double KinVis (const double entalpy)				const override;
	double DynVisc (const double temperature)			const override;
	double Cp (const double temperature)				const override;
	double HeatCond(const double temperature)			const override;
	double Nu(const double Re, const double Pr)			const override;
};

class Lead : public Coolant {
public:
	using type = double;
	double Entalpy (const double temperature)			const override;
	double Volume (const double temperature)			const override;
	double Temperature (const double enthalpy)			const override;
	double Pr (const double temperature)				const override;
	double density (const double temperature)			const override;
	double KinVis (const double entalpy)				const override;
	double DynVisc (const double temperature)			const override;
	double Cp (const double temperature)				const override;
	double HeatCond(const double temperature)			const override;
	double Nu(const double Re, const double Pr)			const override;
};

extern void define_coolant(std::unique_ptr<Coolant>& coolant, std::string_view coolantName);

#endif 
