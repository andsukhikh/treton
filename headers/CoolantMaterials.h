#ifndef COOLANT_MATERIALS_HPP
#define COOLANT_MATERIALS_HPP

#include <memory>
#include "Coolant.h"


template<typename T>
extern void define_coolant(std::unique_ptr<Coolant<T>>& coolant, const std::string coolantName);


template<typename T>
class Sodium : public Coolant<T> {
public:
	using type = T;
	T Entalpy (const T temperature)			const override;
	T Volume (const T temperature)			const override;
	T Temperature (const T enthalpy)		const override;
	T Pr (const T temperature)					const override;
	T density (const T temperature)				const override;
	T KinVis (const T entalpy)					const override;
	T DynVisc (const T temperature)				const override;
	T Cp (const T temperature)					const override;
	T HeatCond(const T temperature)				const override;
	T Nu(const T Re, const T Pr)						const override;
};

template<typename T>
class Lead : public Coolant<T> {
public:
	using type = T;
	T Entalpy (const T temperature)			const override;
	T Volume (const T temperature)			const override;
	T Temperature (const T enthalpy)		const override;
	T Pr (const T temperature)					const override;
	T density (const T temperature)				const override;
	T KinVis (const T entalpy)					const override;
	T DynVisc (const T temperature)				const override;
	T Cp (const T temperature)					const override;
	T HeatCond(const T temperature)				const override;
	T Nu(const T Re, const T Pr)						const override;
};


#endif 
