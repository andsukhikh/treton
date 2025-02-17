#ifndef COOLANT_HPP
#define COOLANT_HPP

template<typename T>
class Coolant {
public:
	using type = T;

	virtual T Entalpy(const T temperature) const = 0;
	virtual T Volume(const T temperature) const = 0;
	virtual T Temperature(const T enthalpy) const = 0;
	virtual T Pr(const T temperature) const = 0;
	virtual T KinVis(const T enthalpy) const = 0;
	virtual T DynVisc(const T temperature) const = 0;
	virtual T density(const T temperature) const = 0;
	virtual T Cp (const T temperature) const = 0;
	virtual T HeatCond(const T temperature) const = 0;
	virtual T Nu (const T Re, const T Pr) const = 0;

	virtual ~Coolant() = default;
};

#endif
