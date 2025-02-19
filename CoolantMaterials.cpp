#include <string>
#include <cmath>
#include <memory>
#include <iostream>
#include <string_view>
#include <exception>

#include "headers/CoolantMaterials.h"

template class Sodium<double>;
template class Sodium<int>;

template class Lead<double>;
template class Lead<int>;

template void define_coolant(std::unique_ptr<Coolant<double>>& coolant, std::string_view coolantName);
template void define_coolant(std::unique_ptr<Coolant<int>>& coolant, std::string_view coolantName);

template<typename T>
void define_coolant(std::unique_ptr<Coolant<T>>& coolant, std::string_view coolantName) {
	if (coolantName == "Lead")					coolant = std::make_unique<Lead<T>>();
	else if (coolantName == "Sodium")			coolant = std::make_unique<Sodium<T>>();

	else { throw std::invalid_argument("This coolant is not in the database"); }
}


template<typename T>
T Sodium<T>::Entalpy(const T temperature) const {
	return -140954 + 1437 * temperature - 0.29 * std::pow(temperature, 2) + 1.54 * 1E-04 * std::pow(temperature, 3);
}

template<typename T>
T Sodium<T>::Volume(const T temperature) const {
	return 1 / density(temperature);
}

template<typename T>
T Sodium<T>::Temperature(const T enthalpy) const {
	return 0.00079 * enthalpy + 123.35;
}

template<typename T>
T Sodium<T>::Pr(const T temperature) const {
	T a = 1e-5 * (6.6951 + 5.264 * 1e-4 * temperature - 2.689 * 1e-6 * std::pow(temperature, 2));
	return DynVisc(temperature) / (density(temperature) * a);
}

template<typename T>
T Sodium<T>::density(const T temperature) const {
	return 949 - 0.223 * temperature - 1.75 * 1e-5 * std::pow(temperature, 2);
}

template<typename T>
T Sodium<T>::KinVis(const T entalpy) const {
	T temperature = Temperature(entalpy);
	T volume = Volume(temperature);

	return DynVisc(temperature) * volume;
}

template<typename T>
T Sodium<T>::DynVisc(const T temperature) const {
	return 1.2162 * 1e-5 * std::cbrt(density(temperature)) * exp(0.6976 * density(temperature) / (temperature + 273));
}

template<typename T>
T Sodium<T>::Cp(const T temperature) const {
	return 1436 - 0.5805 * temperature + 4.62 * 1e-4 * std::pow(temperature, 2);
}

template<typename T>
T Sodium<T>::HeatCond(const T temperature) const {
	return 90.6 - 0.04852 * temperature;
}

template<typename T>
T Sodium<T>::Nu(const T Re, const T Pr) const {
	if (Re <= 2300.0) {
		return 3.66;
	}
	else if (Re > 10000.0) {
		return 0.012 * std::pow(Pr, 0.33) * std::pow(Re, 0.87);
	}
	else {
		return 4.82 + 0.0185 * std::pow(Pr, 0.33) * std::pow(Re, 0.8);
	}
}




template<typename T>
T Lead<T>::Entalpy(const T temperature) const {
	double temperature_in_kelvin = temperature + 273;
	return 7.02 * (temperature_in_kelvin)-0.025 * std::pow(temperature_in_kelvin, 2) + 4e-4 * std::pow(temperature_in_kelvin, 3) - 2.95 * 1e-7 * std::pow(temperature_in_kelvin, 4);
}

template<typename T>
T Lead<T>::Volume(const T temperature) const {
	return 1 / density(temperature);
}

template<typename T>
T Lead<T>::Temperature(const T enthalpy) const {
	return 0.0037 * enthalpy + 214.81 - 273;
}

template<typename T>
T Lead<T>::Pr(const T temperature) const {
	return KinVis(Entalpy(temperature)) * Cp(temperature) / HeatCond(temperature);
}

template<typename T>
T Lead<T>::density(const T temperature) const {
	return (11.42 - 12.42 * 1e-4 * (temperature + 273)) * 1e3;
}

template<typename T>
T Lead<T>::KinVis(const T entalpy) const {
	T temperature = Temperature(entalpy);
	T volume = Volume(temperature);
	return DynVisc(temperature) * volume;
	//return std::pow((43.8 - 7.57 * 1e-2 * temperature + 0.467 * 1e-4 * temperature), 2) * 1e-8;
}

template<typename T>
T Lead<T>::DynVisc(const T temperature) const {
	return density(temperature) * std::pow( (43.8 - 7.57 * 1e-2 * temperature + 0.467 * 1e-4 * temperature), 2) * 1e-8;
}

template<typename T>
T Lead<T>::Cp(const T temperature) const {
	return 147.3;
}

template<typename T>
T Lead<T>::HeatCond(const T temperature) const {
	return 11 + 9.9 * 1e-3 * (temperature + 273 - 397);
}

template<typename T>
T Lead<T>::Nu(const T Re, const T Pr) const {
	if (Re <= 2300.0) {
		return 3.66;
	}
	else if (Re > 10000.0) {
		return 0.012 * std::pow(Pr, 0.33) * std::pow(Re, 0.87);
	}
	else {
		return 4.82 + 0.0185 * std::pow(Pr, 0.33) * std::pow(Re, 0.8);
	}
}

