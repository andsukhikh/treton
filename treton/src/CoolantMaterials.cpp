#include <string>
#include <cmath>
#include <memory>
#include <iostream>
#include <string_view>
#include <exception>
#include <cstdlib>

#include "CoolantMaterials.hpp"


CoolantDecriptor::CoolantDecriptor(std::string_view coolantName) {
	if (coolantName == "Lead") {
		coolant_ = std::make_unique<Lead>();

	} else if (coolantName == "Sodium") {
		coolant_ = std::make_unique<Sodium>();

	} else {
		throw std::invalid_argument("This coolant is not in the database");
	}
}

Coolant& CoolantDecriptor::getCoolant() {
	return *coolant_;
}


double Sodium::Entalpy(const double temperature) const {
	return -140954 + 1437 * temperature - 0.29 * std::pow(temperature, 2) + 1.54 * 1E-04 * std::pow(temperature, 3);
}

double Sodium::Volume(const double temperature) const {
	return 1 / density(temperature);
}

double Sodium::Temperature(const double enthalpy) const {
	return 85.033 + 0.0008 * enthalpy - 2e-11 * std::pow(enthalpy, 2);
}

double Sodium::Pr(const double temperature) const {
	double a = 1e-5 * (6.6951 + 5.264 * 1e-4 * temperature - 2.689 * 1e-6 * std::pow(temperature, 2));
	return DynVisc(temperature) / (density(temperature) * a);
}

double Sodium::density(const double temperature) const {
	return 949 - 0.223 * temperature - 1.75 * 1e-5 * std::pow(temperature, 2);
}

double Sodium::KinVis(const double entalpy) const {
	double temperature = Temperature(entalpy);
	double volume = Volume(temperature);

	return DynVisc(temperature) * volume;
}

double Sodium::DynVisc(const double temperature) const {
	return 1.2162 * 1e-5 * std::cbrt(density(temperature)) * exp(0.6976 * density(temperature) / (temperature + 273));
}

double Sodium::Cp(const double temperature) const {
	return 1436 - 0.5805 * temperature + 4.62 * 1e-4 * std::pow(temperature, 2);
}

double Sodium::HeatCond(const double temperature) const {
	return 90.6 - 0.04852 * temperature;
}

double Sodium::Nu(const double Re, const double Pr) const {
	double Pe = Re * Pr;

	if (Re <= 2300.0) {
		return 4.36;
	}
	if (Re > 10000.0) {
		return 5 + 0.025 * std::pow(Pe, 0.8);
	}
	return 5 + 0.025 * std::pow(Pe, 0.8);
}




double Lead::Entalpy(const double temperature) const {
	double temperature_in_kelvin = temperature + 273;
	return 7.02 * (temperature_in_kelvin)-0.025 * std::pow(temperature_in_kelvin, 2) + 4e-4 * std::pow(temperature_in_kelvin, 3) - 2.95 * 1e-7 * std::pow(temperature_in_kelvin, 4);
}

double Lead::Volume(const double temperature) const {
	return 1 / density(temperature);
}

double Lead::Temperature(const double enthalpy) const {
	return 0.0037 * enthalpy + 214.81 - 273;
}

double Lead::Pr(const double temperature) const {
	return KinVis(Entalpy(temperature)) * Cp(temperature) / HeatCond(temperature);
}

double Lead::density(const double temperature) const {
	return (11.42 - 12.42 * 1e-4 * (temperature + 273)) * 1e3;
}

double Lead::KinVis(const double entalpy) const {
	double temperature = Temperature(entalpy);
	double volume = Volume(temperature);
	return DynVisc(temperature) * volume;
}

double Lead::DynVisc(const double temperature) const {
	return density(temperature) * std::pow( (43.8 - 7.57 * 1e-2 * temperature + 0.467 * 1e-4 * temperature), 2) * 1e-8;
}

double Lead::Cp(const double temperature) const {
	return 147.3;
}

double Lead::HeatCond(const double temperature) const {
	return 11 + 9.9 * 1e-3 * (temperature + 273 - 397);
}

double Lead::Nu(const double Re, const double Pr) const {
	double Pe = Re * Pr;

	if (Re <= 2300.0) {
		return 4.36;
	}
	if (Re > 10000.0) {
		return 3 + 0.014 * std::pow(Pe, 0.8);
	}
	return 3 + 0.014 * std::pow(Pe, 0.8);
}

