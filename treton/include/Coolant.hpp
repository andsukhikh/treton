#ifndef COOLANT_HPP
#define COOLANT_HPP

class Coolant {
public:
	virtual double Entalpy(const double temperature)		const = 0;
	virtual double Volume(const double temperature)			const = 0;
	virtual double Temperature(const double enthalpy)		const = 0;
	virtual double Pr(const double temperature)				const = 0;
	virtual double KinVis(const double enthalpy)			const = 0;
	virtual double DynVisc(const double temperature)		const = 0;
	virtual double density(const double temperature)		const = 0;
	virtual double Cp (const double temperature)			const = 0;
	virtual double HeatCond(const double temperature)		const = 0;
	virtual double Nu (const double Re, const double Pr)	const = 0;

	virtual ~Coolant() = default;
};

#endif
