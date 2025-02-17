#ifndef THECHYCO_HPP
#define THECHYCO_HPP
#include <vector>

#include "CoolantMaterials.h"

template<typename T>
T thehyco(T dt, const Coolant<T>& coolant);

template<typename T>
void HeatHydroOnce(const Coolant<T>& coolant);

void CrossConnection();
void HeatHydroGeometry();
void normal();
double Sodium_KinVis(double, double);
double Sodium_Density(double, double);
double Lead_KinVis();
double Lead_Density(double, double );

template<typename T>
void density(const Coolant<T>& coolant);

void sy(std::vector<double>& , std::vector<double>& , std::vector<double>& , std::vector<double>& , int , int );
int is(int);

template<typename T>
void InOut_f(const Coolant<T>& coolant);

void write_all();
void read_all();
double V_r(int , int );
void Vxy_ij(int , int , double &, double &);
void VrFi_nm(int , int , double , double , double &, double &);
void V_full_calc();
double absV(int , int );

#endif
