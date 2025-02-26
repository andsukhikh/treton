#ifndef THECHYCO_HPP
#define THECHYCO_HPP
#include <vector>

#include "CoolantMaterials.hpp"


double thehyco(double dt, const Coolant& coolant);


void HeatHydroOnce(const Coolant& coolant);

void CrossConnection();
void HeatHydroGeometry();
void normal();


void density(const Coolant& coolant);

void sy(std::vector<double>& , std::vector<double>& , std::vector<double>& , std::vector<double>& , int , int );
int is(int);


void InOut_f(const Coolant& coolant);

double V_r(int , int );
void Vxy_ij(int , int , double &, double &);
void VrFi_nm(int , int , double , double , double &, double &);
void V_full_calc();
double absV(int , int );

#endif
