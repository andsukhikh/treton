#ifndef THECHYCO_HPP
#define THECHYCO_HPP
#include <vector>

double thehyco(double);
void HeatHydroOnce();
void CrossConnection();
void HeatHydroGeometry();
void normal();
double Sodium_KinVis(double, double);
double Sodium_Density(double, double);
double Lead_KinVis();
double Lead_Density(double, double );
void density();
void sy(std::vector<double>& , std::vector<double>& , std::vector<double>& , std::vector<double>& , int , int );
int is(int);
void InOut_f();
void write_all();
void read_all();
double V_r(int , int );
void Vxy_ij(int , int , double &, double &);
void VrFi_nm(int , int , double , double , double &, double &);
void V_full_calc();
double absV(int , int );

#endif
