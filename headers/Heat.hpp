#ifndef HEAT_H
#define HEAT_H

#include "CoolantMaterials.h"

void RodOnce();
void RodGeometry();
void rod_property();
double Rod(double dt, std::vector<double>& old, std::vector<double>& New, double alf, double tf, double Q_v, std::vector<double>& g_left, std::vector<double>& g_right, std::vector<double>& cor, double l, double rc);
void rod_average();
double EneRoDisbalance();

template<typename T>
void heat(T dt, const Coolant<T>& coolant);

double EnerFluiDisbalance();

template<typename T>
void HeatConduction(const Coolant<T>& coolant);

template<typename T>
void alf(const Coolant<T>& coolant);

template<typename T>
T HeatTransfer(T Ux, T Uy, T Uz, T x, T d_hyd, T ent, const Coolant<T>& coolant);


double EnerCoreDisbalance();
void swapOnNewValue(std::vector<double>&, std::vector<std::vector<std::vector<std::vector<double>>>>&, int, int, int);

#endif 