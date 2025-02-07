#ifndef HEAT_H
#define HEAT_H

void RodOnce();
void RodGeometry();
void rod_property();
double Rod(double dt, std::vector<double>& old, std::vector<double>& New, double alf, double tf, double Q_v, std::vector<double>& g_left, std::vector<double>& g_right, std::vector<double>& cor, double l, double rc);
void rod_average();
double EneRoDisbalance();
void heat(double dt);
double EnerFluiDisbalance();
void HeatConduction();
void alf();
double HeatTransfer(double Ux, double Uy, double Uz, double x, double d_hyd, double ent, double pvod);
double EnerCoreDisbalance();
void swapOnNewValue(std::vector<double>&, std::vector<std::vector<std::vector<std::vector<double>>>>&, int, int, int);

#endif 