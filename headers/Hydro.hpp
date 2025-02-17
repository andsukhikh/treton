#ifndef HYDRO_H
#define HYDRO_H
#include "CoolantMaterials.h"

double V_nz(int i, int j, int jf, int jV_n);
double MassDisbalance(double dt);
double Groz(int i, int j);
void V_zBlockade();
double random_number();
void Viter(double dt);
void pes(double dt);
void piter();

template<typename T>
void KinViscosity(const Coolant<T>& coolant);

template<typename T>
void FormFriction(const Coolant<T>& coolant);

#endif