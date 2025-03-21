#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numbers>
#include <type_traits>

#include "GlobalVar.hpp"
#include "NamelistReader.hpp"
#include "Hydro.hpp"
#include "Heat.hpp"
#include "Thechyco.hpp"
#include "CoolantMaterials.hpp"


double thehyco(double dt, const Coolant& coolant) {
    double MassDisb;

    for (int j = 0; j != mf; ++j) {
        for (int i = 0; i != n; ++i) {
            OLDh_f[i][j] = h_f[i][j];
            OLDt_f[i][j] = t_f[i][j];
            OLD_ro[i][j] = ro[i][j];
            OLD_C_p[i][j] = C_p[i][j];
            for (int k = 0; k != type; ++k) {
                for (int l = 0; l !=  n_rod + 2 ; ++l) {
                    OLDt_rod[l][i][j][k] = t_rod[l][i][j][k];
                }
            }
        }
    }


    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n + 1; ++i) {
            oldV_z[i][j] = V_z[i][j];
        }
    }


    for (int j = 0; j < mV_n; ++j) {
        for (int i = 0; i != n; ++i) {
            oldV_n[i][j] = V_n[i][j];
        }
    }


    for (int i = 0; i < iterations; ++i) {
        bool disbalanceSatisfied = false;
        for (int j = 0; j < 500'000; ++j) {
            KinViscosity(coolant);
            FormFriction(coolant);
            Viter(dt);
            pes(dt);
            piter();
            MassDisb = MassDisbalance(dt);

            if (MassDisb < Disbalance) {
                disbalanceSatisfied = true;
                break;
            }
        }

        if(!disbalanceSatisfied) {
            std::cout << "MassDisbalance = " << MassDisb << std::endl;
        }

        alf(coolant);
        HeatConduction(coolant);
        heat(dt, coolant);
        density(coolant);
    }

    double thehyco = 0.0;

    if (EneRoDisbalance() > Disbalance) thehyco = 1;
    if (EnerCoreDisbalance() > Disbalance) thehyco = 1;
    if (EnerFluiDisbalance() > Disbalance) thehyco = 1;

    rod_average();
    rod_property();
    
    return thehyco;
}


void HeatHydroOnce(const Coolant& coolant) {
    static int manager = 0;

    if (manager != 1) {
        CrossConnection();
        HeatHydroGeometry();
        density(coolant);

        rod_property();
    
        manager = 1;
    }
}

void CrossConnection() {
    double n_x;
    double n_y;

    for (int i = 0; i < nbf; ++i) {
        for (int j = 0; j < mf; ++j) {
            bonds[i][j] = -1;
        }
    }


    for (int j = 0; j != mf; ++j) {
        int i = 0;
        for (int jx = 0; jx != mf; ++jx) {

            n_x = crd[0][jx] - crd[0][j];
            n_y = crd[1][jx] - crd[1][j];

            if (n_x == 2 && n_y == 0) {
                bonds[0][j] = jx; 
                ++i;
            } else if (n_x == 1 && n_y == 1) {
                bonds[1][j] = jx;
                ++i;
            } else if (n_x == -1 && n_y == 1) {
                bonds[2][j] = jx;
                ++i;
            } else if (n_x == -2 && n_y == 0) {
                bonds[3][j] = jx;
                ++i;
            } else if (n_x == -1 && n_y == -1) {
                bonds[4][j] = jx;
                ++i;
            } else if (n_x == 1 && n_y == -1) {
                bonds[5][j] = jx;
                ++i;
            }

            if (i == nbf) break; 
        }
    }

    for (int i = 0; i != nbf; ++i) {
        for (int j = 0; j != mf; ++j) {
            onds[i][j] = -1;
        }
    }

    int k = -1;
    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < nbf; ++i) {
            int jf = bonds[i][j];
            
            if (j < jf || jf == -1) {
                ++k;
                onds[i][j] = k;
            } else {

                int i_f;
                int jx;
                for (i_f = 0; i_f < nbf; ++i_f) {
                    jx = bonds[i_f][jf];
                    if (jx == j) break;
                }

                if (jx == j) {
                    onds[i][j] = onds[i_f][jf];
                } else {
                    std::cout << "Error in BONDS; aborting..." << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        }
    }

    if (k + 1 != mV_n) {
        mV_n = k + 1;
        GlobalVariables::Resizing();
    }

    // Заполнение 
    // и NE
    k = 0;
    int kk = 0;

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            // нумерация по высоте
            if (i == 0) {
                NC[k] = 1;
                NE[kk] = i + 1 + j * n; 
                ++kk;
            } else if (i == n - 1) {
                NC[k] = 1;
                NE[kk] = i - 1 + j * n;
                ++kk;
            } else {
                NC[k] = 2;
                NE[kk] = i - 1 + j * n; 
                ++kk;
                NE[kk] = i + 1 + j * n; 
                ++kk;
            } 

            //нумерация в поперечном сечении
            for (int l = 0; l < nbf; ++l) {
                if (bonds[l][j] != -1) {
                    ++NC[k];
                    NE[kk] = (bonds[l][j]) * n + i ;
                    ++kk;
                }
            }
            ++k;
        }
    }

    k = NC[0];
    NC[0] = 0;

    for (int i = 1; i < ( n * mf + 1 ); ++i) {
        int j = NC[i];
        NC[i] = k;
        k += j;
    }
}


void HeatHydroGeometry() {
    fz = 0.5 * std::sqrt(3.0) * dr * dr;
    double v = fz * dz;

    double sum = 0.0;
   
    for (int k = 0; k < type; ++k) {
        sum += n_RodsInTBC[k] * (a_fuel[k] + a_clad[k]);
    }
    
    fz -= 0.25 * std::numbers::pi * D_tube * D_tube + sum;
    vf = fz * dz;
    fr = vf / v * dz * dr / std::sqrt(3.0);
    fr_vf = fr / vf;
    fz_vf = fz / vf;
}


void normal(int n, double &x, double &y) {
    double vns[2][6] = {
        {1.0, 0.0, 0.5, 0.8660254, -0.5, 0.8660254},
        {-1.0, 0.0, -0.5, -0.8660254, 0.5, -0.8660254}
    };
    
    x = vns[0][n - 1];
    y = vns[1][n - 1];
}


void density(const Coolant& coolant) {
    InOut_f(coolant);

    Cp_input = coolant.Cp(t_CoreInput);
    ro_output = coolant.density(coolant.Temperature(h_CoreOutput));
    Cp_output = coolant.Cp(t_CoreOutput);

    for (int j = 0; j < mf; ++j) {
        ro_input[j] = coolant.density(coolant.Temperature(h_HeatExchangerOutput_new[j]));
        for (int i = 0; i < n; ++i) {
            ro[i][j] = coolant.density(coolant.Temperature(h_f[i][j]));
            C_p[i][j] = coolant.Cp(t_f[i][j]);
        }
    }
}


void sy(std::vector<double>& AA, std::vector<double>& BB, std::vector<double>& CC, std::vector<double>& DD, int IL, int IU) {
    //SUBROUTINE SY SOLVES TRIDIAGONAL SYSTEM BY ELIMINA
    //IL=SUBSCRIPT OF FIRST EQUATION
    //IU=SUBSCRIPT OF LAST EQUATION
    //AA=COEFFICIENT BEHIND  DIAGONAL
    //BB=COEFFICIENT ON  DIAGONAL
    //CC=COEFFICIENT AHEAD  DIAGONAL
    //DD=ELEMENT OF CONSTANT VECTOR

    //ESTABLISH UPPER TRIANGULAR MATRIX
    for (int i = IL + 1; i <= IU; ++i) {
        double R = AA[i] / BB[i - 1];
        BB[i] -= R * CC[i - 1];
        DD[i] -= R * DD[i - 1];
    }

    DD[IU] = DD[IU] / BB[IU];
    for (int i = IU - 1; i >= IL; --i) {
        DD[i] = (DD[i] - CC[i] * DD[i + 1]) / BB[i];
    }
}


int is(int iarg) {
    if (iarg > 0) {
        return +1;
    } else {
        return -1;
    }
}


void InOut_f(const Coolant& coolant) {
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    double tmp3 = 0.0;
    double tmp4 = 0.0;
    double VAU = 0.0;
    
    for (int j = 0; j < mf; ++j) {
        double V_down = V_z[0][j];

        tmp1 += h_HeatExchangerOutput_new[j] * coolant.density(coolant.Temperature(h_HeatExchangerOutput_new[j])) * V_down;

        if (V_down < 0) {
            tmp1 -= h_f[0][j] * coolant.density(coolant.Temperature(h_f[0][j])) * V_down;

        } else {
            tmp2 += coolant.density(coolant.Temperature(h_HeatExchangerOutput_new[j])) * V_down;

        }

        double V_up = V_z[n][j];
        if (V_up > 0) {
            tmp3 += h_f[n - 1][j] * coolant.density(coolant.Temperature(h_f[n - 1][j])) * V_up;

            tmp4 += coolant.density(coolant.Temperature(h_f[n - 1][j])) * V_up;

        }
    }

    if (tmp2 != 0) {
        h_CoreInput = tmp1 / tmp2;
        t_CoreInput = coolant.Temperature(h_CoreInput);

    } else {
        std::cerr << "Global reverse flow on the bottom" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (tmp4 != 0) {
        h_CoreOutput = tmp3 / tmp4;
        t_CoreOutput = coolant.Temperature(h_CoreOutput);

    } else {
        std::cerr << "Global reverse flow on the top" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


double V_r(int i, int j) {
    double Vx, Vy;
    Vxy_ij(i, j, Vx, Vy);
    return std::sqrt(Vx * Vx + Vy * Vy);
}


void Vxy_ij(int i, int j, double &Vx_ij, double &Vy_ij) {
    int jf1 = bonds[0][j];
    int jf2 = bonds[1][j];
    int jf3 = bonds[2][j];
    int jf4 = bonds[3][j];                
    int jf5 = bonds[4][j];
    int jf6 = bonds[5][j];

    double Vn1 = 0.0, Vn2 = 0.0, Vn3 = 0.0, Vn4 = 0.0, Vn5 = 0.0, Vn6 = 0.0;

    if (jf1 != -1) Vn1 = is(jf1 - j) * V_n[i][onds[0][j]];
    if (jf2 != -1) Vn2 = is(jf2 - j) * V_n[i][onds[1][j]];
    if (jf3 != -1) Vn3 = is(jf3 - j) * V_n[i][onds[2][j]];
    if (jf4 != -1) Vn4 = is(jf4 - j) * V_n[i][onds[3][j]];
    if (jf5 != -1) Vn5 = is(jf5 - j) * V_n[i][onds[4][j]];
    if (jf6 != -1) Vn6 = is(jf6 - j) * V_n[i][onds[5][j]];

    double V1 = (Vn1 - Vn4);
    double V2 = (Vn2 - Vn5);
    double V3 = (Vn3 - Vn6);

    Vx_ij = (2.0 * V1 + V2 - V3) / 6.0;
    Vy_ij = (V2 + V3) / 2.0 / std::sqrt(3);
}


void VrFi_nm(int i, int j, double Vn, double Vm, double &Vr, double &Fi) {
    Vr = std::sqrt(Vn * Vn + Vm * Vm);
    Fi = std::asin(Vm / Vr);
    
    if (Vn < 0.0) {
        Fi = std::numbers::pi - Fi;
    } else if (Vm < 0.0) {
        Fi += 2.0 * std::numbers::pi;
    }
}


void V_full_calc() {
    double Vx, Vy, Vx1, Vy1;

    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < mf; ++j) {

            if (i == 0) {
                Vxy_ij(i, j, Vx, Vy);
            }
            else if (i == n) {
                Vxy_ij(i - 1, j, Vx, Vy);
            }
            else {
                Vxy_ij(i - 1, j, Vx1, Vy1);
                Vxy_ij(i, j, Vx, Vy);
                Vx = (Vx + Vx1) / 2.0;
                Vy = (Vy + Vy1) / 2.0;
            }

            V_full[0][i][j] = Vx;
            V_full[1][i][j] = Vy;
            V_full[2][i][j] = V_z[i][j];
        }
    }
}


double absV(int i, int j) {

    if (i == 0) {
        return std::sqrt(V_z[i][j] * V_z[i][j] + V_r(i, j) * V_r(i, j));
    }
    if (i == n) {
        return std::sqrt(V_z[i][j] * V_z[i][j] + V_r(i - 1, j) * V_r(i - 1, j));
    }
    return std::sqrt(V_z[i][j] * V_z[i][j] + 0.25 * std::pow(V_r(i - 1, j) + V_r(i, j), 2));
}
