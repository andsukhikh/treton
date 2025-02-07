#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "headers/SodiumProp.hpp"
#include "headers/ThechycoGlobalVar.hpp"
#include "headers/NamelistReader.hpp"
#include "headers/Hydro.hpp"
#include "headers/Heat.hpp"
#include "headers/CpSodium.hpp"
#include "headers/AMUV.hpp"
#include "headers/Thechyco.hpp"


double thehyco(double dt) {
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
        for (int j = 0; j < 50; ++j) {
            KinViscosity();
            //std::cout << "KinViscosity in thehyco is completely p = " << p[i][j] << std::endl;
            FormFriction();
            //std::cout << "FormFriction in thehyco is completely p = " << p[i][j] << std::endl;
            Viter(dt);
            //std::cout << "Viter in thehyco is completely p = "<< p[i][j] << std::endl;
            pes(dt);
            //std::cout << "pes in thehyco is completely p = " << p[i][j]<< std::endl;
            piter();
            //std::cout << "piter in thehyco is completely p[0][0] = "<< p[i][j] << std::endl;
            //std::cout << "V_z["<< i << "][" << j << "] = "<< V_z[i][j] << std::endl;

            MassDisb = MassDisbalance(dt);
            //std::cout << "MassDisbalance in thehyco is completely p = "<< p[i][j] << std::endl;

            if (MassDisb < Disbalance) {
                disbalanceSatisfied = true;
                break;
            }
        }

        if(!disbalanceSatisfied) {
            std::cout << "MassDisbalance = " << MassDisb << std::endl;
        }

        //std::cout << "iteration = " << iterations << std::endl;
        //std::cout << "i = " << i  << std::endl;
        //std::cout << "j = " << i << std::endl;

        alf();
        HeatConduction();
        heat(dt);
        density();
    }
    /*auto sodium_density = ro[n - 1][mf - 1];
    auto pressure = p[n - 1][mf - 1];*/

    double thehyco = 0.0;

    if (EneRoDisbalance() > Disbalance) thehyco = 1;
    if (EnerCoreDisbalance() > Disbalance) thehyco = 1;
    if (EnerFluiDisbalance() > Disbalance) thehyco = 1;

    //std::cout << "thehyco = " << static_cast<int>(thehyco) << std::endl;

    rod_average();
    //std::cout << "rod_average in thehyco is completely" << std::endl;
    rod_property();
    //std::cout << "rod_property in thehyco is completely" << std::endl;
    
    return thehyco;
}

//реализация с помощью namelist программа не умеет корректно считыать массив, значени которого раскиданы по строкам
//следует написать новую функцию, которая будет считыать одномерный и двумерный массив, раскиданный по строкам
void HeatHydroOnce() {
    static int manager = 0;

    // std::vector<std::string> varName = {"n_RodsInTBC", "crd", "D_tube", "dr", "dz", "p_input", "p_output", 
    //                                     "BlockadePorousity", "t_HeatExchangerOutput", "blockade", "PVTerror", 
    //                                     "Disbalance", "iterations"};
    // namelist HEATandHYDROlist("HEATandHYDROlist", varName, n_RodsInTBC, crd, D_tube, dr, dz, p_input, p_output, 
    //                                                 BlockadePorousity, t_HeatExchangerOutput, blockade, PVTerror, 
    //                                                 Disbalance, iterations);
    
    NLReader::NamelistReader nlr("THEHYCO.INI");

    // заполнение HEATandHYDROlist
    nlr.use_namelist("HEATandHYDROlist");
    dr = nlr.get<double>("dr", 1);
    dz = nlr.get<double>("dz", 1);
    D_tube = nlr.get<double>("D_tube", 1);
    t_HeatExchangerOutput = nlr.get<double>("t_HeatExchangerOutput", 1);
    Disbalance = nlr.get<double>("Disbalance", 1);
    PVTerror = nlr.get<double>("PVTerror", 1);
    BlockadePorousity = nlr.get<int>("BlockadePorousity", 1);
    p_input = nlr.get<double>("p_input", 1);
    p_output = nlr.get<double>("p_output", 1);
    iterations = nlr.get<int>("iterations", 1);
    
    for (int i = 0; i < type; ++i) {
        n_RodsInTBC[i] = nlr.get<int>("n_RodsInTBC", 1.0, i);
    }


    //заполнение blockade
    for (auto& var : blockade) {
        var = 0;
    }

    if (manager != 1) {
        // readNamelistFile(HEATandHYDROlist, "THEHYCO.INI");
        CrossConnection();
        HeatHydroGeometry();
        density();

        rod_property();
    
        manager = 1;
    }
}

             //WARNING поменял jx + 1 -> jx и ... = k++ на ... = ++k в bonds//****************//
void CrossConnection() {
    double n_x;
    double n_y;

    for (int i = 0; i < nbf; ++i) {
        for (int j = 0; j < mf; ++j) {
            bonds[i][j] = -1;
        }
    }


    for (int j = 0; j != mf; ++j) {
        //int i = -1;
        int i = 0;
        for (int jx = 0; jx != mf; ++jx) {
            //const auto crd_0jx = crd[0][jx];
            //const auto crd_0j = crd[0][j];
            //const auto crd_1jx = crd[1][jx];
            //const auto crd_1j = crd[1][j];

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

    int k = -1;// k = 0;
    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < nbf; ++i) {
            int jf = bonds[i][j];

            //std::cout << "bonds[" << i << "][" << j << "] = " << bonds[i][j] << std::endl;

            if (j < jf || jf == -1) {
                ++k;
                onds[i][j] = k;
                //std::cout << "onds[" << i << "][" << j << "] = " << onds[i][j] << std::endl;
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
                    exit(1);
                }
                //std::cout << "onds[" << i << "][" << j << "] = " << onds[i][j] << std::endl;
            }
            //std::cout << "bonds[" << i << "][" << j << "] = " << bonds[i][j] << std::endl;
        }
    }

    //const auto var = onds[nbf- 1][mf - 1];

    if (k + 1 != mV_n) {
        std::cout << "Для данного mf неправильно указано mV_n!" << std::endl;
        std::cout << "mV_n должно быть: " << k + 1 << std::endl;
        std::cout << "Error in BONDS; aborting..." << std::endl;
        exit(1);
    }

    // Заполнение NC и NE
    k = 0;
    int kk = 0;

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i == 0) {
                //NC[k] = 1;
                NC[k] = 0;
                NE[kk] = i + 1 + j * n; // +1 для индексации в C++
                ++kk;
            } else if (i == n - 1) {
                //NC[k] = 1;
                NC[k] = 0;
                NE[kk] = i - 1 + j * n; // индексация
                ++kk;
            } else {
                //NC[k] = 2;
                NC[k] = 1;
                NE[kk] = i - 1 + j * n; // индексация
                ++kk;
                NE[kk] = i + 1 + j * n; // индексация
                ++kk;
            }
            for (int l = 0; l < nbf; ++l) {
                if (bonds[l][j] != -1) {
                    ++NC[k];
                    NE[kk] = (bonds[l][j]) * n + i ; // индексация
                    ++kk;
                }
            }
            //std::cout << "NE[" << kk << "] = " << NE[kk] << std::endl;
            ++k;
        }
    }

    k = NC[0];
    NC[0] = 0;
    //k = NC[0] + 1;
    //NC[0] = 1;

    for (int i = 1; i < ( n * mf + 1 ); ++i) {
        int j = NC[i] + 1;
        NC[i] = k + 1;
        //std::cout << "NC[" << i << "] = " << NC[i] << std::endl;
        k += j;
    }

   /* for (int j = 0; j != mf; ++j) {
        for (int i = 0; i != nbf; ++i) {
            std::cout << "onds[" << i + 1 << "][" << j + 1 << "] = " << onds[i][j] << std::endl;
        }
    }*/
}


void HeatHydroGeometry() {
    fz = 0.5 * std::sqrt(3.0) * dr * dr;

    double v = fz * dz;

    double sum = 0.0;
   
    for (int k = 0; k < type; ++k) {
        sum += n_RodsInTBC[k] * (a_fuel[k] + a_clad[k]);
    }
    
    fz -= 0.25 * Pi * D_tube * D_tube - sum;
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


double Sodium_KinVis(double pvod, double ent) {
    double t_r = 0.0, v_r = 0.0;
    
    SodiumTV(pvod, ent, t_r, v_r);
    return AMUV(t_r) * v_r;
}



double Sodium_Density(double pvod, double ent) {

    double t_r = 0.0, v_r = 0.0;
    
    SodiumTV(pvod, ent, t_r, v_r);
    double ro = 1 / v_r;

    return ro;
}


void density() {
    InOut_f();

    double Cp_input = CpSodium(t_CoreInput);
    double ro_output = Sodium_Density(p_output, h_CoreOutput);
    double Cp_output = CpSodium(t_CoreOutput);

    for (int j = 0; j < mf; ++j) {
        ro_input[j] = Sodium_Density(p_input, h_HeatExchangerOutput_new[j]);
        for (int i = 0; i < n; ++i) {
            ro[i][j] = Sodium_Density(p[i][j], h_f[i][j]);
            C_p[i][j] = CpSodium(t_f[i][j]);
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
        //std::cout << "\tBB[" << i << "] = " << BB[i] << std::endl;
        //std::cout << "\tDD[" << i << "] = " << DD[i] << std::endl;
    }

    DD[IU] = DD[IU] / BB[IU];
    //std::cout << "\tDD[" << IU - 1 << "] = " << DD[IU - 1] << std::endl;
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


void InOut_f() {
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    double tmp3 = 0.0;
    double tmp4 = 0.0;
    double VAU = 0.0;
    
    for (int j = 0; j < mf; ++j) {
        double V_down = V_z[0][j];

        tmp1 += h_HeatExchangerOutput_new[j] *  Sodium_Density(p_input, h_HeatExchangerOutput_new[j]) * V_down;

        if (V_down < 0) {
            tmp1 -= h_f[0][j] * Sodium_Density(p_input, h_f[0][j]) * V_down;

        } else {
            tmp2 += Sodium_Density(p_input, h_HeatExchangerOutput_new[j]) * V_down;
        }

        double V_up = V_z[n][j];
        if (V_up > 0) {
            tmp3 += h_f[n - 1][j] * Sodium_Density(p_output, h_f[n - 1][j]) * V_up;

            tmp4 += Sodium_Density(p_output, h_f[n - 1][j]) * V_up;
            
        }
    }

    if (tmp2 != 0) {
        h_CoreInput = tmp1 / tmp2;
        SodiumTV(p_input, h_CoreInput, t_CoreInput, VAU);

    } else {
        std::cerr << "Global reverse flow on the bottom" << std::endl;
        exit(1);
    }

    if (tmp4 != 0) {
        h_CoreOutput = tmp3 / tmp4;
        SodiumTV(p_output, h_CoreOutput, t_CoreOutput, VAU);

    } else {
        std::cerr << "Global reverse flow on the top" << std::endl;
        exit(EXIT_FAILURE);
    }
}


void write_all() {
    std::ofstream file("THEHYCO.DAT", std::ios::binary);
    if (file.is_open()) {
        file.write(reinterpret_cast<char*>(&p), sizeof(p));
        file.write(reinterpret_cast<char*>(&V_z), sizeof(V_z));
        file.write(reinterpret_cast<char*>(&V_n), sizeof(V_n));
        file.write(reinterpret_cast<char*>(&h_f), sizeof(h_f));
        file.write(reinterpret_cast<char*>(&t_f), sizeof(t_f));
        file.write(reinterpret_cast<char*>(&t_rod), sizeof(t_rod));
        file.write(reinterpret_cast<char*>(&t_fuel), sizeof(t_fuel));
        file.write(reinterpret_cast<char*>(&t_clad), sizeof(t_clad));
        file.close();
    } else {
        std::cout << "Error in opening THEHYCO.DAT" << std::endl;
    }
}


void read_all() {
    std::ifstream file("THEHYCO.DAT", std::ios::binary);

    if (!file.eof()){
        file.read(reinterpret_cast<char*>(&p), sizeof(p));
        file.read(reinterpret_cast<char*>(&V_z), sizeof(V_z));
        file.read(reinterpret_cast<char*>(&V_n), sizeof(V_n));
        file.read(reinterpret_cast<char*>(&h_f), sizeof(h_f));
        file.read(reinterpret_cast<char*>(&t_f), sizeof(t_f));
        file.read(reinterpret_cast<char*>(&t_rod), sizeof(t_rod));
        file.read(reinterpret_cast<char*>(&t_fuel), sizeof(t_fuel));
        file.read(reinterpret_cast<char*>(&t_clad), sizeof(t_clad));
        file.close();
    } else { 
        std::cout << "THEHYCO.DAT is empty" << std::endl;
    }
}


double V_r(int i, int j) {
    double Vx, Vy;
    Vxy_ij(i, j, Vx, Vy);
    return sqrt(Vx * Vx + Vy * Vy);
}


void Vxy_ij(int i, int j, double &Vx_ij, double &Vy_ij) {
    const double sqrt3 = 1.7320508;
    int jf1 = bonds[0][j];
    int jf2 = bonds[1][j];
    int jf3 = bonds[2][j];
    int jf4 = bonds[3][j];                
    int jf5 = bonds[4][j];
    int jf6 = bonds[5][j];

    double Vn1 = 0.0, Vn2 = 0.0, Vn3 = 0.0, Vn4 = 0.0, Vn5 = 0.0, Vn6 = 0.0;

    if (jf1 != -1) Vn1 = is(jf1 - j) * V_n[i][bonds[0][j]];
    if (jf2 != -1) Vn2 = is(jf2 - j) * V_n[i][bonds[1][j]];
    if (jf3 != -1) Vn3 = is(jf3 - j) * V_n[i][bonds[2][j]];
    if (jf4 != -1) Vn4 = is(jf4 - j) * V_n[i][bonds[3][j]];
    if (jf5 != -1) Vn5 = is(jf5 - j) * V_n[i][bonds[4][j]];
    if (jf6 != -1) Vn6 = is(jf6 - j) * V_n[i][bonds[5][j]];

    double V1 = (Vn1 - Vn4);
    double V2 = (Vn2 - Vn5);
    double V3 = (Vn3 - Vn6);

    Vx_ij = (2.0 * V1 + V2 - V3) / 6.0;
    Vy_ij = (V2 + V3) / (2.0 * sqrt3);
}


void VrFi_nm(int i, int j, double Vn, double Vm, double &Vr, double &Fi) {
    Vr = std::sqrt(Vn * Vn + Vm * Vm);
    Fi = std::asin(Vm / Vr);

    //std::cout << "Vm = " << Vm << std::endl;
    //std::cout << "Vr = " << Vr << std::endl;

    
    if (Vn < 0.0) {
        Fi = Pi - Fi;
    } else if (Vm < 0.0) {
        Fi += 2.0 * Pi;
    }
}


void V_full_calc() {
    const double sqrt3 = 1.7320508;
    double Vx, Vy, Vx1, Vy1;

    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < mf; ++j) {
            switch (i) {
                case 0:
                    Vxy_ij(i, j, Vx, Vy);
                    //std::cout << "case 0: Vx = " << Vx << ", Vy = " << Vy << std::endl;
                    break;
                case n:
                    Vxy_ij(i - 1, j, Vx, Vy);
                    //std::cout << "case n: Vx = " << Vx << ", Vy = " << Vy << std::endl;

                    break;
                default:
                    Vxy_ij(i - 1, j, Vx1, Vy1);
                    Vxy_ij(i, j, Vx, Vy);
                    Vx = (Vx + Vx1) / 2.0;
                    Vy = (Vy + Vy1) / 2.0;
                    //std::cout << "case default: Vx = " << Vx << ", Vy = " << Vy << std::endl;
                    break;
            }

            V_full[0][i][j] = Vx;
            V_full[1][i][j] = Vy;
            V_full[2][i][j] = V_z[i][j];
        }
    }
}


double absV(int i, int j) {
    double result;
    
    switch (i) {
        case 0:
            result = std::sqrt(V_z[i][j] * V_z[i][j] + V_r(i, j) * V_r(i, j));
            break;
        case n:
            result = std::sqrt(V_z[i][j] * V_z[i][j] + V_r(i - 1, j) * V_r(i - 1, j));
            break;
        default:
            result = std::sqrt(V_z[i][j] * V_z[i][j] + 0.25 * std::pow(V_r(i - 1, j) + V_r(i, j), 2));
            break;
    }

    return result;
}
