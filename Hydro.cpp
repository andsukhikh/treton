#include <random>
#include <numbers>
#include <cmath>
#include <numbers>
#include <ctime>
#include <iostream>

#include "headers/Thechyco.hpp"
#include "headers/ThechycoGlobalVar.hpp"
#include "headers/CoolantMaterials.h"
#include "headers/Hydro.hpp"


double V_nz(int i, int j, int jf, int jV_n) {
    double avgV_z = (V_z[i][j] + V_z[i + 1][j] + V_z[i][jf] + V_z[i + 1][jf]) / 4.0;
    double result = std::sqrt(avgV_z * avgV_z + (std::numbers::pi / 2.0 * V_n[i][jV_n]) * (std::numbers::pi / 2.0 * V_n[i][jV_n]));
    return result;
}


double MassDisbalance(double dt) {
    
    double MassDisbalance = 0.0;
    
    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            aPV[i] = 0.0;
            bPV[i] = 0.0;

            double dro_dt = vf * (ro[i][j] - OLD_ro[i][j]) / dt;
            if (dro_dt > 0) {
                bPV[i] += dro_dt;
            } else {
                aPV[i] += dro_dt;
            }

            double V_up = fz * V_z[i + 1][j];
            double V_down = -fz * V_z[i][j];
            double ro_up, ro_down;

            if (i == 0) {
                ro_up = (ro[i][j] + ro[i + 1][j]) / 2.0;
                ro_down = ro[i][j];
            } else if (i == n - 1) {
                ro_up = ro[i][j];
                ro_down = (ro[i][j] + ro[i - 1][j]) / 2.0;
            } else {
                ro_up = (ro[i][j] + ro[i + 1][j]) / 2.0;
                ro_down = (ro[i][j] + ro[i - 1][j]) / 2.0;
            }

            if (V_up > 0) {
                bPV[i] += V_up * ro_up;
            } else {
                aPV[i] += V_up * ro_up;
            }

            if (V_down > 0) {
                bPV[i] += V_down * ro_down;
            } else {
                aPV[i] += V_down * ro_down;
            }
        }

        for (int k = 0; k < nbf; ++k) {
            int jf = bonds[k][j];
            if (jf != -1) {
                int jV_n = onds[k][j];
                for (int i = 0; i < n; ++i) {
                    double Vn = is(jf - j) * V_n[i][jV_n] * fr;
                    double roV_n = (ro[i][j] + ro[i][jf]) / 2.0;

                    if (Vn > 0) {
                        bPV[i] += Vn * roV_n;
                    } else {
                        aPV[i] += Vn * roV_n;
                    }
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            if (bPV[i] != 0.0) {
                double error = std::abs(1.0 + aPV[i] / bPV[i]);
                if (MassDisbalance < error) {
                    MassDisbalance = error;
                    int i_error = i;
                    int j_error = j;
                }
            }
        }
    }
    
    return MassDisbalance;
}



double Groz(int i, int j) {
    double Groz = 0.0;
    int iz = i;
    
    while (iz < n) {
        Groz += ro[iz][j];
        iz++;
    }

    Groz = Groz * Gravity * dz;
    return Groz;
}


void V_zBlockade() {
    for (int j = 0; j < mf; j++) { 
        if (blockade[j] == 1) {
            V_z[0][j] *= BlockadePorousity; 
        }
    }
}


double random_number() {

    static std::default_random_engine generator(static_cast<unsigned int>(time(0))); 
    std::uniform_real_distribution<double> distribution(0.0, 1.0); 

    return distribution(generator); 
}



void Viter(double dt) 
{
    std::vector<int> targ(mf, 0);
    double x;
    int sum;
    double error = 1.0;
    int iii = 0;
    double V_up, V_down, jf, jV_n, const_Vn, cr_M, roV_n, Vn, e;   

    while ((error > PVTerror) && (iii < 1500)) 
    {
        error = 0.0;
        iii++;  
                    /**************V_z ******************/
        std::fill(targ.begin(), targ.end(), 0);

        sum = 0;
        while (sum < mf) 
        {
            x = random_number();
            int j = static_cast<int>(mf * x);

            if (targ[j] == 1) {
                continue;
            } else {
                targ[j] = 1;
                sum++;

                for (int i = 0; i < n + 1; ++i)
                {
                    aPV[i] = 0.0;
                    cPV[i] = 0.0;
                    ePV[i] = V_z[i][j];
                    fPV[i] = std::abs(V_z[i][j]);
                }

                int i = 0;

                switch(blockade[j])
                {
                    case 1: 
                        bPV[i] = 1.0;
                        dPV[i] = V_z[i][j];
                    break;

                    case 0:
                        bPV[i] = ro[i][j] * (1.0 / dt + effK_z[i][j] * fPV[i]);
                        dPV[i] = OLD_ro[i][j] * oldV_z[i][j] / dt - 2.0 * (p[i][j] - p_input) / dz - ro[i][j] * Gravity;
                        double V_up = fz_vf * (V_z[i][j] + V_z[i + 1][j]);
                        double V_down = -fz_vf * V_z[i][j] * 2.0;
                        
                        if(V_up > 0.0) {
                            bPV[i] += V_up * ro[i][j];
                        } else {
                            cPV[i] += V_up * ro[i][j];
                        }

                        if(V_down > 0.0) {
                            bPV[i] += V_down * ro[i][j];
                        } else {
                            bPV[i] += V_down * ro_input[j];
                        }

                        for (int k = 0; k < nbf; ++k) 
                        {
                            int jf = bonds[k][j];
                            if (jf != -1)
                            {
                                int jV_n = onds[k][j];
                                double Vn = fr_vf * is(jf - j) * V_n[i][jV_n];
                                double cr_M = fr_vf * 0.5 * (effM[i][j] + effM[i][jf]) / dr;
                                bPV[i] += cr_M;
                                dPV[i] += cr_M * V_z[i][jf];

                                if (Vn > 0.0) {
                                    bPV[i] += Vn * ro[i][j];
                                } else {
                                    dPV[i] -= Vn * ro[i][jf] * V_z[i][jf];
                                }
                            }
                        }
                    break;
                }


                for(int i = 1; i < n; ++i) 
                {
                    bPV[i] = 0.5 * (ro[i][j] + ro[i - 1][j]) * (1.0 / dt + effK_z[i][j] * fPV[i]);
                    dPV[i] = 0.5 * (OLD_ro[i][j] + OLD_ro[i - 1][j]) / dt * oldV_z[i][j] 
                            -(p[i][j] - p[i - 1][j]) / dz - 0.5 * (ro[i][j] + ro[i - 1][j]) * Gravity;
                    double V_up = 0.5 * fz_vf * (V_z[i][j] + V_z[i + 1][j]);  
                    double V_down = -0.5 * fz_vf * (V_z[i][j] + V_z[i - 1][j]);

                    if(V_up > 0.0) {
                        bPV[i] += V_up * ro[i][j];
                    } else {
                        cPV[i] += V_up * ro[i][j];
                    }

                    if(V_down > 0.0) {
                        bPV[i] += V_down * ro[i - 1][j];
                    } else {
                        aPV[i] += V_down * ro[i - 1][j];
                    }
                }

                i = n;
                
                bPV[i] = ro[i - 1][j] * (1.0 / dt + effK_z[i][j] * fPV[i]);
                dPV[i] = OLD_ro[i - 1][j] * oldV_z[i][j] / dt - 2.0 * (p_output - p[i - 1][j]) / dz 
                        - ro[i - 1][j] * Gravity;
                double V_up = 2.0 * fz_vf * V_z[i][j];
                double V_down = -fz_vf * (V_z[i][j] + V_z[i - 1][j]);

                if (V_up > 0.0) {
                    bPV[i] += V_up * ro[i - 1][j];
                } else {
                    bPV[i] += V_up * ro_output;
                }

                if(V_down > 0.0) {
                    bPV[i] += V_down * ro[i - 1][j];
                } else {
                    aPV[i] += V_down * ro[i - 1][j];
                }

                for(int k = 0; k < nbf; k++) 
                {
                    double jf = bonds[k][j];
                    if (jf != -1){
                        jV_n = onds[k][j];
                        const_Vn = 0.5 * fr_vf * is(jf - j);

                        for (int i = 1; i < n; ++i) {
                            Vn = const_Vn * (V_n[i - 1][jV_n] + V_n[i][jV_n]);
                            cr_M = fr_vf * (effM[i][j] + effM[i][jf]) / 2.0 / dr;
                            bPV[i] += cr_M;
                            dPV[i] += cr_M * V_z[i][jf];

                            if (Vn > 0.0) {
                                bPV[i] += Vn * (ro[i - 1][j] + ro[i][j]) / 2.0;
                            } else {
                                dPV[i] -= Vn * (ro[i - 1][jf] + ro[i][j]) / 2 * V_z[i][jf];
                            }
                        }

                        i = n;

                        Vn = fr_vf * is(jf - j) * V_n[i - 1][jV_n];
                        cr_M = 0.5 * fr_vf * (effM[i][j] + effM[i][jf]) / dr;
                        bPV[i] += cr_M;
                        dPV[i] += cr_M * V_z[i][jf];

                        if (Vn > 0.0) {
                            bPV[i] += Vn * ro[i - 1][j];
                        } else {
                            dPV[i] -= Vn * ro[i - 1][jf] * V_z[i][jf];
                        }

                    }
                }
                sy(aPV, bPV, cPV, dPV, 0, n);
                for(int i = 0; i < n + 1; i++) 
                {
                    V_z[i][j] = dPV[i];
                    if(ePV[i] != 0.0)
                    {
                        e = abs(1.0 - V_z[i][j] / ePV[i]);
                        if(e > error) error = e;
                    }
                }
            }
        } 

                    /*****V_n****** */

        std::fill(V_nMap.begin(), V_nMap.end(), 0);

        std::fill(targ.begin(), targ.end(), 0);

        sum = 0;
        while (sum < mf) {
            x = random_number();
            int j = static_cast<int>(mf * x);
            
            if (targ[j] == 1) {
                continue;
            } else {
                targ[j] = 1;
                sum++;

                for (int k = 0; k < nbf; ++k) {
                    int jf = bonds[k][j];
                    int jV_n = onds[k][j];
                    
                    if (V_nMap[jV_n] == 0) {
                        if(jf == -1) {
                            for (int i = 0; i < n; ++i) {
                                V_n[i][jV_n] = 0.0;
                            }
                            V_nMap[jV_n] = 1;
                        } else {
                            for (int i = 0; i < n; ++i) {
                                aPV[i] = 0.0;
                                cPV[i] = 0.0;
                                ePV[i] = V_n[i][jV_n];
                                roV_n = 0.5 * (ro[i][j] + ro[i][jf]);
                                bPV[i] = roV_n * (1.0 / dt + effK_r[i][jV_n] * V_nz(i, j, jf, jV_n));
                                dPV[i] = 0.5 * (OLD_ro[i][j] + OLD_ro[i][jf]) / dt * oldV_n[i][jV_n] - 
                                        is(jf - j) * (p[i][jf] - p[i][j]) / dr;
                                V_up = 0.5 * fz_vf * (V_z[i + 1][j] + V_z[i + 1][jf]);
                                V_down = -0.5 * fz_vf * (V_z[i][j] + V_z[i][jf]);

                                if (V_up > 0.0) {
                                    bPV[i] += V_up * roV_n;
                                } else {
                                    if(i == n - 1) {
                                        dPV[i] -= V_up * ro_output * V_n[i][jV_n];
                                    } else {
                                        cPV[i] += V_up * (ro[i + 1][j] + ro[i + 1][jf]) * 0.5;
                                    }
                                }

                                if(V_down > 0) {
                                    bPV[i] += V_down * roV_n;
                                } else {
                                    if(i == 0) {
                                        dPV[i] -= V_down * ro_input[j] * V_n[i][jV_n];
                                    } else {
                                        aPV[i] += V_down * (ro[i - 1][j] + ro[i - 1][jf]) / 2;
                                    }
                                }
                            }
                            int j1 = bonds[0][j];
                            int j2 = bonds[1][j];
                            int j3 = bonds[2][j];
                            int j4 = bonds[3][j];
                            int j5 = bonds[4][j];
                            int j6 = bonds[5][j];

                            int jf1 = bonds[0][jf];
                            int jf2 = bonds[1][jf];
                            int jf3 = bonds[2][jf];
                            int jf4 = bonds[3][jf];
                            int jf5 = bonds[4][jf];
                            int jf6 = bonds[5][jf];
                            double V3;

                            switch (k) {
                                case 0:
                                    for(int i = 0; i < n; i++) {
                                        //a
                                        V3 = (is(j3-j) * V_n[i][onds[2][j]] - is(j6-j) * V_n[i][onds[5][j]] 
                                            + is(jf3-jf) * V_n[i][onds[2][jf]]) * fr_vf;
                                        if (V3 > 0.0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            if (j3 != -1) {
                                                Vx += is(jf - j) * is(j2 - j3) * V_n[i][onds[k][j3]] / 2.0;
                                            }
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }
                                        
                                        //b
                                        V3 = (is(j5 - j) * V_n[i][onds[4][j]] - is(j2 - j) * V_n[i][onds[1][j]] 
                                            + is(jf5 - jf) * V_n[i][onds[4][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] = bPV[i] + V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j5 != -1) {
                                                Vx += is(jf - j) * is(j6 - j5) * V_n[i][onds[k][j5]] / 2.0;
                                            }
                                            
                                            dPV[i] = dPV[i] - V3 * ro[i][j] * Vx;
                                        }

                                        //c
                                        V3 = (is(jf2 - jf) * V_n[i][onds[1][jf]] - is(jf5 - jf) * V_n[i][onds[4][jf]] 
                                            + is(j2 - j) * V_n[i][onds[1][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] = bPV[i] + V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf3 != -1) {
                                                Vx += is(jf - j) * is(jf2 - jf3) * V_n[i][onds[k][jf3]] / 2.0;
                                            }
                                            
                                            dPV[i] = dPV[i] - V3 * ro[i][jf] * Vx;
                                        }

                                        //d
                                        V3 = (is(jf6 - jf) * V_n[i][onds[5][jf]] - is(jf3 - jf) * V_n[i][onds[2][jf]] 
                                            + is(j6 - j) * V_n[i][onds[5][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] = bPV[i] + V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf5 != -1) {
                                                Vx += is(jf - j) * is(jf6 - jf5) * V_n[i][onds[k][jf5]] / 2.0;
                                            }
                                            
                                            dPV[i] = dPV[i] - V3 * ro[i][jf] * Vx;
                                        }
                                    }
                                    break;

                                case 1:
                                    for(int i = 0; i < n; ++i) {
                                        //a
                                        V3 = (is(j4 - j) * V_n[i][onds[3][j]] - is(j1 - j) * V_n[i][onds[0][j]] 
                                            + is(jf4 - jf) * V_n[i][onds[3][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j4 != -1) {
                                                Vx += is(jf - j) * is(j3 - j4) * V_n[i][onds[k][j4]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //b
                                        V3 = (is(j6 - j) * V_n[i][onds[5][j]] - is(j3 - j) * V_n[i][onds[2][j]] 
                                            + is(jf6 - jf) * V_n[i][onds[5][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j6 != -1) {
                                                Vx += is(jf - j) * is(j1 - j6) * V_n[i][onds[k][j6]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //c
                                        V3 = (is(jf3 - jf) * V_n[i][onds[2][jf]] - is(jf6 - jf) * V_n[i][onds[5][jf]] 
                                            + is(j3 - j) * V_n[i][onds[2][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf4 != -1) {
                                                Vx += is(jf - j) * is(jf3 - jf4) * V_n[i][onds[k][jf4]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }

                                        //d
                                        V3 = (is(jf1 - jf) * V_n[i][onds[0][jf]] - is(jf4 - jf) * V_n[i][onds[3][jf]] 
                                            + is(j1 - j) * V_n[i][onds[0][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf6 != -1) {
                                                Vx += is(jf - j) * is(jf1 - jf6) * V_n[i][onds[k][jf6]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }
                                    }
                                    break;

                                case 2:
                                    for(int i = 0; i < n; i++) {
                                        //a
                                        V3 = (is(j5 - j) * V_n[i][onds[4][j]] - is(j2 - j) * V_n[i][onds[1][j]] 
                                            + is(jf5 - jf) * V_n[i][onds[4][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j5 != -1) {
                                                Vx += is(jf - j) * is(j4 - j5) * V_n[i][onds[k][j5]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //b
                                        V3 = (is(j1 - j) * V_n[i][onds[0][j]] - is(j4 - j) * V_n[i][onds[3][j]] 
                                            + is(jf1 - jf) * V_n[i][onds[0][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j1 != -1) {
                                                Vx += is(jf - j) * is(j2 - j1) * V_n[i][onds[k][j1]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //c
                                        V3 = (is(jf4 - jf) * V_n[i][onds[3][jf]] - is(jf1 - jf) * V_n[i][onds[0][jf]] 
                                            + is(j4 - j) * V_n[i][onds[3][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf5 != -1) {
                                                Vx += is(jf - j) * is(jf4 - jf5) * V_n[i][onds[k][jf5]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }

                                        //d
                                        V3 = (is(jf2 - jf) * V_n[i][onds[1][jf]] - is(jf5 - jf) * V_n[i][onds[4][jf]] 
                                            + is(j2 - j) * V_n[i][onds[1][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf1 != -1) {
                                                Vx += is(jf - j) * is(jf2 - jf1) * V_n[i][onds[k][jf1]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }
                                    }
                                    break;

                                case 3:
                                    for(int i = 0; i < n; i++) {
                                        //a
                                        V3 = (is(j6 - j) * V_n[i][onds[5][j]] - is(j3 - j) * V_n[i][onds[2][j]] 
                                            + is(jf6 - jf) * V_n[i][onds[5][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j6 != -1) {
                                                Vx += is(jf - j) * is(j5 - j6) * V_n[i][onds[k][j6]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //b
                                        V3 = (is(j2 - j) * V_n[i][onds[1][j]] - is(j5 - j) * V_n[i][onds[4][j]] 
                                            + is(jf2 - jf) * V_n[i][onds[1][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j2 != -1) {
                                                Vx += is(jf - j) * is(j3 - j2) * V_n[i][onds[k][j2]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //c
                                        V3 = (is(jf5 - jf) * V_n[i][onds[4][jf]] - is(jf2 - jf) * V_n[i][onds[1][jf]] 
                                            + is(j5 - j) * V_n[i][onds[4][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf6 != -1) {
                                                Vx += is(jf - j) * is(jf5 - jf6) * V_n[i][onds[k][jf6]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }

                                        //d
                                        V3 = (is(jf3 - jf) * V_n[i][onds[2][jf]] - is(jf6 - jf) * V_n[i][onds[5][jf]] 
                                            + is(j3 - j) * V_n[i][onds[2][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf2 != -1) {
                                                Vx += is(jf - j) * is(jf3 - jf2) * V_n[i][onds[k][jf2]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }
                                    }
                                    break;

                                case 4:
                                    for(int i = 0; i < n; i++) {
                                        //a
                                        V3 = (is(j1 - j) * V_n[i][onds[0][j]] - is(j4 - j) * V_n[i][onds[3][j]] 
                                            + is(jf1 - jf) * V_n[i][onds[0][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j1 != -1) {
                                                Vx += is(jf - j) * is(j6 - j1) * V_n[i][onds[k][j1]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //b
                                        V3 = (is(j3 - j) * V_n[i][onds[2][j]] - is(j6 - j) * V_n[i][onds[5][j]] 
                                            + is(jf3 - jf) * V_n[i][onds[2][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j3 != -1) {
                                                Vx += is(jf - j) * is(j4 - j3) * V_n[i][onds[k][j3]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //c
                                        V3 = (is(jf6 - jf) * V_n[i][onds[5][jf]] - is(jf3 - jf) * V_n[i][onds[2][jf]] 
                                            + is(j6 - j) * V_n[i][onds[5][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf1 != -1) {
                                                Vx += is(jf - j) * is(jf6 - jf1) * V_n[i][onds[k][jf1]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }

                                        //d
                                        V3 = (is(jf4 - jf) * V_n[i][onds[3][jf]] - is(jf1 - jf) * V_n[i][onds[0][jf]] 
                                            + is(j4 - j) * V_n[i][onds[3][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf3 != -1) {
                                                Vx += is(jf - j) * is(jf4 - jf3) * V_n[i][onds[k][jf3]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }
                                    }
                                    break;
                                
                                case 5:
                                    for(int i = 0; i < n; i++) {
                                        //a
                                        V3 = (is(j2 - j) * V_n[i][onds[1][j]] - is(j5 - j) * V_n[i][onds[4][j]] 
                                            + is(jf2 - jf) * V_n[i][onds[1][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j2 != -1) {
                                                Vx += is(jf - j) * is(j1 - j2) * V_n[i][onds[k][j2]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //b
                                        V3 = (is(j4 - j) * V_n[i][onds[3][j]] - is(j1 - j) * V_n[i][onds[0][j]] 
                                            + is(jf4 - jf) * V_n[i][onds[3][jf]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][j];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (j4 != -1) {
                                                Vx += is(jf - j) * is(j5 - j4) * V_n[i][onds[k][j4]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][j] * Vx;
                                        }

                                        //c
                                        V3 = (is(jf1 - jf) * V_n[i][onds[0][jf]] - is(jf4 - jf) * V_n[i][onds[3][jf]] 
                                            + is(j1 - j) * V_n[i][onds[0][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf2 != -1) {
                                                Vx += is(jf - j) * is(jf1 - jf2) * V_n[i][onds[k][jf2]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }

                                        //d
                                        V3 = (is(jf5 - jf) * V_n[i][onds[4][jf]] - is(jf2 - jf) * V_n[i][onds[1][jf]] 
                                            + is(j5 - j) * V_n[i][onds[4][j]]) * fr_vf;

                                        if (V3 > 0) {
                                            bPV[i] += V3 * ro[i][jf];
                                        } else {
                                            double Vx = V_n[i][onds[k][j]] / 2.0;
                                            
                                            if (jf4 != -1) {
                                                Vx += is(jf - j) * is(jf5 - jf4) * V_n[i][onds[k][jf4]] / 2.0;
                                            }
                                            
                                            dPV[i] -= V3 * ro[i][jf] * Vx;
                                        }
                                    }
                                    break;
                            }

                            sy(aPV, bPV, cPV, dPV, 0, n - 1);
                            for(int i = 0; i < n; i++) {
                                V_n[i][jV_n] = dPV[i];
                                if (ePV[i] > 0.000002 || ePV[i] < -0.000002) {
                                    e = abs(1.0 - V_n[i][jV_n] / ePV[i]);
                                    if (e > error) error = e;
                                }
                            }
                        }
                        V_nMap[jV_n] = 1;
                    }
                }
            }
        }
    }

    if(iii == 1500) std::cout << "Piter. PVTerror > eps!" << std::endl;
}

void pes(double dt) {

    for(int jV_n = 0; jV_n < mV_n; jV_n++) {
        V_nMap[jV_n] = 0;
    }
    

    for(int j = 0; j < mf; j++) {
                                    /*********dp_dz**********/
        int i = 0; 
        switch (blockade[j]) {
            case 1:
                dp_dz[i][j] = ro[i][j] * V_z[i][j] / dt;
                break;
            case 0:
                dp_dz[i][j] = OLD_ro[i][j] * oldV_z[i][j] / dt -
                            ro[i][j] * V_z[i][j] * effK_z[i][j] * std::abs(V_z[i][j]);
                double V_up = fz_vf * (V_z[i][j] + V_z[i + 1][j]);
                double V_down = -2.0 * fz_vf * V_z[i][j];

                if (V_up > 0) {
                    dp_dz[i][j] -= V_up * ro[i][j] * V_z[i][j];
                } else {
                    dp_dz[i][j] -= V_up * ro[i][j] * V_z[i + 1][j];
                }

                if (V_down > 0) {
                    dp_dz[i][j] -= V_down * ro[i][j] * V_z[i][j];
                } else {
                    dp_dz[i][j] -= V_down * ro_input[j] * V_z[i][j];
                }

                for (int k = 0; k < nbf; ++k) {
                    int jf = bonds[k][j];
                    if (jf != -1) {
                        double Vn = fr_vf * is(jf - j) * V_n[i][onds[k][j]];
                        double cr_M =  fr_vf * (effM[i][j] + effM[i][jf]) / 2.0 / dr;
                        dp_dz[i][j] += cr_M * (V_z[i][jf] - V_z[i][j]);

                        if (Vn > 0) {
                            dp_dz[i][j] -= Vn * ro[i][j] * V_z[i][j];
                        } else {
                            dp_dz[i][j] -= Vn * ro[i][jf] * V_z[i][jf];
                        }
                    }
                }
                break;
        }
            //модуль полной скорости
        for (int i = 1; i < n; ++i) {
            dp_dz[i][j] = 0.5 * (OLD_ro[i][j] + OLD_ro[i - 1][j]) * oldV_z[i][j] / dt -
                        0.5 * (ro[i][j] + ro[i - 1][j]) * V_z[i][j] * effK_z[i][j] *
                        std::abs(V_z[i][j]);

            double V_up = fz_vf * (V_z[i][j] + V_z[i + 1][j]) / 2.0;
            double V_down = -fz_vf * (V_z[i][j] + V_z[i - 1][j]) / 2.0;

            if (V_up > 0) {
                dp_dz[i][j] -= V_up * ro[i][j] * V_z[i][j];
            } else {
                dp_dz[i][j] -= V_up * ro[i][j] * V_z[i + 1][j];
            }

            if (V_down > 0) {
                dp_dz[i][j] -= V_down * ro[i - 1][j] * V_z[i][j];
            } else {
                dp_dz[i][j] -= V_down * ro[i - 1][j] * V_z[i - 1][j];
            }
        }

        for (int k = 0; k < nbf; ++k) {
            int jf = bonds[k][j];
            if (jf != -1) {
                int jV_n = onds[k][j];
                for (int i = 1; i < n; ++i) {
                    double Vn = fr_vf * is(jf - j) * (V_n[i - 1][jV_n] + V_n[i][jV_n]) / 2.0;
                    double cr_M = fr_vf * ((effM[i][j] + effM[i][jf]) / 2.0 )/ dr;
                    dp_dz[i][j] += cr_M * (V_z[i][jf] - V_z[i][j]);

                    if (Vn > 0) {
                        dp_dz[i][j] -= Vn * (ro[i - 1][j] + ro[i][j]) / 2.0  * V_z[i][j];
                    } else {
                        dp_dz[i][j] -= Vn * (ro[i - 1][jf] + ro[i][jf]) / 2.0 * V_z[i][jf];
                    }
                }
            }
        }
            //модуль полной скорости
        i = n;
        dp_dz[i][j] = OLD_ro[i - 1][j] * oldV_z[i][j] / dt -
                    ro[i - 1][j] * V_z[i][j] * effK_z[i][j] * std::abs(V_z[i][j]);

        double V_up = 2.0 * fz_vf * V_z[i][j];
        double V_down = -fz_vf * (V_z[i][j] + V_z[i - 1][j]);

        if (V_up > 0) {
            dp_dz[i][j] -= V_up * ro[i - 1][j] * V_z[i][j];
        } else {
            dp_dz[i][j] -= V_up * ro_output * V_z[i][j];
        }

        if (V_down > 0) {
            dp_dz[i][j] -= V_down * ro[i - 1][j] * V_z[i][j];
        } else {
            dp_dz[i][j] -= V_down * ro[i - 1][j] * V_z[i - 1][j];
        }

        for (int k = 0; k < nbf; ++k) {
            int jf = bonds[k][j];
            if (jf != -1) {
                int jV_n = onds[k][j];
                double Vn = fr_vf * is(jf - j) * V_n[i - 1][jV_n];
                double cr_M = fr_vf * (effM[i][j] + effM[i][jf]) / 2.0 / dr;
                dp_dz[i][j] += cr_M * (V_z[i][jf] - V_z[i][j]);

                if (Vn > 0) {
                    dp_dz[i][j] -= Vn * ro[i - 1][j] * V_z[i][j];
                } else {
                    dp_dz[i][j] -= Vn * ro[i - 1][jf] * V_z[i][jf];
                }
            }
        }

                                /******************dp_dn*****************/
        for(int k = 0; k < nbf; ++k) {
            double jf = bonds[k][j];
            if(jf == -1) continue;
            double jV_n = onds[k][j];

            if(V_nMap[jV_n] == 0) {
                for (int i = 0; i < n; ++i) {
                    double roV_n = 0.5 * (ro[i][j] + ro[i][jf]);
                    dp_dn[i][jV_n] = 0.5 * (OLD_ro[i][j] + OLD_ro[i][jf]) * oldV_n[i][jV_n] / dt -
                                    roV_n * V_n[i][jV_n] * effK_r[i][jV_n] * V_nz(i, j, jf, jV_n);

                    double V_up = fz_vf * (V_z[i + 1][j] + V_z[i + 1][jf]) / 2.0;
                    double V_down = -fz_vf * (V_z[i][j] + V_z[i][jf]) / 2.0;

                    if (V_up > 0) {
                        dp_dn[i][jV_n] -= V_up * roV_n * V_n[i][jV_n];
                    } else {
                        if (i == n - 1) {
                            dp_dn[i][jV_n] -= V_up * ro_output * V_n[i][jV_n];
                        } else {
                            dp_dn[i][jV_n] -= V_up * (ro[i + 1][j] + ro[i + 1][jf]) / 2.0 * V_n[i + 1][jV_n];
                        }
                    }

                    if (V_down > 0) {
                        dp_dn[i][jV_n] -= V_down * roV_n * V_n[i][jV_n];
                    } else {
                        if (i == 0) {
                            dp_dn[i][jV_n] -= V_down * ro_input[jV_n] * V_n[i][jV_n];
                        } else {
                            dp_dn[i][jV_n] -= V_down * (ro[i - 1][j] + ro[i - 1][jf]) / 2.0 * V_n[i - 1][jV_n];
                        }
                    }
                }


                int j1 = bonds[0][j];
                int j2 = bonds[1][j];
                int j3 = bonds[2][j];
                int j4 = bonds[3][j];
                int j5 = bonds[4][j];
                int j6 = bonds[5][j];  // Ограничение: нельзя задавать V_n на вогнутой границе.

                int jf1 = bonds[0][jf];
                int jf2 = bonds[1][jf];
                int jf3 = bonds[2][jf];
                int jf4 = bonds[3][jf];
                int jf5 = bonds[4][jf];
                int jf6 = bonds[5][jf];
                double V3;

                switch (k) {
                    case 0:
                        for (int i = 0; i < n; ++i) {
                            // a:
                            V3 = (is(j3 - j) * V_n[i][onds[2][j]] - is(j6 - j) * V_n[i][onds[5][j]] 
                                    + is(jf3 - jf) * V_n[i][onds[2][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j3 != -1) Vx += is(jf - j) * is(j2 - j3) * V_n[i][onds[k][j3]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            // b:
                            V3 = (is(j5 - j) * V_n[i][onds[4][j]] - is(j2 - j) * V_n[i][onds[1][j]] 
                                + is(jf5 - jf) * V_n[i][onds[4][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j5 != -1) Vx += is(jf - j) * is(j6 - j5) * V_n[i][onds[k][j5]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            // c:
                            V3 = (is(jf2 - jf) * V_n[i][onds[1][jf]] - is(jf5 - jf) * V_n[i][onds[4][jf]] 
                                + is(j2 - j) * V_n[i][onds[1][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf3 != -1) Vx += is(jf - j) * is(jf2 - jf3) * V_n[i][onds[k][jf3]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }

                            // d:
                            V3 = (is(jf6 - jf) * V_n[i][onds[5][jf]] - is(jf3 - jf) * V_n[i][onds[2][jf]] 
                                + is(j6 - j) * V_n[i][onds[5][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf5 != -1) Vx += is(jf - j) * is(jf6 - jf5) * V_n[i][onds[k][jf5]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }
                        }
                        break;

                    case 1:
                        for (int i = 0; i < n; ++i) {
                            // a:
                            V3 = (is(j4 - j) * V_n[i][onds[3][j]] - is(j1 - j) * V_n[i][onds[0][j]] 
                                    + is(jf4 - jf) * V_n[i][onds[3][jf]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j4 != -1) Vx += is(jf - j) * is(j3 - j4) * V_n[i][onds[k][j4]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            // b:
                            V3 = (is(j6 - j) * V_n[i][onds[5][j]] - is(j3 - j) * V_n[i][onds[2][j]] 
                                + is(jf6 - jf) * V_n[i][onds[5][jf]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j6 != -1) Vx += is(jf - j) * is(j1 - j6) * V_n[i][onds[k][j6]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            // c:
                            V3 = (is(jf3 - jf) * V_n[i][onds[2][jf]] - is(jf6 - jf) * V_n[i][onds[5][jf]] 
                                + is(j3 - j) * V_n[i][onds[2][j]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf4 != -1) Vx += is(jf - j) * is(jf3 - jf4) * V_n[i][onds[k][jf4]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }

                            // d:
                            V3 = (is(jf1 - jf) * V_n[i][onds[0][jf]] - is(jf4 - jf) * V_n[i][onds[3][jf]] 
                                + is(j1 - j) * V_n[i][onds[0][j]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf6 != -1) Vx += is(jf - j) * is(jf1 - jf6) * V_n[i][onds[k][jf6]] / 2.0;
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }
                        }
                        break;

                    case 2:
                        for (int i = 0; i < n; ++i) {
                            // a:
                            V3 = (is(j5 - j) * V_n[i][onds[4][j]] - is(j2 - j) * V_n[i][onds[1][j]]
                                    + is(jf5 - jf) * V_n[i][onds[4][jf]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j5 != -1) Vx += is(jf - j) * is(j4 - j5) * V_n[i][onds[k][j5]] / 2.0;
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][j] * Vx;
                            }

                            // b:
                            V3 = (is(j1 - j) * V_n[i][onds[0][j]] - is(j4 - j) * V_n[i][onds[3][j]]
                                + is(jf1 - jf) * V_n[i][onds[0][jf]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j1 != -1) Vx += is(jf - j) * is(j2 - j1) * V_n[i][onds[k][j1]] / 2.0;
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][j] * Vx;
                            }

                            // c:
                            V3 = (is(jf4 - jf) * V_n[i][onds[3][jf]] - is(jf1 - jf) * V_n[i][onds[0][jf]]
                                + is(j4 - j) * V_n[i][onds[3][j]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf5 != -1) Vx += is(jf - j) * is(jf4 - jf5) * V_n[i][onds[k][jf5]] / 2.0;
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][jf] * Vx;
                            }

                            // d:
                            V3 = (is(jf2 - jf) * V_n[i][onds[1][jf]] - is(jf5 - jf) * V_n[i][onds[4][jf]]
                                + is(j2 - j) * V_n[i][onds[1][j]]) * fr_vf;
                            if (V3 > 0) {
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf1 != -1) Vx += is(jf - j) * is(jf2 - jf1) * V_n[i][onds[k][jf1]] / 2.0;
                                dp_dn[i][jV_n] = dp_dn[i][jV_n] - V3 * ro[i][jf] * Vx;
                            }
                        }
                        break;

                    case 3:
                        for (int i = 0; i < n; ++i) {
                            //a
                            V3 = (is(j6 - j) * V_n[i][onds[5][j]] - is(j3 - j) * V_n[i][onds[2][j]]
                                + is(jf6 - jf) * V_n[i][onds[5][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j6 != -1) {
                                    Vx += is(jf - j) * is(j5 - j6) * V_n[i][onds[k][j6]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            //b
                            V3 = (is(j2 - j) * V_n[i][onds[1][j]] - is(j5 - j) * V_n[i][onds[4][j]]
                                + is(jf2 - jf) * V_n[i][onds[1][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j2 != -1) {
                                    Vx += is(jf - j) * is(j3 - j2) * V_n[i][onds[k][j2]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            //c
                            V3 = (is(jf5 - jf) * V_n[i][onds[4][jf]] - is(jf2 - jf) * V_n[i][onds[1][jf]]
                                + is(j5 - j) * V_n[i][onds[4][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf6 != -1) {
                                    Vx += is(jf - j) * is(jf5 - jf6) * V_n[i][onds[k][jf6]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }

                            //d
                            V3 = (is(jf3 - jf) * V_n[i][onds[2][jf]] - is(jf6 - jf) * V_n[i][onds[5][jf]]
                                + is(j3 - j) * V_n[i][onds[2][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf2 != -1) {
                                    Vx += is(jf - j) * is(jf3 - jf2) * V_n[i][onds[k][jf2]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }
                        }
                        break;
                    
                    case 4:
                        for (int i = 0; i < n; ++i) {
                            V3;

                            //a
                            V3 = (is(j1 - j) * V_n[i][onds[0][j]] - is(j4 - j) * V_n[i][onds[3][j]]
                                + is(jf1 - jf) * V_n[i][onds[0][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j1 != -1) {
                                    Vx += is(jf - j) * is(j6 - j1) * V_n[i][onds[k][j1]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            //b
                            V3 = (is(j3 - j) * V_n[i][onds[2][j]] - is(j6 - j) * V_n[i][onds[5][j]]
                                + is(jf3 - jf) * V_n[i][onds[2][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j3 != -1) {
                                    Vx += is(jf - j) * is(j4 - j3) * V_n[i][onds[k][j3]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            //c
                            V3 = (is(jf6 - jf) * V_n[i][onds[5][jf]] - is(jf3 - jf) * V_n[i][onds[2][jf]]
                                + is(j6 - j) * V_n[i][onds[5][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf1 != -1) {
                                    Vx += is(jf - j) * is(jf6 - jf1) * V_n[i][onds[k][jf1]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }

                            //d
                            V3 = (is(jf4 - jf) * V_n[i][onds[3][jf]] - is(jf1 - jf) * V_n[i][onds[0][jf]]
                                + is(j4 - j) * V_n[i][onds[3][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf3 != -1) {
                                    Vx += is(jf - j) * is(jf4 - jf3) * V_n[i][onds[k][jf3]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }
                        }
                        break;

                    case 5:
                        for (int i = 0; i < n; ++i) {
                            V3;

                            //a
                            V3 = (is(j2 - j) * V_n[i][onds[1][j]] - is(j5 - j) * V_n[i][onds[4][j]]
                                + is(jf2 - jf) * V_n[i][onds[1][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j2 != -1) {
                                    Vx += is(jf - j) * is(j1 - j2) * V_n[i][onds[k][j2]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            //b
                            V3 = (is(j4 - j) * V_n[i][onds[3][j]] - is(j1 - j) * V_n[i][onds[0][j]]
                                + is(jf4 - jf) * V_n[i][onds[3][jf]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (j4 != -1) {
                                    Vx += is(jf - j) * is(j5 - j4) * V_n[i][onds[k][j4]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][j] * Vx;
                            }

                            //c
                            V3 = (is(jf1 - jf) * V_n[i][onds[0][jf]] - is(jf4 - jf) * V_n[i][onds[3][jf]]
                                + is(j1 - j) * V_n[i][onds[0][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf2 != -1) {
                                    Vx += is(jf - j) * is(jf1 - jf2) * V_n[i][onds[k][jf2]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }

                            //d
                            V3 = (is(jf5 - jf) * V_n[i][onds[4][jf]] - is(jf2 - jf) * V_n[i][onds[1][jf]]
                                + is(j5 - j) * V_n[i][onds[4][j]]) * fr_vf;

                            if (V3 > 0) {
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * V_n[i][jV_n];
                            } else {
                                double Vx = V_n[i][onds[k][j]] / 2.0;
                                if (jf4 != -1) {
                                    Vx += is(jf - j) * is(jf5 - jf4) * V_n[i][onds[k][jf4]] / 2.0;
                                }
                                dp_dn[i][jV_n] -= V3 * ro[i][jf] * Vx;
                            }
                        }
                        break;
                }
                V_nMap[jV_n] = 1;
            }
        }
    }

//   ф                  ф
//   і div(grad P) dV = і grad(P) dF = source :
//   х                  х

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            source[i][j] = (dp_dz[i + 1][j] - dp_dz[i][j]) * fz + vf * (ro[i][j] - OLD_ro[i][j]) / (dt * dt);
        }

        
        source[0][j] = source[0][j] - Gravity * (ro[0][j] - ro_input[j]) * fz / dz;
        
        for (int i = 1; i < n; ++i) {
            source[i][j] = source[i][j] - Gravity * (ro[i][j] - ro[i - 1][j]) * fz / dz;
        }
        
        for (int k = 0; k < nbf; ++k) {
            int jf = bonds[k][j];
            if (jf == -1) continue;
            
            int jV_n = onds[k][j];
            for (int i = 0; i < n; ++i) {
                source[i][j] += is(jf - j) * dp_dn[i][jV_n] * fr;
            }
        }
    }
}


void piter() {
    double fr_dr = fr / dr;
    double fz_dz = fz / dz;
    int iii = 0;
    double error = 1.0;
    double x;
    int sum;
    std::vector<int> targ(mf);

    while (error > PVTerror && iii < 1500) {
        error = 0.0;
        iii++;

        std::fill(targ.begin(), targ.end(), 0);

        int sum = 0;

        while (sum < mf) {
            x = random_number();
            int j = static_cast<int>(mf * x);

            if (targ[j] == 1) {
                continue;
            } else {
                targ[j] = 1;
                sum++;

                for (int i = 0; i < n; i++) {
                    ePV[i] = p[i][j];
                }

                int i = 0;
                switch (blockade[j]) {
                    case 0:
                        bPV[i] = -fz_dz * 3.0;
                        cPV[i] = fz_dz;
                        dPV[i] = -fz_dz * 2.0 * p_input + source[i][j];
                        break;
                    case 1:
                        bPV[i] = -fz_dz;
                        cPV[i] = fz_dz;
                        dPV[i] = source[i][j];
                        break;
                }

                for (int i = 1; i < n - 1; i++) {
                    aPV[i] = fz_dz;
                    bPV[i] = -fz_dz * 2.0;
                    cPV[i] = fz_dz;
                    dPV[i] = source[i][j];
                }

                i = n - 1;
                aPV[i] = fz_dz;
                bPV[i] = -fz_dz * 3.0;
                dPV[i] = -fz_dz * 2.0 * p_output + source[i][j];

                for(int k = 0; k < nbf; k++) {
                    int jf = bonds[k][j];
                    if (jf == -1) continue;

                    for (int i = 0; i < n; i++) {
                        bPV[i] -= fr_dr;
                        dPV[i] -= fr_dr * p[i][jf];
                    }
                }

                sy(aPV, bPV, cPV, dPV, 0, n - 1);

                for (int i = 0; i < n; i++) {
                    p[i][j] = dPV[i];
                    if (ePV[i] == 0.0) continue;

                    double e = std::abs(1.0 - p[i][j] / ePV[i]);
                    if (e > error) error = e;
                }
            }
        }
    }

    if (iii == 1500) {
        std::cout << "Piter. PVTerror > eps!" << std::endl;
    }
}

template void KinViscosity(const Coolant<double>& coolant);

template<typename T>
void KinViscosity(const Coolant<T>& coolant) {
    double const_term = 2.0 * d_mesh * std::pow(x_mesh, 2) / std::numbers::pi * (2.0 + 0.115 / (x_mesh - 1.0)) * (x_mesh - 1.0);
    
    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n + 1; ++i) {
            double vel = std::abs(V_z[i][j]); 
            double temperature, pression;

            if (i == 0) {
                temperature = h_f[i][j];
                pression = p[i][j];
            } else if (i == n) {
                temperature = h_f[i - 1][j];
                pression = p[i - 1][j];
            } else {
                temperature = (h_f[i - 1][j] + h_f[i][j]) / 2.0;
                pression = (p[i - 1][j] + p[i][j]) / 2.0;
            }

            double viscosity = coolant.KinVis(temperature);

            double re = vel * d_mesh / viscosity;

            if (vel != 0.0) {
                effM[i][j] = coolant.KinVis(temperature) * const_term * vel / std::pow(re, 0.1);
            } else {
                effM[i][j] = coolant.KinVis(temperature) * viscosity;
            }
        }
    }
}

template void FormFriction(const Coolant<double>& coolant);

template<typename T>
void FormFriction(const Coolant<T>& coolant) {
    double a_mesh = 0.58 + 9.2 * (x_mesh - 1.0);   
    double formula = 0.57 + 0.18 * (x_mesh - 1.0) + 0.53 * (1.0 - exp(-a_mesh));
    double rows = dr / (0.6830127 * d_mesh * x_mesh);
    double Constant = 3.586;
    double zKoeff1 = 3.0;
    double zKoeff2 = 1.5;
    double d_hydraulic = 0.01;

    for (int jV_n = 0; jV_n < mV_n; ++jV_n) {
        V_nMap[jV_n] = 0;
    }

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n + 1; ++i) {
            double vel = std::abs(V_z[i][j]);
            if (vel == 0.0) {
                effK_z[i][j] = 0.0;
            } else {
                double temperature, pression;
                if (i == 0) {
                    temperature = h_f[i][j];
                    pression = p[i][j];
                } else if (i == n) {
                    temperature = h_f[i - 1][j];
                    pression = p[i - 1][j];
                } else {
                    temperature = (h_f[i - 1][j] + h_f[i][j]) / 2.0;
                    pression = (p[i - 1][j] + p[i][j]) / 2.0;
                }

                double viscosity = coolant.KinVis(temperature);
                double re = vel * d_hydraulic / viscosity;

                effK_z[i][j] = zKoeff1 * formula * std::pow(100 * re, -0.25) / (2.0 * d_hydraulic);
            }
        }

        for (int k = 0; k < nbf; ++k) {
            int jf = bonds[k][j];
            if (jf == -1) continue;

            int jV_n = onds[k][j];
            if (V_nMap[jV_n] == 0) {
                for (int i = 0; i < n; ++i) {
                    double vel = V_nz(i, j, jf, jV_n);
                    if (vel == 0.0) {
                        effK_r[i][jV_n] = 0.0;
                    } else {
                        double re = vel * d_mesh / coolant.KinVis((h_f[i][j] + h_f[i][jf]) / 2.0);
                        effK_r[i][jV_n] = Constant * (rows + 1.0) * std::pow(re, -0.27) / (2.0 * dr);
                    }
                }
                V_nMap[jV_n] = 1;
            }
        }

    }
}
