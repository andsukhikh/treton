#include <iostream>
#include <cmath>
#include <numbers>
#include <iomanip>

#include "GlobalVar.hpp"
#include "NamelistReader.hpp"
#include "Thechyco.hpp"
#include "Matrix.hpp"
#include "Heat.hpp"
#include "CoolantMaterials.hpp"


void RodOnce() {
    static int manager = 1;

    if (manager) {

        NLReader::NamelistReader nlr(input_dir + "//THEHYCO.INI");

        nlr.use_namelist("RodList");

        s_mesh =            nlr.get<double>("s_mesh", 1);
        d_mesh =            nlr.get<double>("d_mesh", 1);
        R_contact =         nlr.get<double>("R_contact", 1);
        clad_rc =           nlr.get<double>("clad_rc", 1);
        clad_l =            nlr.get<double>("clad_l", 1);
        wireStep =          nlr.get<double>("wireStep", 1);


        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < type; ++j)
            {
                tmp[i][j] = nlr.get<double>("tmp", 1.0, i);
            }
        }

        x_mesh = s_mesh / d_mesh;
        d_mesh /= 1000.0;

        for (int k = 0; k < type; ++k) {
            for (int l = 0; l < 4; ++l) {
                tmp[l][k] /= 1000.0;

            }
        }

        RodGeometry();
        manager = 0;
    }
}


void RodGeometry() {
    for (int j = 0; j < type; ++j) {
        cvd[0][j] = tmp[0][j];
        cvd[n_rod][j] = tmp[1][j];
        cvd[n_rod + 1][j] = tmp[2][j];
        cvd[n_rod + 2][j] = tmp[3][j];

        double zero = std::pow(cvd[0][j], 2);
        double eq_area = (std::pow(cvd[n_rod][j], 2) - zero) / (2 * (n_rod - 1));

        for (int i = 1; i < n_rod; ++i) {
            cvd[i][j] = std::sqrt(zero + (2 * (i + 1) - 3) * eq_area);
        }
    }

    for (int j = 0; j < type; ++j) {
        int i = 0;
        bundle[i][j] = cvd[i][j];
        for (int i = 1; i < n_rod - 1; ++i) {
            bundle[i][j] = (cvd[i][j] + cvd[i + 1][j]) / 2;
        }

        bundle[n_rod - 1][j] = cvd[n_rod][j];
        bundle[n_rod][j] = cvd[n_rod][j];
        bundle[n_rod + 1][j] = cvd[n_rod + 2][j];

        for (int i = 0; i < n_rod + 2; ++i) {
            area[i][j] = std::pow(cvd[i + 1][j], 2) - std::pow(cvd[i][j], 2);
            geo_left[i][j] = 2 * cvd[i][j] / area[i][j];
            geo_right[i][j] = 2 * cvd[i + 1][j] / area[i][j];
            area[i][j] = std::numbers::pi * area[i][j];
        }

        a_fuel[j] = std::numbers::pi * (std::pow(cvd[n_rod][j], 2) - std::pow(cvd[0][j], 2));
        a_clad[j] = std::numbers::pi * (std::pow(cvd[n_rod + 2][j], 2) - std::pow(cvd[n_rod][j], 2));
        P_rod[j] = 2 * std::numbers::pi * cvd[n_rod + 2][j];
    }
}


void rod_property() {
    // CLADDING

    for (int k = 0; k < type; ++k) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                double t = t_fuel[i][j][k] + 273;

                if (k == 0) {
                    // UO2
                    if (t > (2850 + 273)) {
                        std::cout << "ATTENTION! At (" << i + 1 << "," << j + 1 << ") UO2 fuel is melting." << std::endl;
                    }
                    double tau = t / 1000.0;

                    double density = 10960 * (0.99672 + 1.179 * 1e-5 * t - 2.429 * 1e-9 * std::pow(t, 2) + 1.219 * 1e-12 * std::pow(t, 3));          // ~ 7200-10620 [kg/m**3]
                    double heatcap_in_mol = 52.17 + 87.95 * tau - 84.24 * std::pow(tau, 2) + 31.54 * std::pow(tau, 3)
                                        - 2.63 * std::pow(tau, 4) + 0.71 * std::pow(tau, -2);            // Heat capacity

                    double heatcap = heatcap_in_mol / 0.270;
                    double conductivity = (100.0 / (7.5408 + 17.692 * tau + 3.6142 * std::pow(tau, 2)) +
                        (6400.0 / std::pow(tau, 2.5)) * exp(-16.35 / tau));

                    //double conductivity = 2.5;
                    //double density = 10960;
                    //double heatcap = 410;

                    fuel_rc[i][j][k] = density * heatcap;
                    fuel_l[i][j][k] = conductivity;             // Conductivity

                } else if (k == 1) {
                    // UC
                    double density = 10200.0;          // ~[kg/m**3]
                    double heatcap = 80.076 + 0.4202 * (t - 273) - 0.000307 * (t - 273) * (t - 273);  // Heat capacity
                    fuel_rc[i][j][k] = density * heatcap;
                    fuel_l[i][j][k] = 22.0 * cosh(1.0 - (t - 273) / 750.0);           // Conductivity
                }
            }
        }
    }
}


// ����� ������������ ����� ���� ������ � ��������� ��-�� �������� ������������ ������� �� �������
// �������� �������� ���������������� ������� ��� ����������� ������������� ��������
 
double Rod(double dt, std::vector<double>& old, std::vector<double>& New, double alf, 
           double tf, double Q_v, std::vector<double>& g_left, std::vector<double>& g_right, 
           std::vector<double>& cor, double l, double rc) {
    
    double Rod = 0.0;

    // FUEL
    int i = 0;
    a[i] = 0.0;
    c[i] = -l * g_right[i] / (cor[i+1] - cor[i]);
    
    for (i = 1; i < n_rod - 1; ++i) {
        a[i] = -l * g_left[i] / (cor[i] - cor[i-1]);
        c[i] = -l * g_right[i] / (cor[i+1] - cor[i]);
    }
    
    i = n_rod - 1;
    a[i] = -l * g_left[i] / (cor[i] - cor[i-1]);
    c[i] = -g_right[i] / R_contact;
    
    for (i = 0; i < n_rod; ++i) {
        b[i] = rc / dt - a[i] - c[i];
        d[i] = rc / dt * old[i] + Q_v;
    }
    
    // CLADDING
    i = n_rod;
    a[i] = -g_left[i] / R_contact;
    c[i] = -clad_l * g_right[i] / (cor[i+1] - cor[i]);
    b[i] = clad_rc / dt - a[i] - c[i];
    d[i] = clad_rc / dt * old[i];
    
    i = n_rod + 1;
    a[i] = -clad_l * g_left[i] / (cor[i] - cor[i-1]);
    c[i] = 0.0;
    b[i] = clad_rc / dt - a[i] + g_right[i] * alf;
    d[i] = clad_rc / dt * old[i] + g_right[i] * alf * tf;
    
    sy(a, b, c, d, 0, n_rod + 1);
    
    for (int i = 0; i < n_rod + 2; ++i) {
        double e = std::abs(1.0 - New[i] / d[i]);
        New[i] = d[i];
        Rod = std::max(Rod, e);
    }
    
    return Rod;
}


void rod_average() {
    RodOnce();  

    for (int k = 0; k < type; ++k) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                double fuel = 0.0;
                double clad = 0.0;

                for (int l = 0; l < n_rod; ++l) {
                    fuel += t_rod[l][i][j][k] * area[l][k];
                }

                for (int l = n_rod; l < n_rod + 2; ++l) {
                    clad += t_rod[l][i][j][k] * area[l][k];
                }

                t_fuel[i][j][k] = fuel / a_fuel[k];
                t_clad[i][j][k] = clad / a_clad[k];
            }
        }


    }
}


double EneRoDisbalance() {
    double EneRoDisbalance = 0.0; 
    int i_er = 0, j_er = 0, k_er = 0; 

    for (int k = 0; k < type; ++k) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                double Qf = alfa[i][j] * dz * P_rod[k] * (t_rod[n_rod + 1][i][j][k] - t_f[i][j]) * n_RodsInTBC[k];
                double error = std::abs(1 - Qf / Q_neutron[i][j][k]);

                if (EneRoDisbalance < error) {
                    EneRoDisbalance = error;
                    i_er = i;
                    j_er = j;
                    k_er = k;
                }
            }
        }
    }

    std::cout << '\n' <<  std::endl; 
    std::cout << "EneRoDisbalance in (" << i_er + 1 << "," << j_er + 1 << "," << k_er + 1 << "): " << 
                    std::setprecision(6) << std::fixed << EneRoDisbalance << std::endl;

    return EneRoDisbalance;
}


void separationMatrix2D(std::vector<std::vector<double>>& before, std::vector<double>& after, int numCol)
{
    for (int i = 0; i < n_rod + 2; ++i)
    {
        after[i] = before[i][numCol];
    }
}


void separationMatrix4D(std::vector<std::vector<std::vector<std::vector<double>>>>& before, std::vector<double>& after, int numCol, int numLayer, int numLayer4d)
{
    for (int i = 0; i < n_rod + 2; ++i)
    {
        after[i] = before[i][numCol][numLayer][numLayer4d];
    }
}


void swapOnNewValue(std::vector<double>& arrayValue, std::vector<std::vector<std::vector<std::vector<double>>>>& swapableValue, int numCol, int numLayer, int numLayer4d)
{
    for (int i = 0; i < n_rod + 2; ++i)
    {
        swapableValue[i][numCol][numLayer][numLayer4d] = arrayValue[i];
    }
}


void heat(double dt, const Coolant& coolant) {
    double tmp1 = 0.0;
    double error = 0.0;
    int iii = 0;
    int NM = n * mf;
    double p_r, t_r, h_r, v_r;

    for (int k = 0; k < type; ++k) {
        tmp1 += P_rod[k] * n_RodsInTBC[k];
    }

    for (int k = 0; k < ((mf + mV_n) * n - mf) * 2 + n * mf; ++k) {
        Mat_A[k] = 0;
    }

    InOut_f(coolant);

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            int kk = j * n + i;
            double alfadz_vf = alfa[i][j] * dz / vf;
            double tmp2 = 0.0;

            t_f[i][j] = coolant.Temperature(h_f[i][j]);
            Vect_X[kk] = h_f[i][j];

            for (int k = 0; k < type; ++k) {
                double Qv = Q_neutron[i][j][k] / (a_fuel[k] * dz * n_RodsInTBC[k]);

                separationMatrix4D(OLDt_rod, transposed_OLDt_rod, i, j, k);
                separationMatrix4D(t_rod, transposed_t_rod, i, j, k);
                separationMatrix2D(geo_left, transposed_geo_left, k);
                separationMatrix2D(geo_right, transposed_geo_right, k);
                separationMatrix2D(bundle, transposed_bundle, k);

                double RodError = Rod(dt, transposed_OLDt_rod, transposed_t_rod, alfa[i][j], t_f[i][j],
                                        Qv, transposed_geo_left, transposed_geo_right, transposed_bundle, fuel_l[i][j][k], fuel_rc[i][j][k]);
                swapOnNewValue(transposed_t_rod, t_rod, i, j, k);

                if (RodError > error) error = RodError;
                tmp2 += P_rod[k] * (t_rod[n_rod + 1][i][j][k] - t_f[i][j]) * n_RodsInTBC[k];
            }

            Mat_A[kk] = ro[i][j] / dt;
            Vect_B[kk] = OLDh_f[i][j] * OLD_ro[i][j] / dt + tmp2 * alfadz_vf;
        }

        int i = 0;
        int kk = j * n + i;
        double V_up = fz_vf * V_z[i + 1][j];
        double V_down = -fz_vf * V_z[i][j];

        if (V_up > 0) {
            Mat_A[kk] += V_up * ro[i][j];
        } else {
            Mat_A[NM + NC[kk]] += V_up * ro[i + 1][j];
        }
        if (V_down > 0) {
            Mat_A[kk] += V_down * ro[i][j];
        } else {
            Vect_B[kk] -= V_down * ro_input[j] * h_HeatExchangerOutput_new[j];
        }

        for (i = 1; i < n - 1; ++i) {
            kk = j * n + i;
            V_up = fz_vf * V_z[i + 1][j];
            V_down = -fz_vf * V_z[i][j];
            if (V_up > 0) {
                Mat_A[kk] += V_up * ro[i][j];
            } else {
                Mat_A[NM + NC[kk] + 1] += V_up * ro[i + 1][j];
            }
            if (V_down > 0) {
                Mat_A[kk] += V_down * ro[i][j];
            } else {
                Mat_A[NM + NC[kk]] += V_down * ro[i - 1][j];
            }
        }

        i = n - 1;
        kk = j * n + i;
        V_up = fz_vf * V_z[i + 1][j];
        V_down = -fz_vf * V_z[i][j];
        if (V_up > 0) {
            Mat_A[kk] += V_up * ro[i][j];
        } else {
            Vect_B[kk] -= V_up * ro_output * h_CoreOutput;
        }
        if (V_down > 0) {
            Mat_A[kk] += V_down * ro[i][j];
        } else {
            Mat_A[NM + NC[kk]] += V_down * ro[i - 1][j];
        }

        for (i = 0; i < n; ++i) {
            kk = j * n + i;
            int nk = (i == 0 || i == n - 1) ? 1 : 2;
            for (int k = 0; k < nbf; ++k) {
                int jf = bonds[k][j];
                if (jf != -1) {
                    int jV_n = onds[k][j];
                    double const_V = fr_vf * is(jf - j);
                    double const_L = 0.5 * fr_vf / dr;
                    double cr_V = const_V * V_n[i][jV_n];
                    double cr_L = const_L * (effL[i][j] + effL[i][jf]);
                    Vect_B[kk] += cr_L * (t_f[i][jf] - t_f[i][j]);
                    if (cr_V > 0) {
                        Mat_A[kk] += cr_V * ro[i][j];
                    } else {
                        Mat_A[NM + NC[kk] + nk] += cr_V * ro[i][jf];
                        nk++;
                    }
                }
            }
        }
    }

    ResCalc(n * mf, Mat_A, NC, NE, Vect_B, Vect_X, Vect_X1, PVTerror);

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            int kk = j * n + i;
            h_f[i][j] = Vect_X1[kk];

            t_f[i][j] = coolant.Temperature(h_f[i][j]);
        }
    }

}


double EnerFluiDisbalance() {
    double EnerFluiDisbalance = 0.0;
    int i_er = 0, j_er = 0;

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            double tmp1 = 0.0;
            for (int k = 0; k < type; ++k) {
                tmp1 -= P_rod[k] * (t_rod[n_rod + 1][i][j][k] - t_f[i][j]) * n_RodsInTBC[k];
            }
            aa[i] = alfa[i][j] * dz * tmp1;
            bb[i] = 0.0;
        }

        int i = 0;
        double V_up = fz * V_z[i + 1][j];
        double V_down = -fz * V_z[i][j];

        if (V_up > 0) {
            bb[i] += V_up * ro[i][j] * h_f[i][j];
        } else {
            aa[i] += V_up * ro[i + 1][j] * h_f[i + 1][j];
        }
        
        if (V_down > 0) {
            bb[i] += V_down * ro[i][j] * h_f[i][j];
        } else {
            aa[i] += V_down * ro_input[j] * h_HeatExchangerOutput_new[j];
        }

        for (i = 1; i < n - 1; ++i) {
            V_up = fz * V_z[i + 1][j];
            V_down = -fz * V_z[i][j];
            if (V_up > 0) {
                bb[i] += V_up * ro[i][j] * h_f[i][j];
            } else {
                aa[i] += V_up * ro[i + 1][j] * h_f[i + 1][j];
            }
            if (V_down > 0) {
                bb[i] += V_down * ro[i][j] * h_f[i][j];
            } else {
                aa[i] += V_down * ro[i - 1][j] * h_f[i - 1][j];
            }
        }

        i = n - 1;
        V_up = fz * V_z[i + 1][j];
        V_down = -fz * V_z[i][j];
        if (V_up > 0) {
            bb[i] += V_up * ro[i][j] * h_f[i][j];
        } else {
            aa[i] += V_up * ro_output * h_CoreOutput;
        }
        if (V_down > 0) {
            bb[i] += V_down * ro[i][j] * h_f[i][j];
        } else {
            aa[i] += V_down * ro[i - 1][j] * h_f[i - 1][j];
        }

        for (int k = 0; k < nbf; ++k) {
            int jf = bonds[k][j];
            if (jf != -1) {
                int jV_n = onds[k][j];
                for (i = 0; i < n; ++i) {
                    double cr_V = fr * is(jf - j) * V_n[i][jV_n];
                    if (cr_V > 0) {
                        bb[i] += cr_V * ro[i][j] * h_f[i][j];
                    } else {
                        aa[i] += cr_V * ro[i][jf] * h_f[i][jf];
                    }

                    double cr_L = -0.5 * fr * (effL[i][j] + effL[i][jf]) * (t_f[i][jf] - t_f[i][j]) / dr;
                    if (cr_L > 0) {
                        bb[i] += cr_L;
                    } else {
                        aa[i] += cr_L;
                    }
                    
                }
            }
        }

        for (i = 0; i < n; ++i) {
            double error = std::abs(1.0 + aa[i] / bb[i]);
            if (EnerFluiDisbalance < error) {
                EnerFluiDisbalance = error;
                i_er = i;
                j_er = j;
            }
        }
    }

    std::cout << "EnerFluiDisbalance in (" << i_er << "," << j_er << "): " << 
        std::setprecision(7) << std::fixed << EnerFluiDisbalance << std::endl;

    return EnerFluiDisbalance;
}


void HeatConduction(const Coolant& coolant) {
    double const_val = 2.0 * d_mesh * std::pow(x_mesh, 2) / std::numbers::pi * (2.0 + 0.115 / (x_mesh - 1.0)) * (x_mesh - 1.0);

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            double vel = std::abs(V_z[i][j]);
            double viscosity = coolant.KinVis(h_f[i][j]);
            double re = vel * d_mesh / viscosity;

            if (vel != 0) {
                effL[i][j] = C_p[i][j] * ro[i][j] * const_val * vel / std::pow(re, 0.1);
            } else {
                effL[i][j] = coolant.HeatCond(coolant.Temperature(h_f[i][j]));
            }
        }
    }
}


void alf(const Coolant& coolant) {
    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            double Uz = (V_z[i][j] + V_z[i + 1][j]) / 2.0;
            double Ux = 0.0;
            double Uy = 0.0;
            alfa[i][j] = HeatTransfer(Ux, Uy, Uz, x_mesh, d_mesh, h_f[i][j], coolant);
        }
    }
}


double EnerCoreDisbalance() {
    double Energy = 0.0;
    double bottom = 0.0;
    double top = 0.0;

    for (int j = 0; j < mf; ++j) {
        if (V_z[0][j] > 0) {
            bottom += V_z[0][j] * ro_input[j] * h_HeatExchangerOutput_new[j];
        } else {
            bottom += V_z[0][j] * ro[0][j] * h_f[0][j];
        }

        if (V_z[n][j] > 0) {
            top += V_z[n][j] * ro[n - 1][j] * h_f[n - 1][j];
        } else {
            top += V_z[n][j] * ro_output * h_CoreOutput;
        }
    }

    for (int k = 0; k < type; ++k) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                Energy += Q_neutron[i][j][k];
            }
        }
    }

    double EnerCoreDisbalance = std::abs(1.0 - (top - bottom) * fz / Energy);

    std::cout << "EnerCoreDisbalance :  " << std::fixed << std::setprecision(6) << EnerCoreDisbalance << std::endl;

    return EnerCoreDisbalance;
}


double HeatTransfer(double Ux, double Uy, double Uz, double x, double d_hyd, double ent, const Coolant& coolant) {
    double t_r = coolant.Temperature(ent);
    double heat_cond = coolant.HeatCond(t_r);
    
    double d_hydro = d_hyd * (1.103 * std::pow(x_mesh, 2) - 1.0);
    double ro_Sodium = coolant.density(t_r);

    double absU = std::sqrt(Ux * Ux + Uy * Uy + Uz * Uz);

    double Re = absU * d_hydro * ro_Sodium / coolant.DynVisc(t_r);

    double Pr = coolant.Pr(t_r);
    
    // double Pe_along = absU * d_hydro / a_Pb;
    // double NU_along = NUl(x) + (0.041 / std::pow(x, 2)) * std::pow(Pe_along, 0.56 + 0.19 * x);
    // double along = NU_along * l_Pb / d_hydro;

    // // ���������� ���������
    // double sinfi = std::sqrt(Ux * Ux + Uy * Uy) / absU;
    // double Pe_cross = absU * d / a_Pb;
    // double cross = 0.0;

    // if (sinfi != 0) {
    //     double NUfi = 2 * std::sqrt(Pe_cross * (1 + sinfi) / (1 / sinfi + sinfi));
    //     cross = NUfi * l_Pb / d;
    // }

    // // ��������� ����������� � ����������� ���������
    // if (along > cross) {
    //     return along;
    // } else {
    //     return cross;
    // }

    return coolant.Nu(Re, Pr) * heat_cond / d_hydro;

}

