#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "headers/ThechycoGlobalVar.hpp"
#include "headers/SodiumProp.hpp"
#include "headers/Heat.hpp"
#include "headers/Thechyco.hpp"
#include "headers/Hydro.hpp"
#include "headers/NamelistReader.hpp"


int main() {

    double icall = 0.0;
    std::vector<std::vector<double>> bes(mf, std::vector<double>(n));
    std::vector<double> xx(mf);
    std::vector<double> yy(mf);
    std::vector<double> z(31);

    double p_r, t_r, h_r, v_r;

    //дописать код, использующий namelist 

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

    //std::cout << p_output << std::endl;


    crd = 
    {
        {
            -5, -3, -1, 1, 3, 5,
            -8, -6, -4, -2, 0, 2, 4, 6, 8,
            -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 
            -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10,
            -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 
            -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 
            -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 
            -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12,
            -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 
            -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12,
            -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 
            -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 
            -9, -7, -5, -3, -1, 1, 3, 5, 7, 9,  
            -8, -6, -4, -2, 0, 2, 4, 6, 8, 
            -5, -3, -1, 1, 3, 5
        },
        {
            7, 7, 7, 7, 7, 7, 
            6, 6, 6, 6, 6, 6, 6, 6, 6,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 
            -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 
            -6, -6, -6, -6, -6, -6, -6, -6, -6, 
            -7, -7, -7, -7, -7, -7
        }
    };


    p_r = p_input;
    //double h_HeatExchangerOutput2 = 1268500.0;

    std::ifstream file_in("T_in.txt");  
    if(file_in.is_open()) {
        for(int j = 0; j < mf; ++j) {
            for(int i = 0; i < 1; ++i) {
                file_in >> bes[j][0];
            }

            t_r = bes[j][0];
            SodiumEV(p_r, h_r, t_r, v_r);
            // WODAT(p_r, h_r, t_r, v_r);
        
            // h_HeatExchangerOutput_new[j] = h_r - h_HeatExchangerOutput2;
            h_HeatExchangerOutput_new[j] = h_r;

        }

        file_in.close();
    } 

    double VAU = v_r;

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            p[i][j] = p_input - (p_input - p_output) * (i + 0.5) / n;
            t_f[i][j] = bes[j][0]; 
            
            p_r =p[i][j];
            t_r = t_f[i][j];
      
            // WODAT(p_r, h_r, t_r, v_r);
            SodiumEV(p_r, h_r, t_r, v_r);
        
            // h_f[i][j] = h_r - h_HeatExchangerOutput2;
            h_f[i][j] = h_r;

            VAU = v_r;
        }

        
        for (int i = 0; i < n + 1; ++i) {
            V_z[i][j] = 5.0;
        }
    }

     for (int j = 0; j < mV_n; ++j) {
        for (int i = 0; i < n; ++i) {
            V_n[i][j] = 0.0;
        }
    }

    for (int k = 0; k < type; ++k) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                for (int l = 0; l < n_rod + 2; ++l) {
                    t_rod[l][i][j][k] = 1400.0 - 7.5 * std::pow(l, 2);
                }
                t_fuel[i][j][k] = 1100.0;
                t_clad[i][j][k] = 340.0;
            }
        }
    }

    z[0] = -dz * n / 2;

     for (int i = 0; i < n; ++i) {
        z[i + 1] = z[i] + dz;
    }

    //for (int i = 0; i < n - 1; ++i) {
    //    z[i + 1] = z[i] + dz;
    //}


    for (int j = 0; j < mf; ++j) {
        yy[j] = crd[1][j] * dr * 0.5 * sqrt(3.0);
        xx[j] = crd[0][j] * dr / 2;
    }

    std::ifstream Q6_file("Q6.txt");
    if (Q6_file.is_open()) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                Q6_file >> bes[j][i];  
            }
        }

        Q6_file.close();
    }

    double Q = 0.0;

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            Q_neutron[i][j][0] = bes[j][i];
            Q += Q_neutron[i][j][0];
        }
    }

    int icont_key;
    std::cout << "Enter 0 for new or 1 for continue: ";
     std::cin >> icont_key;
    //icont_key = 0;
    //icont_key = 0;

    if (icont_key == 1) {
        read_all();
        std::cout << "Main.cpp function read_all have completed successfully, maby)))" << std::endl;
    }

    RodOnce();
    //std::cout << "RodOnce in Main.cpp is completely" << std::endl;
    HeatHydroOnce();
    //std::cout << "HeatHydroOnce in Main.cpp is completely" << std::endl;
    V_zBlockade();
    //std::cout << "V_zBlockade in Main.cpp is completely" << std::endl;
    //std::cout << "PRESSURE = " << p[0][0] <<  std::endl;


    double time = 0.0;

    double dt;
    std::cout << "Enter dt = ";
    std::cin >> dt;
    //dt = 0.005;
    int kk = 100;

    for (int k = 1; k <= kk; ++k) {
        if (k % 20 == 0) {
            write_all();
        }
        std::cout << "                                    " << k << " of " << kk << std::endl;

        for (int i = 0; i < 100; ++i) {
            icall = thehyco(dt);
            // std::cout << "thehyco in Main.cpp is completely i = "<< i << std::endl;
            if (icall < 1) break;
            time += dt;
        }

        double av = 0.0;
        for (int j = 0; j < mf; ++j) {
            av += V_z[0][j];
        }
        av /= mf;

        // –асчет распределени€ компонент полного вектора скорости (Vx, Vy, Vz)
        V_full_calc();


        // ѕол€ скоростей в поперечном сечении, подготовленные дл€ векторного представлени€
        std::ofstream VxVy_z_file("VxVy_z.dat");
        for (int ii = 0; ii < n + 1; ++ii) {
            for (int j = 0; j < mf; ++j) {
                double Vr_nm;
                double Fi_nm;
                VrFi_nm(ii, j, V_full[0][ii][j], V_full[1][ii][j], Vr_nm, Fi_nm);

                VxVy_z_file.setf(std::ios::scientific | std::ios::right);
                VxVy_z_file.precision(6);

                VxVy_z_file << std::setw(14) << xx[j] 
                            << std::setw(14) << yy[j]
                            << std::setw(14) << Fi_nm
                            << std::setw(14) << Vr_nm << std::endl;
            }
        }
        VxVy_z_file.close();
      

        // ѕол€ скоростей в центральном продольном сечении, подготовленные дл€ векторного представлени€
        std::ofstream VxVz_file("VxVz.dat");
        for (int j = 76; j < 88; ++j) {
            for (int ii = 0; ii < n + 1; ++ii) {
                double Vr_nm;
                double Fi_nm;
                VrFi_nm(ii, j, V_full[0][ii][j], V_full[2][ii][j], Vr_nm, Fi_nm);

                VxVz_file.setf(std::ios::scientific | std::ios::right);
                VxVz_file.precision(6);

                VxVz_file << std::setw(14) << xx[j]
                        << std::setw(14) << z[ii]
                        << std::setw(14) << Fi_nm
                        << std::setw(14) << Vr_nm << std::endl;
            }
        }
        VxVz_file.close();

        // roVz_xz.dat
        std::ofstream roVz_xz_file("roVz_xz.dat");
        for (int j = 76; j < 88; ++j) {

            roVz_xz_file.setf(std::ios::scientific | std::ios::right);
            roVz_xz_file.precision(6);

            roVz_xz_file << std::setw(14) << xx[j] 
                         << std::setw(14) << z[0] 
                         << std::setw(14) << ro_input[j] * V_z[0][j] << std::endl;

            for (int ii = 1; ii < n; ++ii) {
                roVz_xz_file << std::setw(14) << xx[j]  
                             << std::setw(14) << z[ii]
                             << std::setw(14) << (ro[ii][j] + ro[ii - 1][j]) / 2 * V_z[ii][j] << std::endl;
            }
            roVz_xz_file << std::setw(14) << xx[j] 
                        << std::setw(14) 
                        << std::setw(14) << z[n] 
                        << std::setw(14) << ro_output * V_z[n][j] << std::endl;
        }
        roVz_xz_file.close();

        // T_xz.dat
        std::ofstream T_xz_file("T_xz.dat");
        for (int j = 76; j < 88; ++j) {
            for (int ii = 0; ii < n; ++ii) {

                T_xz_file.setf(std::ios::scientific | std::ios::right);
                T_xz_file.precision(6);

                T_xz_file << std::setw(14) << xx[j]
                            << std::setw(14) << z[ii] + 0.5 * dz  
                            << std::setw(14)  << t_f[ii][j] << std::endl;  // »спользуетс€ dz/2, но используйте реальный dz
            }
        }
        T_xz_file.close();

        // Vz77.dat
        std::ofstream Vz77_file("Vz77.dat");
        for (int ii = 0; ii < n + 1; ++ii) {
            Vz77_file << std::setw(5) << std::right << ii << " "
                        << std::setw(8) << std::setprecision(6) << V_z[ii][77] << std::endl;
        }
        Vz77_file.close();

        // Kord.dat
        std::ofstream Kord_file("Kord.dat");
        for (int j = 0; j < mf; ++j) {
            Kord_file << std::setw(5) << std::right << j << " "
                    << std::setw(8) << std::right << std::setprecision(5) << xx[j] 
                    << std::setw(8) << std::right << std::setprecision(5) << yy[j] << std::endl;
        }
        Kord_file.close();

        // V_n.dat
        std::ofstream V_n_file("V_n.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n + 1; ++ii) {

                V_n_file.setf(std::ios::right);

                V_n_file << std::setw(10) << std::setprecision(5) << " " << V_n[ii][j];
            }
            V_n_file << std::endl;
        }
        V_n_file.close();

        // Vz.dat
        std::ofstream Vz_file("Vz.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n + 1; ++ii) {

                Vz_file.setf(std::ios::right);

                Vz_file << std::setw(12) << std::setprecision(6) << xx[j] 
                        << std::setw(12) << std::setprecision(6) << yy[j] 
                        << std::setw(12) << std::setprecision(6) << V_z[ii][j];
            }
            Vz_file << std::endl;
        }
        Vz_file.close();

        // G.dat
        std::ofstream G_file("G.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n + 1; ++ii) {

                G_file.precision(6);
                G_file.setf(std::ios::right);

                G_file << std::setw(12) << xx[j] << std::setw(12) << yy[j]<< std::setw(12) << V_z[ii][j] * fz;
            }
            G_file << std::endl;
        }
        G_file.close();

        //G_inp.dat
        std::ofstream file("G_inp.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(6) << xx[j] 
                << std::setw(10) << std::setprecision(6) << yy[j]
                << std::setw(10) << std::setprecision(6) << V_z[0][j] * fz << std::endl;
        }
        file.close();

        // G_out.dat
        file.open("G_out.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(6) << xx[j]
                << std::setw(10) << std::setprecision(6) << yy[j]
                << std::setw(10) << std::setprecision(6) << V_z[n][j] * fz << std::endl;
        }
        file.close();

        // p_tepl.dat
        file.open("p_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(15) << std::setprecision(4) << xx[j]
                << std::setw(15) << std::setprecision(4) << yy[j];
            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(15) << std::setprecision(4) << p[ii][j];
            }
            file << std::endl;
        }
        file.close();

        // ro_tepl.dat
        file.open("ro_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(4) << xx[j]
                << std::setw(10) << std::setprecision(4) << yy[j];
            for (int ii= 0; ii < n; ++ii) {
                file << std::setw(10) << std::setprecision(4) << ro[ii][j];
            }
            file << std::endl;
        }
        file.close();

        // p_tepl_poln.dat
        file.open("p_tepl_poln.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(15) << std::setprecision(4) << xx[j] 
                << std::setw(15) << std::setprecision(4) << yy[j];
            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(15) << std::setprecision(4) << 3.47 * (n - ii) * ro[ii][j];
            }
            file << std::endl;
        }
        file.close();


        // t_tepl.dat
        file.open("t_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(4) << xx[j]
                << std::setw(10) << std::setprecision(4) << yy[j];
            for (int i = 0; i < n; ++i) {
                file << std::left << std::setw(10) << std::setprecision(4) << t_f[i][j];
            }
            file << std::endl;
        }
        file.close();

        // h_tepl.dat
        file.open("h_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(4) << xx[j]
                << std::setw(10) << std::setprecision(4) << yy[j];
            for (int i = 0; i < n; ++i) {
                file << std::setw(10) << std::setprecision(4) << h_f[i][j] / 1000.0;
            }
            file << std::endl;
        }
        file.close();

        // alfa.dat
        file.open("alfa.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(2) << 0.0
                << std::setw(10) << std::setprecision(2) << 0.0;
            for (int i = 0; i < n; ++i) {
                file << std::setw(10) << std::setprecision(2) << alfa[i][j];
            }
            file << std::endl;
        }
        file.close();

        // t_clad.dat
        file.open("t_clad.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(2) << 0
                << std::setw(10) << std::setprecision(2) << 0;
            for (int i = 0; i < n; ++i) {
                file << std::setw(10) << std::setprecision(2) << t_clad[i][j][0];
            }
            file << std::endl;
        }
        file.close();

        // t_fuel.dat
        file.open("t_fuel.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(2) << 0
                << std::setw(10) << std::setprecision(2) << 0;
            for (int i = 0; i < n; ++i) {
                file << std::setw(10) << std::setprecision(4) << t_fuel[i][j][0];
            }
            file << std::endl;
        }
        file.close();

        // n_rod.dat
        file.open("n_rod.dat");
        for (int ii = 0; ii < n_rod + 2; ++ii) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(4) << 1000 * bundle[ii][0] << std::endl;
        }
        file.close();

        // Q.dat
        file.open("Q.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(15) << std::setprecision(4) << xx[j]
                << std::setw(15) << std::setprecision(4) << yy[j];
            for (int i = 0; i < n; ++i) {
                file << std::setw(15) << std::setprecision(4) << Q_neutron[i][j][0];
            }
            file << std::endl;
        }
        file.close();

        // t_r77.dat
        file.open("t_r77.dat");
        for (int ii = 0; ii < n; ++ii) {
            file.setf(std::ios::right);
            for (int l = 0; l < n_rod + 2; ++l) {
                file << std::setw(10) << std::setprecision(4) << t_rod[l][ii][77][0];
            }
            file << std::endl;
        }
        file.close();

        // t_r_centr.dat
        file.open("t_r_centr.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::right);
            file << std::setw(10) << std::setprecision(4) << 0.0
                << std::setw(10) << std::setprecision(4) << 0.0;
            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(10) << std::setprecision(4) << t_rod[0][ii][j][0];
            }
            file << std::endl;
        }
        file.close();


        double tmax = 300.0;
        double t_tepl_max = 287.0;
        double t_obol_max = 287.0;
        int iimax = -1, jjmax = -1;
        int iiimax = -1, jjjmax = -1;
        int iiiimax = -1, jjjjmax = -1;

        file.open("t_r_centr.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n; ++ii) {
                if (t_rod[0][ii][j][0] > tmax) {
                    tmax = t_rod[0][ii][j][0];
                    iimax = ii;
                    jjmax = j;
                }

                if (t_f[ii][j] > t_tepl_max) {
                    t_tepl_max = t_f[ii][j];
                    iiimax = ii;
                    jjjmax = j;
                }

                if (t_rod[n_rod + 1][ii][j][0] > t_obol_max) {
                    t_obol_max = t_rod[n_rod + 1][ii][j][0];
                    iiiimax = ii;
                    jjjjmax = j;
                }
            }
        }

        file.close();

        std::ofstream XXX_file("XXX.dat");
        XXX_file << "Max T topl: " << tmax << "\n"
                 << "V kassete # " << jjmax << "\n"
                 << "V el-te # " << iimax << "\n"
                 << "Max T tepl: " << t_tepl_max << "\n"
                 << "V kassete # " << jjjmax << "\n"
                 << "V el-te # " << iiimax << "\n"
                 << "Max T obol: " << t_obol_max << "\n"
                 << "V kassete # " << jjjjmax << "\n"
                 << "V el-te # " << iiiimax << std::endl;
        
        XXX_file.close();

        // не пон€тно дл€ чего
        // for (int j = 0; j < mf; ++j) {
        //     for (int i = 0; i < n; ++i) {
        //         p_r = p[i][j];
        //         h_r = h_f[i][j];
        //         // LNAS(p_r, h_r);
        //     }
        // }


        // 1. t_rXXX.dat
        file.open("t_rXXX.dat");
        for (int ii = 0; ii < n; ++ii) {
            for (int l = 0; l < n_rod + 2; ++l) {
                file.setf(std::ios::right);
                file << std::setw(10) << std::setprecision(4) << t_rod[l][ii][jjmax][0];
            }
            file << std::endl;
        }
        file.close();

        // 2. t_r_obl_vnutr.dat
        file.open("t_r_obl_vnutr.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n; ++ii) {
                file.setf(std::ios::right);
                file << std::setw(10) << std::setprecision(4) << 0.0 
                    << std::setw(10) << std::setprecision(4) << 0.0 
                    << std::setw(10) << std::setprecision(4) << t_rod[n_rod][ii][j][0];
            }
            file << "\n";
        }
        file.close();

        // 3. t_r_obl_naruj.dat
        file.open("t_r_obl_naruj.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n; ++ii) {
                file.setf(std::ios::right);
                file << std::setw(10) << std::setprecision(4) << 0.0
                    << std::setw(10) << std::setprecision(4) << 0.0 
                    << std::setw(10) << std::setprecision(4) << t_rod[n_rod + 1][ii][j][0];
            }
            file << std::endl;
        }
        file.close();

        // 4. koeff_z.dat
        file.open("koeff_z.dat");
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                file.setf(std::ios::right);
                file << std::setw(10) << std::setprecision(4) << 0.0
                    << std::setw(10) << std::setprecision(4) << 0.0 
                    << std::setw(10) << std::setprecision(4) << effK_z[i][j];
            }
            file << std::endl;
        }
        file.close();

        // 5. Cp.dat
        file.open("Cp.dat");
        file << std::setw(10) << std::setprecision(4) << Cp_input 
            << std::setw(10) << std::setprecision(4) << Cp_output << std::endl;
        for (int j = 0; j < mf; ++j) {
            file << std::setw(10) << std::setprecision(4) << xx[j]
                << std::setw(10) << std::setprecision(4) << yy[j];
            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(10) << std::setprecision(4) << C_p[ii][j];
            }
            file << std::endl;
        }
        file.close();

        double G = 0;
        double G1 = 0;

        file.open("GQT.dat");
        for (int j = 0; j < mf; ++j) {
            G += V_z[0][j] * fz; 
            G1 += V_z[n][j] * fz;
        }

        file << std::setprecision(6) 
                << std::setw(12) << G << " "
                << std::setw(12) << G1 << " "
                << std::setw(12) << Q << " "
                << std::setw(12) << t_CoreOutput << " "
                << std::setw(12) << t_CoreInput << " "
                << std::setw(12) << (t_CoreOutput - t_CoreInput) << " "
                << std::setw(12) << ro_input[76] << " " 
                << std::setw(12) << ro_output << std::endl;
        file.close();
    
        if (icall < 1) {
            std::cout << "icall < 1 ==> exit(1)" << std::endl;
            exit(1);
        }
    }

    write_all();

    std::cout << "End of program" << std::endl;
    std::cin.get();

    return 0;
}