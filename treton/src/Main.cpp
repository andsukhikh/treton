#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include <exception>
#include <memory>
#include <unordered_map>
#include <filesystem>
#include <algorithm>
#include <chrono>

#include "GlobalVar.hpp"
#include "Heat.hpp"
#include "Heat.hpp"
#include "Thechyco.hpp"
#include "Hydro.hpp"
#include "NamelistReader.hpp"
#include "CoolantMaterials.hpp"
#include "BinFileTools.hpp"
#include "SolverEffCoef.hpp"

#include "CRD_test.cpp"


int main() {

    input_dir = "..//input";
    if (!std::filesystem::exists(input_dir)) {
        std::filesystem::create_directories(input_dir);
    }


    NLReader::NamelistReader nlr(input_dir + "//THEHYCO.INI");

    nlr.use_namelist("PartitionList");

    mf =                            nlr.get<int>("mf", 1);
    n =                             nlr.get<int>("n", 1);
    n_rod =                         nlr.get<int>("n_rod", 1);
    type =                          nlr.get<int>("type", 1);
    nbf =                           nlr.get<int>("nbf", 1);

    GlobalVariables::Resizing();

    nlr.use_namelist("HEATandHYDROlist");

    dr =                            nlr.get<double>("dr", 1);
    Height =                        nlr.get<double>("Height", 1);
    D_tube =                        nlr.get<double>("D_tube", 1);
    Disbalance =                    nlr.get<double>("Disbalance", 1);
    PVTerror =                      nlr.get<double>("PVTerror", 1);
    BlockadePorousity =             nlr.get<int>("BlockadePorousity", 1);
    p_input =                       nlr.get<double>("p_input", 1);
    p_output =                      nlr.get<double>("p_output", 1);
    iterations =                    nlr.get<int>("iterations", 1);
    coolantName =                   nlr.get<std::string>("coolant", "non-existent");
    resist_multiplier =             nlr.get<double>("resist_multiplier", 1);

    dz = static_cast<double>(Height / n);
    BlockadeThickness = nlr.get<double>("BlockadeThickness", dz);

    double order = std::pow(10, std::floor(std::log10(std::fabs(dz))));
    eps = std::round(dz / order) * order;

    if (BlockadeThickness <= dz) BlockadeThickness = (order + eps);


    for (int i = 0; i < type; ++i) {
        n_RodsInTBC[i] = nlr.get<int>("n_RodsInTBC", 1.0, i);
    }

    for (int i = 0; i < mf; ++i) {
        blockade[i] = nlr.get<int>("blockade", 1.0, i);
    }

    for (int i = 0; i < mf; ++i) {
        blockade_coord[i] = nlr.get<int>("blockade_coord", 1.0, i);
    }

    for (size_t i = 0; i < 2 * mf; ++i) {
        crd[i % 2][i / 2] = nlr.get<int>("crd", 0.0, i);
    }
	
    #ifdef TEST_CRD
        testCRD();
    #endif // TEST_CRD

    double icall = 0.0;
    
    try {
        CoolantDecriptor decript(coolantName);
        Coolant& coolant = decript.getCoolant();
    } catch (std::exception& exception) {
        std::cout << "Program was terminated due to: " << exception.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    CoolantDecriptor decript(coolantName);
    Coolant& coolant = decript.getCoolant();


    std::ifstream file_in( input_dir + "//T_in.txt");
    if(file_in.is_open()) {
        for(int j = 0; j < mf; ++j) {
            for(int i = 0; i < 1; ++i) {
                file_in >> bes[j][0];
            }
            h_HeatExchangerOutput_new[j] = (coolant).Entalpy(bes[j][0]);
        }
        file_in.close();
    } 

    for (int j = 0; j < mf; ++j) {
        for (int i = 0; i < n; ++i) {
            p[i][j] = p_input - (p_input - p_output) * (i + 0.5) / n;
            t_f[i][j] = bes[j][0];
            h_f[i][j] = (coolant).Entalpy(t_f[i][j]);
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
                    t_rod[l][i][j][k] = bes[j][0];
                }
                t_fuel[i][j][k] = bes[j][0];
                t_clad[i][j][k] = bes[j][0];
            }
        }
    }

    z[0] = -dz * static_cast<double>(n / 2);

     for (int i = 0; i < n; ++i) {
        z[i + 1] = z[i] + dz;
    }

    for (int j = 0; j < mf; ++j) {
        yy[j] = crd[1][j] * dr * 0.5 * std::sqrt(3.0);
        xx[j] = crd[0][j] * dr / 2;
    }


#ifdef ANALYTICS
    std::cout << std::endl;
    std::cout << "***********analitic energy distribution function mode is turned on***********" << std::endl;

    double sum_Q = 0.0;

    double K_z = nlr.get<double>("K_z", 1);
    double K_r = nlr.get<double>("K_r", 1);
    double Power = nlr.get<double>("Power", 1);
    auto PowerInCentre = Power * K_z * K_r / (1e-6 * n * mf);

    auto Radius = *(std::max_element(xx.begin(), xx.end()));

    double H_eff = EquationSolver(K_z, Height, EquationSolver::EquationType::Height).solve();
    double R_eff = EquationSolver(K_r, Radius, EquationSolver::EquationType::Radius).solve();

    std::cout << "                        ***********Power = " << Power << " ***********              " << "\n" << std::endl;

    auto analitic_distr = [&](double norm_coeff = 1.0)
        {
            sum_Q = 0.0;

            for (int j = 0; j < mf; ++j) {
                for (int i = 0; i < n; ++i) {
                    auto r = std::sqrt(std::pow(xx[j], 2) + std::pow(yy[j], 2));
                    bes[j][i] = norm_coeff * PowerInCentre * std::cyl_bessel_j(0, 2.41 * r / R_eff) * std::cos(std::numbers::pi * z[i] / H_eff);
                    sum_Q += bes[j][i];
                }
            }
        };

    analitic_distr();
    double eps = sum_Q * 1e-6 - Power;
    while (std::abs(eps) > 0.01)
    {
        analitic_distr(Power/(sum_Q * 1e-6));
        eps = sum_Q * 1e-6 - Power;
    }

    std::ofstream Q6(input_dir + "//Q6.txt");
    for (auto&& val1 : bes) {
        for (auto&& val2 : val1) {
            Q6 << std::setw(10) << std::right << std::fixed << std::setprecision(0) << val2;
        }
        Q6 << "\n";
    }
    Q6.close();
#else
    std::ifstream Q6_file(input_dir + "//Q6.txt");
    if (Q6_file.is_open()) {
        for (int j = 0; j < mf; ++j) {
            for (int i = 0; i < n; ++i) {
                Q6_file >> bes[j][i];
            }
        }
        Q6_file.close();
    }
#endif

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

    if (icont_key == 1) {
        Reader(std::ifstream(input_dir + "//data.dat", std::ios::binary)).read(p, V_z, V_n, h_f, t_f, t_rod, t_fuel, t_clad, time_);
    }

    RodOnce();
    HeatHydroOnce(coolant);
    V_zBlockade();

    double dt;
    std::cout << "Enter dt = ";
    std::cin >> dt;
    int k = 0;

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    while (true) {
        ++k;
        if (k % 20 == 0) {
            Writer(std::ofstream(input_dir + "//data.dat", std::ios::binary)).write(p, V_z, V_n, h_f, t_f, t_rod, t_fuel, t_clad, time_);
        }
        std::cout << "\n";
        std::cout << "                                    loop - " << k << std::endl;

        for (int i = 0; i < 100; ++i) {
            icall = thehyco(dt, coolant);
            if (icall < 1) break;
            time_ += dt;
        }

        end = std::chrono::high_resolution_clock::now();

        std::cout << "Calculation time = " << std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << " min" << std::endl;
        std::cout << "Reactor time = " << time_ << " sec" << std::endl;

        double av = 0.0;
        for (int j = 0; j < mf; ++j) {
            av += V_z[0][j];
        }
        av /= mf;

        // Расчет распределения компонент полного вектора скорости (Vx, Vy, Vz)
        V_full_calc();


        output_dir = "..//output";
        if (!std::filesystem::exists(output_dir)) {
            std::filesystem::create_directories(output_dir);
        }



        std::ofstream VxVy_z_file(output_dir + "//VxVy_z.dat");
        for (int ii = 0; ii < n + 1; ++ii) {
            for (int j = 0; j < mf; ++j) {
                double Vr_nm;
                double Fi_nm;
                VrFi_nm(ii, j, V_full[0][ii][j], V_full[1][ii][j], Vr_nm, Fi_nm);

                VxVy_z_file.setf(std::ios::scientific);
                VxVy_z_file.precision(6);

                VxVy_z_file << std::setw(14) << xx[j] 
                            << std::setw(14) << yy[j]
                            << std::setw(14) << Fi_nm
                            << std::setw(14) << Vr_nm << std::endl;
            }
        }
        VxVy_z_file.close();

        std::ofstream VxVz_file(output_dir + "//VxVz.dat");
        for (int j = 75; j < 88; ++j) {
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
        std::ofstream roVz_xz_file(output_dir + "//roVz_xz.dat");
        for (int j = 75; j < 88; ++j) {

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
        std::ofstream T_xz_file(output_dir + "//T_xz.dat");
        for (int j = 75; j < 88; ++j) {
            for (int ii = 0; ii < n; ++ii) {

                T_xz_file.setf(std::ios::scientific | std::ios::right);
                T_xz_file.precision(6);

                T_xz_file << std::setw(14) << xx[j]
                          << std::setw(14) << z[ii] + 0.5 * dz  
                          << std::setw(14)  << t_f[ii][j] << std::endl;  // Используется dz/2, но используйте реальный dz
            }
        }
        T_xz_file.close();

        // Vz77.dat
        std::ofstream Vz77_file(output_dir + "//Vz77.dat");
        for (int ii = 0; ii < n + 1; ++ii) {
            Vz77_file << std::setw(5) << std::right << ii + 1 << " "
                      << std::setw(8) << std::setprecision(6) << V_z[ii][77] << std::endl;
        }
        Vz77_file.close();

        // Kord.dat
        std::ofstream Kord_file(output_dir + "//Kord.dat");
        for (int j = 0; j < mf; ++j) {
            Kord_file.setf(std::ios::fixed);

            Kord_file << std::setw(5) << std::right << j + 1 << " "
                      << std::setw(9) << std::right << std::setprecision(5) << xx[j] 
                      << std::setw(9) << std::right << std::setprecision(5) << yy[j] << std::endl;
        }
        Kord_file.close();

        // V_n.dat
        std::ofstream V_n_file(output_dir + "//V_n.dat");
        for (int j = 0; j < mf; ++j) {
            for (int ii = 0; ii < n + 1; ++ii) {

                V_n_file.setf(std::ios::fixed);
                V_n_file.precision(5);
                V_n_file.width(10);

                V_n_file << V_n[ii][j];
            }
            V_n_file << std::endl;
        }
        V_n_file.close();

        // Vz.dat
        std::ofstream Vz_file(output_dir + "//Vz.dat");
        for (int j = 0; j < mf; ++j) {
            Vz_file.precision(6);
            Vz_file.setf(std::ios::fixed);
            Vz_file << std::setw(12) << xx[j]
                    << std::setw(12) << yy[j];

            for (int ii = 0; ii < n + 1; ++ii) {

                Vz_file << std::setw(12) << V_z[ii][j];
            }
            Vz_file << std::endl;
        }
        Vz_file.close();

        // G.dat
        std::ofstream G_file(output_dir + "//G.dat");
        for (int j = 0; j < mf; ++j) {
            G_file.precision(6);
            G_file.setf(std::ios::fixed);
            G_file << std::setw(12) << xx[j]
                   << std::setw(12) << yy[j];
            for (int ii = 0; ii < n + 1; ++ii) {
                G_file << std::setw(12) << V_z[ii][j] * fz;
            }
            G_file << std::endl;
        }
        G_file.close();

        //G_inp.dat
        std::ofstream file(output_dir + "//G_inp.dat");
        for (int j = 0; j < mf; ++j) {
            file.precision(6);
            file.setf(std::ios::fixed);

            file << std::setw(10) << xx[j]
                 << std::setw(10) << yy[j]
                 << std::setw(10) << V_z[0][j] * fz << std::endl;
        }
        file.close();

        // G_out.dat
        file.open(output_dir + "//G_out.dat");
        for (int j = 0; j < mf; ++j) {
            file.precision(6);
            file.setf(std::ios::fixed);

            file << std::setw(10) << xx[j]
                << std::setw(10) << yy[j]
                << std::setw(10) << V_z[n][j] * fz << std::endl;
        }
        file.close();

        // p_tepl.dat
        file.open(output_dir + "//p_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.precision(4);
            file.setf(std::ios::fixed);

            file << std::setw(15) << xx[j]
                 << std::setw(15) << yy[j];

            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(15) << p[ii][j];
            }
            file << std::endl;
        }
        file.close();

        // ro_tepl.dat
        file.open(output_dir + "//ro_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(12) << xx[j]
                 << std::setw(12) << yy[j];

            for (int ii= 0; ii < n; ++ii) {
                file << std::setw(12) << ro[ii][j];
            }
            file << std::endl;
        }
        file.close();

        // p_tepl_poln.dat
        file.open(output_dir + "//p_tepl_poln.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(15) << xx[j] 
                << std::setw(15) << yy[j];

            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(15) << 3.47 * (n - ii) * ro[ii][j];
            }
            file << std::endl;
        }
        file.close();


        // t_tepl.dat
        file.open(output_dir + "//t_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(12) << xx[j]
                << std::setw(12) << yy[j];

            for (int i = 0; i < n; ++i) {
                file << std::setw(12) << t_f[i][j];
            }
            file << std::endl;
        }
        file.close();

        // h_tepl.dat
        file.open(output_dir + "//h_tepl.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(12) << xx[j]
                 << std::setw(12) << yy[j];

            for (int i = 0; i < n; ++i) {
                file << std::setw(12) << h_f[i][j] / 1000.0;
            }
            file << std::endl;
        }
        file.close();

        // alfa.dat
        file.open(output_dir + "//alfa.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(2);

            file << std::setw(10) << 0.0
                 << std::setw(10) << 0.0;

            for (int i = 0; i < n; ++i) {
                file << std::setw(10) << alfa[i][j];
            }
            file << std::endl;
        }
        file.close();

        // t_clad.dat
        file.open(output_dir + "//t_clad.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(2);

            file << std::setw(10)<< 0.0
                 << std::setw(10) << 0.0;

            for (int i = 0; i < n; ++i) {
                file << std::setw(10) << t_clad[i][j][0];
            }
            file << std::endl;
        }
        file.close();

        // t_fuel.dat
        file.open(output_dir + "//t_fuel.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(2);

            file << std::setw(10) << 0.0
                 << std::setw(10) << 0.0;

            for (int i = 0; i < n; ++i) {
                file.precision(4);
                file << std::setw(12) << t_fuel[i][j][0];
            }
            file << std::endl;
        }
        file.close();

        // n_rod.dat
        file.open(output_dir + "//n_rod.dat");
        for (int ii = 0; ii < n_rod + 2; ++ii) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(10) << 1000 * bundle[ii][0] << std::endl;
        }
        file.close();

        // Q.dat
        file.open(output_dir + "//Q.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(15) << xx[j]
                 << std::setw(15) << yy[j];

            for (int i = 0; i < n; ++i) {
                file << std::setw(15) << Q_neutron[i][j][0];
            }
            file << std::endl;
        }
        file.close();

        // t_r77.dat
        file.open(output_dir + "//t_r77.dat");
        for (int ii = 0; ii < n; ++ii) {
            file.setf(std::ios::fixed);
            file.precision(4);

            for (int l = 0; l < n_rod + 2; ++l) {
                file << std::setw(12) << t_rod[l][ii][76][0];
            }
            file << std::endl;
        }
        file.close();

        // t_r_centr.dat
        file.open(output_dir + "//t_r_centr.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);

            file << std::setw(12) << 0.0
                << std::setw(12) << 0.0;

            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(12) << t_rod[0][ii][j][0];
            }
            file << std::endl;
        }
        file.close();


        double tmax = 300.0;
        double t_tepl_max = 287.0;
        double t_obol_max = 287.0;
        int iimax = 0, jjmax = 0;
        int iiimax = 0, jjjmax = 0;
        int iiiimax = 0, jjjjmax = 0;

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

        std::ofstream XXX_file(output_dir + "//XXX.dat");
        XXX_file << "Max T topl: " << tmax << "\n"
                 << "V kassete # " << jjmax + 1 << "\n"
                 << "V el-te # " << iimax + 1<< "\n"
                 << "Max T tepl: " << t_tepl_max << "\n"
                 << "V kassete # " << jjjmax + 1 << "\n"
                 << "V el-te # " << iiimax + 1<< "\n"
                 << "Max T obol: " << t_obol_max << "\n"
                 << "V kassete # " << jjjjmax + 1 << "\n"
                 << "V el-te # " << iiiimax + 1 << std::endl;
        XXX_file.close();

        // 1. t_rXXX.dat
        file.open(output_dir + "//t_rXXX.dat");
        for (int ii = 0; ii < n; ++ii) {
            for (int l = 0; l < n_rod + 2; ++l) {
                file.setf(std::ios::fixed);
                file.precision(4);

                file << std::setw(12) << t_rod[l][ii][jjmax][0];
            }
            file << std::endl;
        }
        file.close();

        // 2. t_r_obl_vnutr.dat
        file.open(output_dir + "//t_r_obl_vnutr.dat");
        for (int j = 0; j < mf; ++j) {

            file.setf(std::ios::fixed);
            file.precision(4);
            file << std::setw(12) << 0.0
                 << std::setw(12) << 0.0;

            for (int ii = 0; ii < n; ++ii) {

                file << std::setw(12) << t_rod[n_rod][ii][j][0];
            }
            file << "\n";
        }
        file.close();

        // 3. t_r_obl_naruj.dat
        file.open(output_dir + "//t_r_obl_naruj.dat");
        for (int j = 0; j < mf; ++j) {

            file.setf(std::ios::fixed);
            file.precision(4);
            file << std::setw(12) << 0.0
                 << std::setw(12) << 0.0;

            for (int ii = 0; ii < n; ++ii) {

                file << std::setw(12) << t_rod[n_rod + 1][ii][j][0];
            }
            file << std::endl;
        }
        file.close();

        // 4. koeff_z.dat
        file.open(output_dir + "//koeff_z.dat");
        for (int j = 0; j < mf; ++j) {
            file.setf(std::ios::fixed);
            file.precision(4);
            file << std::setw(12) << 0.0
                 << std::setw(12) << 0.0;

            for (int i = 0; i < n; ++i) {
                file << std::setw(10)<< effK_z[i][j];
            }
            file << std::endl;
        }
        file.close();

        // 5. Cp.dat
        file.open(output_dir + "//Cp.dat");
        file.setf(std::ios::fixed);
        file.precision(4);

        file << std::setw(11) << Cp_input 
             << std::setw(11) << Cp_output << std::endl;

        for (int j = 0; j < mf; ++j) {
            file << std::setw(10) << xx[j]
                 << std::setw(10) << yy[j];
            for (int ii = 0; ii < n; ++ii) {
                file << std::setw(11) << C_p[ii][j];
            }
            file << std::endl;
        }
        file.close();


        double G = 0;
        double G1 = 0;

        file.open(output_dir + "//GQT.dat");
        for (int j = 0; j < mf; ++j) {
            G += V_z[0][j] * fz; 
            G1 += V_z[n][j] * fz;
        }

        file.precision(4);
        file << "G = " << G << "\n"
             << "G1 = " << G1 << "\n"
             << "Q = " << Q << "\n"
             << "t_CoreOutput = " << t_CoreOutput << "\n"
             << "t_CoreInput = " << t_CoreInput << "\n"
             << "(t_CoreOutput - t_CoreInput) = " << (t_CoreOutput - t_CoreInput) << "\n"
             << "ro_input[76] = " << ro_input[76] << "\n"
             << "ro_output = " << ro_output << std::endl;
        file.close();
    
        if (icall < 1) {
            std::cout << "****//no calculation is required since the imbalance is less than the specified error//****" << std::endl;

            end = std::chrono::high_resolution_clock::now();

            std::cout << "\n";
            std::cout << "Calculation time = " << std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << " min" << std::endl;
            std::cout << "Steady state time = " << time_ << " sec" << std::endl;

            std::cin.get();
            std::exit(EXIT_FAILURE);
        }
    }

    Writer(std::ofstream("data.dat", std::ios::binary)).write(p, V_z, V_n, h_f, t_f, t_rod, t_fuel, t_clad, time_);


    end = std::chrono::high_resolution_clock::now();

    std::cout << "\n";
    std::cout << "Steady state time = " << time_ << " sec" << std::endl;
    std::cout << "Calculation time = " << std::chrono::duration_cast<std::chrono::minutes>(end - start).count() << " min" << std::endl;
    std::cout << "\n" << "End of program" << std::endl;
    std::cin.get();

    return 0;
}