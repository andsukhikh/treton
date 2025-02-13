#include <vector>
#include <string>

template <typename T, typename Head, typename... Tail>
void resize_vector(std::vector<T>& vec, Head head, Tail... tail) {
    vec.resize(head);
    for (auto&& var : vec) {
        if constexpr (!std::is_scalar_v<std::remove_cvref_t<decltype(var)>>)
        {
            resize_vector(var, tail...);
        }
    }
}

namespace GlobalVariables {

    double Gravity = 9.8154;    //ускорение свободного падения
    int mf;                     //число твс
    int mV_n;                   //число взаимодействующих граней
    int n;                      //число разбиений по высоте
    int n_rod;                  //число разбиений по радиусу твэлаod
    int type;                   //число типов твэловe
    int nbf;                    //число граней в элементе

    namespace Power_per_TBC {
        std::vector<std::vector<std::vector<double>>> Q_neutron;

        void initialisation() {
            resize_vector(Q_neutron, n, mf, type);
        }
    }

    namespace temperature {
        std::vector<std::vector<double>> t_f;
        std::vector<std::vector<std::vector<double>>> t_fuel;
        std::vector<std::vector<std::vector<double>>> t_clad;

        void initialisation() {
            resize_vector(t_f, n, mf);
            resize_vector(t_fuel, n, mf, type);
            resize_vector(t_clad, n, mf, type);
        }
    }

    namespace rod_geometry {
        double x_mesh;
        double s_mesh;
        double d_mesh;
        std::vector<std::vector<double>> tmp;
        std::vector<double> P_rod;
        std::vector<double> a_fuel;
        std::vector<double> a_clad;
        std::vector<std::vector<double>> area;
        std::vector<std::vector<double>> geo_left;
        std::vector<std::vector<double>> geo_right;

        void initialisation() {
            resize_vector(tmp, 4, type);
            resize_vector(P_rod, type);
            resize_vector(a_fuel, type);
            resize_vector(a_clad, type);
            resize_vector(area, n_rod + 2, type);
            resize_vector(geo_left, n_rod + 2, type);
            resize_vector(geo_right, n_rod + 2, type);
        }
    }

    namespace rod_temperature {
        std::vector<std::vector<std::vector<std::vector<double>>>> t_rod;

        void initialisation() {
            resize_vector(t_rod, n_rod + 2, n, mf, type);
        }
    }

    namespace Rod_property {
        double clad_l;
        double clad_rc;
        std::vector<std::vector<std::vector<double>>> fuel_l;
        std::vector<std::vector<std::vector<double>>> fuel_rc;

        void initialisation() {
            resize_vector(fuel_l, n, mf, type);
            resize_vector(fuel_rc, n, mf, type);
        }
    }

    namespace bundles {
        std::vector<std::vector<double>> bundle;

        void initialisation() {
            resize_vector(bundle, n_rod + 2, type);
        }
    }

    namespace contact_resistance {
        double R_contact;
    }

    namespace cv_distribution {
        std::vector<std::vector<double>> cvd;

        void initialisation() {
            resize_vector(cvd, n_rod + 3, type);
        }
    }

    namespace Density {
        double ro_output;
        std::vector<double> ro_input;
        std::vector<std::vector<double>> ro;

        void initialisation() {
            resize_vector(ro_input, mV_n);
            resize_vector(ro, n, mf);
        }
    }

    namespace pressure {
        double p_output;
        double p_input;
        std::vector<std::vector<double>> p;

        void initialisation() {
            resize_vector(p, n, mf);
        }
    }

    namespace capacity {
        double Cp_input;
        double Cp_output;
        std::vector<std::vector<double>> C_p;
        std::vector<std::vector<double>> OLD_C_p;

        void initialisation() {
            resize_vector(C_p, n, mf);
            resize_vector(OLD_C_p, n, mf);
        }
    }

    namespace velocity {
        std::vector<std::vector<double>> V_z;
        std::vector<std::vector<double>> V_n;
        std::vector<std::vector<std::vector<double>>> V_full;

        void initialisation() {
            resize_vector(V_z, n + 1, mf);
            resize_vector(V_n, n + 1, mV_n);
            resize_vector(V_full, 3, n + 1, mf);
        }
    }

    namespace nodes {
        std::vector<std::vector<int>> crd;

        void initialisation() {
            resize_vector(crd, 2, mf);
        }
    }

    namespace bondaries {
        std::vector<std::vector<double>> bonds;
        std::vector<std::vector<double>> onds;

        void initialisation() {
            resize_vector(bonds, nbf, mf);
            resize_vector(onds, nbf, mf);
        }
    }

    namespace blockages {
        double BlockadePorousity;
        std::vector<int> blockade;

        void initialisation() {
            resize_vector(blockade, mf);
        }
    }

    namespace errors {
        double PVTerror;
        double Disbalance;
    }

    namespace boundary_temperature {
        double h_HeatExchangerOutput;
        double h_CoreInput;
        double h_CoreOutput;
        double t_HeatExchangerOutput;
        double t_CoreInput;
        double t_CoreOutput;
        std::vector<double> h_HeatExchangerOutput_new;

        void initialisation() {
            resize_vector(h_HeatExchangerOutput_new, mf);
        }
    }

    namespace core_geometry {
        double D_tube;
        double dz;
        double dr;
        double fz;
        double fr;
        double vf;
        double fr_vf;
        double fz_vf;
        std::vector<int> n_RodsInTBC;

        void initialisation() {
            resize_vector(n_RodsInTBC, type);
        }
    }

    namespace loop {
        int iterations;
    }

    namespace info {
        std::string HydroStr;
        std::string HeatStr;
        std::string RodStr;
        std::string HeatStr1;
    }

    namespace mat {
        std::vector<double> Mat_A;
        std::vector<double> Vect_X;
        std::vector<double> Vect_X1;
        std::vector<double> Vect_B;
        std::vector<int> NC;
        std::vector<int> NE;

        void initialisation() {
            resize_vector(Mat_A, ((mf + mV_n) * n - mf) * 2 + n * mf + 1);
            resize_vector(Vect_X, n * mf + 1);
            resize_vector(Vect_X1, n * mf + 1);
            resize_vector(Vect_B, n * mf + 1);
            resize_vector(NC, (n * mf) + 1 + 1);
            resize_vector(NE, ((mf + mV_n) * n - mf) * 2 + 1);
        }
    }

    namespace supporting_var {
        std::vector<std::vector<std::vector<std::vector<double>>>> OLDt_rod;
        std::vector<std::vector<double>> OLDt_f, oldV_z, oldV_n, OLD_ro, OLDh_f, h_f;

        std::vector<std::vector<double>> effL, effM, effK_r, effK_z, alfa;

        std::vector<double> aa, bb, cc, dd;
        std::vector<double> a, b, c, d;
        std::vector<double> aPV, bPV, cPV, dPV, ePV, fPV;

        std::vector<std::vector<double>> source, dp_dz, dp_dn;
        std::vector<double> V_nMap;

        std::vector<std::vector<double>> bes;
        std::vector<double> xx;
        std::vector<double> yy;
        std::vector<double> z;

        void initialisation() {
            resize_vector(OLDt_rod, n_rod + 2, n, mf, type);

            resize_vector(OLDt_f, n, mf);
            resize_vector(oldV_z, n + 1, mf);
            resize_vector(oldV_n, n, mV_n);
            resize_vector(OLD_ro, n, mf);
            resize_vector(OLDh_f, n, mf);
            resize_vector(h_f, n, mf);

            resize_vector(effL, n, mf);
            resize_vector(effM, n + 1, mf);
            resize_vector(effK_r, n, mV_n);
            resize_vector(effK_z, n + 1, mf);
            resize_vector(alfa, n, mf);

            resize_vector(aa, n);  resize_vector(a, n_rod + 2);
            resize_vector(bb, n);  resize_vector(b, n_rod + 2);
            resize_vector(cc, n);  resize_vector(c, n_rod + 2);
            resize_vector(dd, n);  resize_vector(d, n_rod + 2);

            resize_vector(aPV, n + 1);
            resize_vector(bPV, n + 1);
            resize_vector(cPV, n + 1);
            resize_vector(dPV, n + 1);
            resize_vector(ePV, n + 1);
            resize_vector(fPV, n + 1);

            resize_vector(source, n, mf);
            resize_vector(dp_dz, n + 1, mf);
            resize_vector(dp_dn, n, mV_n);
            resize_vector(V_nMap, mV_n);

            resize_vector(bes, mf, n);
            resize_vector(xx, mf);
            resize_vector(yy, mf);
            resize_vector(z, n + 1);


        }
    }

    namespace transosed_matrix_elem {
        std::vector<double> transposed_OLDt_rod;
        std::vector<double> transposed_t_rod;
        std::vector<double> transposed_geo_left;
        std::vector<double> transposed_geo_right;
        std::vector<double> transposed_bundle;

        void initialisation() {
            resize_vector(transposed_OLDt_rod, n_rod + 2);
            resize_vector(transposed_t_rod, n_rod + 2);
            resize_vector(transposed_geo_left, n_rod + 2);
            resize_vector(transposed_geo_right, n_rod + 2);
            resize_vector(transposed_bundle, n_rod + 2);
        }
    }

    void Initialisation() {
        Power_per_TBC::             initialisation();
        temperature::               initialisation();
        rod_geometry::              initialisation();
        rod_temperature::           initialisation();
        Rod_property::              initialisation();
        bundles::                   initialisation();
        cv_distribution::           initialisation();
        Density::                   initialisation();
        pressure::                  initialisation();
        capacity::                  initialisation();
        velocity::                  initialisation();
        nodes::                     initialisation();
        bondaries::                 initialisation();
        blockages::                 initialisation();
        boundary_temperature::      initialisation();
        core_geometry::             initialisation();
        mat::                       initialisation();
        transosed_matrix_elem::     initialisation();
        supporting_var::            initialisation();
    }
}