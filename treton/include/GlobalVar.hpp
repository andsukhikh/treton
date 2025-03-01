#ifndef THECHYCO_GLOBAL_VAR_H
#define THECHYCO_GLOBAL_VAR_H

#include <vector>
#include <string>

template <typename T, typename Head, typename... Tail>
extern void resize_vector(std::vector<T>& vec, Head head, Tail... tail);

namespace GlobalVariables {

    extern double Gravity;
    extern int mf;
    extern int mV_n;
    extern int n;
    extern int n_rod;
    extern int type;
    extern int nbf;

    namespace Power_per_TBC {
        extern std::vector<std::vector<std::vector<double>>> Q_neutron;

        void initialisation();
    }

    namespace temperature {
        extern std::vector<std::vector<double>> t_f;
        extern std::vector<std::vector<std::vector<double>>> t_fuel;
        extern std::vector<std::vector<std::vector<double>>> t_clad;

        void initialisation();
    }

    namespace rod_geometry {
        extern double x_mesh;
        extern double s_mesh;
        extern double d_mesh;
        extern std::vector<std::vector<double>> tmp;
        extern std::vector<double> P_rod;
        extern std::vector<double> a_fuel;
        extern std::vector<double> a_clad;
        extern std::vector<std::vector<double>> area;
        extern std::vector<std::vector<double>> geo_left;
        extern std::vector<std::vector<double>> geo_right;

        void initialisation();
    }

    namespace rod_temperature {
        extern std::vector<std::vector<std::vector<std::vector<double>>>> t_rod;

        void initialisation();
    }

    namespace Rod_property {
        extern double clad_l;
        extern double clad_rc;
        extern std::vector<std::vector<std::vector<double>>> fuel_l;
        extern std::vector<std::vector<std::vector<double>>> fuel_rc;

        void initialisation();
    }

    namespace bundles {
        extern std::vector<std::vector<double>> bundle;

        void initialisation();
    }

    namespace contact_resistance {
        extern double R_contact;
    }

    namespace cv_distribution {
        extern std::vector<std::vector<double>> cvd;

        void initialisation();
    }

    namespace Density {
        extern double ro_output;
        extern std::vector<double> ro_input;
        extern std::vector<std::vector<double>> ro;

        void initialisation();
    }

    namespace pressure {
        extern double p_output;
        extern double p_input;
        extern std::vector<std::vector<double>> p;

        void initialisation();
    }

    namespace capacity {
        extern double Cp_input;
        extern double Cp_output;
        extern std::vector<std::vector<double>> C_p;
        extern std::vector<std::vector<double>> OLD_C_p;

        void initialisation();
    }

    namespace velocity {
        extern std::vector<std::vector<double>> V_z;
        extern std::vector<std::vector<double>> V_n;
        extern std::vector<std::vector<std::vector<double>>> V_full;

        void initialisation();
    }

    namespace nodes {
        extern std::vector<std::vector<int>> crd;

        void initialisation();
    }

    namespace bondaries {
        extern std::vector<std::vector<double>> bonds;
        extern std::vector<std::vector<double>> onds;

        void initialisation();
    }

    namespace blockages {
        extern double BlockadePorousity;
        extern std::vector<int> blockade;

        void initialisation();
    }

    namespace errors {
        extern double PVTerror;
        extern double Disbalance;
    }

    namespace boundary_temperature {
        extern double h_HeatExchangerOutput;
        extern double h_CoreInput;
        extern double h_CoreOutput;
        extern double t_CoreInput;
        extern double t_CoreOutput;
        extern std::vector<double> h_HeatExchangerOutput_new;

        void initialisation();
    }

    namespace core_geometry {
        extern double D_tube;
        extern double dz;
        extern double Height;
        extern double dr;
        extern double fz;
        extern double fr;
        extern double vf;
        extern double fr_vf;
        extern double fz_vf;

        extern std::vector<int> n_RodsInTBC;

        void initialisation();
    }

    namespace loop {
        extern int iterations;
    }

    namespace info {
        extern std::string HydroStr;
        extern std::string HeatStr;
        extern std::string RodStr;
        extern std::string HeatStr1;
        extern std::string coolantName;
        extern std::string input_dir;
        extern std::string output_dir;
    }

    namespace mat {
        extern std::vector<double> Mat_A;
        extern std::vector<double> Vect_X;
        extern std::vector<double> Vect_X1;
        extern std::vector<double> Vect_B;
        extern std::vector<int> NC;
        extern std::vector<int> NE;

        void initialisation();
    }

    namespace supporting_var {
        extern std::vector<std::vector<std::vector<std::vector<double>>>> OLDt_rod;
        extern std::vector<std::vector<double>> OLDt_f, oldV_z, oldV_n, OLD_ro, OLDh_f, h_f;

        extern std::vector<std::vector<double>> effL, effM, effK_r, effK_z, alfa;

        extern std::vector<double> aa, bb, cc, dd;
        extern std::vector<double> a, b, c, d;
        extern std::vector<double> aPV, bPV, cPV, dPV, ePV, fPV;

        extern std::vector<std::vector<double>> source, dp_dz, dp_dn;
        extern std::vector<double> V_nMap;

        extern std::vector<std::vector<double>> bes;
        extern std::vector<double> xx;
        extern std::vector<double> yy;
        extern std::vector<double> z;

        void initialisation();
    }

    namespace transosed_matrix_elem {
        extern std::vector<double> transposed_OLDt_rod;
        extern std::vector<double> transposed_t_rod;
        extern std::vector<double> transposed_geo_left;
        extern std::vector<double> transposed_geo_right;
        extern std::vector<double> transposed_bundle;

        void initialisation();
    }

    using namespace capacity;
    using namespace loop;
    using namespace temperature;
    using namespace Density;
    using namespace Power_per_TBC;
    using namespace temperature;
    using namespace rod_geometry;
    using namespace rod_temperature;
    using namespace Rod_property;
    using namespace bundles;
    using namespace contact_resistance;
    using namespace cv_distribution;
    using namespace pressure;
    using namespace velocity;
    using namespace nodes;
    using namespace bondaries;
    using namespace blockages;
    using namespace errors;
    using namespace boundary_temperature;
    using namespace core_geometry;
    using namespace info;
    using namespace mat;
    using namespace transosed_matrix_elem;
    using namespace supporting_var;

    void Resizing();
}

using namespace GlobalVariables;

#endif 
