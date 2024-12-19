#ifndef THECHYCO_GLOBAL_VAR_H
#define THECHYCO_GLOBAL_VAR_H

#include <vector>
#include <string>

constexpr auto mf = 163;
constexpr auto mV_n = 534;//должно быть 534, но при нем не работает 
constexpr auto n = 30;
constexpr auto n_rod = 10;
constexpr auto type = 1;
constexpr auto nbf = 6;
constexpr auto Pi = 3.1415926;
constexpr auto Gravity = 9.8;
constexpr auto c_p1 = 4800.0;
constexpr auto dt = 0.005;

namespace Power_per_TBC {
    extern std::vector<std::vector<std::vector<double>>> Q_neutron;
}

namespace temperature {
    extern std::vector<std::vector<double>> t_f;
    extern std::vector<std::vector<std::vector<double>>> t_fuel;
    extern std::vector<std::vector<std::vector<double>>> t_clad;
}

namespace rod_geometry {
    extern double x_mesh;
    extern double d_mesh;
    extern std::vector<std::vector<double>> tmp;
    extern std::vector<double> P_rod;
    extern std::vector<double> a_fuel;
    extern std::vector<double> a_clad;
    extern std::vector<std::vector<double>> area;
    extern std::vector<std::vector<double>> geo_left;
    extern std::vector<std::vector<double>> geo_right;
}

namespace rod_temperature {
    extern std::vector<std::vector<std::vector<std::vector<double>>>> t_rod;
}

namespace rod_property1 {
    extern double clad_l;
    extern double clad_rc;
    extern std::vector<std::vector<std::vector<double>>> fuel_l;
    extern std::vector<std::vector<std::vector<double>>> fuel_rc;
}

namespace bundles {
    extern std::vector<std::vector<double>> bundle;
}

namespace contact_resistance {
    extern double R_contact;
}

namespace cv_distribution {
    extern std::vector<std::vector<double>> cvd;
}

namespace density1 {
    extern double ro_output;
    extern std::vector<double> ro_input;
    extern std::vector<std::vector<double>> ro;
}

namespace pressure {
    extern double p_output;
    extern double p_input;
    extern std::vector<std::vector<double>> p;
}

namespace capacity {
    extern double Cp_input;
    extern double Cp_output;
    extern std::vector<std::vector<double>> C_p;
    extern std::vector<std::vector<double>> OLD_C_p;
}

namespace velocity {
    extern std::vector<std::vector<double>> V_z;
    extern std::vector<std::vector<double>> V_n;
    extern std::vector<std::vector<std::vector<double>>> V_full;
}

namespace nodes {
    extern std::vector<std::vector<int>> crd;
}

namespace bondaries {
    extern std::vector<std::vector<double>> bonds;
    extern std::vector<std::vector<double>> onds;
}

namespace blockages {
    extern double BlockadePorousity;
    extern std::vector<int> blockade;
}

namespace errors {
    extern double PVTerror;
    extern double Disbalance;
}

namespace boundary_temperature {
    extern double h_HeatExchangerOutput;
    extern double h_CoreInput;
    extern double h_CoreOutput;
    extern double t_HeatExchangerOutput;
    extern double t_CoreInput;
    extern double t_CoreOutput;
    extern std::vector<double> h_HeatExchangerOutput_new;
}

namespace core_geometry {
    extern double D_tube;
    extern double dz;
    extern double dr;
    extern double fz;
    extern double fr;
    extern double vf;
    extern double fr_vf;
    extern double fz_vf;
    extern std::vector<int> n_RodsInTBC;
}

namespace loop {
    extern int iterations;
}

namespace info {
    extern std::string HydroStr;
    extern std::string HeatStr;
    extern std::string RodStr;
    extern std::string HeatStr1;
}

namespace mat {
    extern std::vector<double> Mat_A;
    extern std::vector<double> Vect_X;
    extern std::vector<double> Vect_X1;
    extern std::vector<double> Vect_B;
    extern std::vector<int> NC;
    extern std::vector<int> NE;
}

extern std::vector<std::vector<std::vector<std::vector<double>>>> OLDt_rod;
extern std::vector<double> aa, bb, cc, dd;
extern std::vector<std::vector<double>> effL, OLDt_f, alfa;
extern std::vector<double> a, b, c, d, aPV, bPV, cPV, dPV, ePV, fPV;
extern std::vector<std::vector<double>> effM, effK_r, effK_z, oldV_z, oldV_n, source, dp_dz, dp_dn, OLD_ro, OLDh_f, h_f;
extern std::vector<double> V_nMap;

using namespace capacity;
    using namespace loop;
    using namespace temperature;
    using namespace density1;
    using namespace Power_per_TBC;
    using namespace temperature;
    using namespace rod_geometry;
    using namespace rod_temperature;
    using namespace rod_property1;
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

#endif 
