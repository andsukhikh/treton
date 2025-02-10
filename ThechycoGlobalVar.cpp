#include <vector>
#include <string>

constexpr auto mf = 163; //число твс
constexpr auto mV_n = 534; //число взаимодействующих граней
constexpr auto n = 30; // число разбиений по высоте
constexpr auto n_rod = 10;// число разбиений по радиусу твэла
constexpr auto type = 1;// число типов твэлов
constexpr auto nbf = 6;//число граней в элементе
constexpr auto Pi = 3.1415926;
constexpr auto Gravity = 9.8;
constexpr auto c_p1 = 1269.0;//уельная теплоемкость
constexpr auto dt = 0.005;


namespace Power_per_TBC {
    std::vector<std::vector<std::vector<double>>> Q_neutron(
        n, 
        std::vector<std::vector<double>>(
            mf, 
            std::vector<double>(type) 
        ) 
    ); 
}


namespace temperature {
    std::vector<std::vector<double>> t_f(n, std::vector<double>(mf)); 

    std::vector<std::vector<std::vector<double>>> t_fuel(
        n, 
        std::vector<std::vector<double>>(
            mf, 
            std::vector<double>(type) 
        ) 
    ); 

    std::vector<std::vector<std::vector<double>>> t_clad(
        n, 
        std::vector<std::vector<double>>(
            mf, 
            std::vector<double>(type)
        )
    ); 
}


namespace rod_geometry {
    double x_mesh;
    double d_mesh;
    std::vector<std::vector<double>> tmp(4, std::vector<double>(type) );
    std::vector<double> P_rod(type);
    std::vector<double> a_fuel(type);
    std::vector<double> a_clad(type);
    std::vector<std::vector<double>> area(n_rod + 2, std::vector<double>(type));
    std::vector<std::vector<double>> geo_left(n_rod + 2, std::vector<double>(type));
    std::vector<std::vector<double>> geo_right(n_rod + 2, std::vector<double>(type));
}


namespace rod_temperature {
    std::vector<std::vector<std::vector<std::vector<double>>>> t_rod(
        n_rod + 2, 
        std::vector<std::vector<std::vector<double>> >(
            n, 
            std::vector<std::vector<double>>(
                mf, 
                std::vector<double>(type) 
            ) 
        ) 
    ); 
}


namespace rod_property1{
    double clad_l, clad_rc;
    std::vector<std::vector<std::vector<double>>> fuel_l(n, std::vector<std::vector<double>>(mf, std::vector<double>(type) ) ); 
    std::vector<std::vector<std::vector<double>>> fuel_rc(n, std::vector<std::vector<double>> (mf, std::vector<double>(type) ) ); 
}

namespace bundles {
    std::vector<std::vector<double>> bundle(n_rod + 2, std::vector<double> (type));
}


namespace contact_resistance {
    double R_contact;
}


namespace cv_distribution {
    std::vector<std::vector<double>> cvd(n_rod + 3, std::vector<double> (type));
}


namespace density1 {
    double ro_output; 
    std::vector<double> ro_input(mV_n);
    std::vector<std::vector<double>> ro(n, std::vector<double>(mf)); 
}


namespace pressure {
    double p_output;
    double p_input;
    std::vector<std::vector<double>> p(n, std::vector<double>(mf)); 
}


namespace capacity {
    double Cp_input, Cp_output;
    std::vector<std::vector<double>> C_p(n, std::vector<double>(mf)); 
    std::vector<std::vector<double>> OLD_C_p(n, std::vector<double>(mf)); 
}

namespace velocity {
    std::vector<std::vector<double>> V_z(n + 1, std::vector<double>(mf));
    std::vector<std::vector<double>> V_n(n + 1, std::vector<double>(mV_n));
    std::vector<std::vector<std::vector<double>>> V_full(3, std::vector<std::vector<double>>(n + 1, std::vector<double> (mf))); 
}

namespace nodes {
    std::vector<std::vector<int>> crd(2, std::vector<int>(mf));
}

namespace bondaries {
    std::vector<std::vector<double>> bonds(nbf, std::vector<double>(mf)); 
    std::vector<std::vector<double>> onds(nbf, std::vector<double>(mf)); 
}

namespace blockages {
    double BlockadePorousity; 
    std::vector<int> blockade(mf);
}

namespace errors {
    double PVTerror, Disbalance; 
}

namespace boundary_temperature {
    double h_HeatExchangerOutput; 
    double h_CoreInput; 
    double h_CoreOutput; 
    double t_HeatExchangerOutput; 
    double t_CoreInput; 
    double t_CoreOutput; 
    std::vector<double> h_HeatExchangerOutput_new(mf); 
}

namespace core_geometry {
    double D_tube, dz, dr, fz, fr, vf, fr_vf, fz_vf; 
    std::vector<int> n_RodsInTBC(type);
}

namespace loop {
    int iterations; 
}

namespace info {
    std::string HydroStr, HeatStr, RodStr, HeatStr1; 
}


    std::vector<std::vector<std::vector<std::vector<double>>>> OLDt_rod(
        n_rod + 2,                                              
        std::vector<std::vector<std::vector<double>>>(
            n,                                                   
            std::vector<std::vector<double>>(
                mf,                                                
                std::vector<double>(type)                  
            )
        )
    );

    std::vector<double> aa(n), bb(n), cc(n), dd(n);
    std::vector<std::vector<double>> effL(n, std::vector<double>(mf));
    std::vector<std::vector<double>> OLDt_f(n, std::vector<double>(mf));
    std::vector<std::vector<double>> alfa(n, std::vector<double>(mf));                      
    std::vector<double> a(n_rod + 2), b(n_rod + 2), c(n_rod + 2), d(n_rod + 2);
    std::vector<double> aPV(n + 1), bPV(n + 1), cPV(n + 1), dPV(n + 1), ePV(n + 1), fPV(n + 1);
    std::vector<std::vector<double>> effM(n + 1, std::vector<double>(mf));
    std::vector<std::vector<double>> effK_r(n, std::vector<double>(mV_n));
    std::vector<std::vector<double>> effK_z(n + 1, std::vector<double>(mf));
    std::vector<std::vector<double>> oldV_z(n + 1, std::vector<double>(mf));
    std::vector<std::vector<double>> oldV_n(n, std::vector<double>(mV_n));
    std::vector<std::vector<double>> source(n, std::vector<double>(mf));
    std::vector<std::vector<double>> dp_dz(n + 1, std::vector<double>(mf));
    std::vector<std::vector<double>> dp_dn(n, std::vector<double>(mV_n));
    std::vector<std::vector<double>> OLD_ro(n, std::vector<double>(mf));
    std::vector<std::vector<double>> OLDh_f(n, std::vector<double>(mf));
    std::vector<std::vector<double>> h_f(n, std::vector<double>(mf));
    std::vector<double> V_nMap(mV_n);

    std::vector<double> transposed_OLDt_rod(n_rod + 2);
    std::vector<double> transposed_t_rod(n_rod + 2);
    std::vector<double> transposed_geo_left(n_rod + 2);
    std::vector<double> transposed_geo_right(n_rod + 2);
    std::vector<double> transposed_bundle(n_rod + 2);



namespace mat {
    std::vector<double> Mat_A(((mf + mV_n) * n - mf) * 2 + n * mf + 1);
    std::vector<double> Vect_X(n * mf + 1);
    std::vector<double> Vect_X1(n * mf + 1);
    std::vector<double> Vect_B(n * mf + 1);
    std::vector<int> NC((n * mf) + 1 + 1);
    std::vector<int> NE(((mf + mV_n) * n - mf) * 2 + 1);
}
