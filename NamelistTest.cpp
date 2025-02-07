#include "headers/ThechycoGlobalVar.hpp"
#include "headers/NamelistReader.hpp"
#include "headers/Thechyco.hpp"
#include <vector>
#include <fstream>
#include <iomanip>

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


int main()
{
    NLReader::NamelistReader nlr("THEHYCO.INI"); // All file reading is done here

    nlr.use_namelist("RodList");

    double s_mesh = nlr.get<double>("s_mesh", 1);

    d_mesh = nlr.get_required<double>("d_mesh");

    R_contact = nlr.get_required<double>("R_contact");

    std::vector<double> row;

    for (int i = 0; i<4; i++){
        row.push_back(nlr.get<double>("tmp",1.0,i));
        // row[i] = nlr.get<double>("tmp",1.0,i);
    }

    tmp[0] = row;

    std::cout << s_mesh <<" \n"<< d_mesh << " \n" << R_contact << std::endl;
    for (auto var : tmp[0]) 
    {
        std::cout << var << " ";
    }

    std::cout << "\n";

    nlr.use_namelist("HEATandHYDROlist");

    for (int i = 0; i<9; i++)
    {
        blockade[i] = nlr.get<int>("blockade",1.0,i);
        std::cout << blockade[i] << " ";
    }

}