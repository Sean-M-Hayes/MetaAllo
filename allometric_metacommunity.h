/**************************************************************************
**    Copyright 2017 Sean M. Hayes
**    This file is part of MetaAllo.
**
**    MetaAllo is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.

**    MetaAllo is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.

**    You should have received a copy of the GNU General Public License
**    along with MetaAllo.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#ifndef ALLOMETRIC_METACOMMUNITY
#define ALLOMETRIC_METACOMMUNITY

#include "parameter_generator_functions.h"
#include "measure_time_series.h"

// <string> & "food_web_methods.h" from "parameter_generator_functions.h"
// <boost/multi_array.hpp> & <vector> from "food_web_methods.h"
// boost_matrix definition from "food_web_methods.h";

std::vector<double> vSequence(double initial, double final, double by);
std::vector<int> vSequence(int initial, int final, int by);

struct Eqn_Parameters{

    int nSpecies;
    int nPatches;

    double Extinction_Threshold;

    int n_Community_Parameters;
    int n_Species_Parameters;
    int n_Patches_Parameters;

    std::vector<std::string> Community_Data_Labels;
    std::vector<std::string> Species_Data_Labels;
    std::vector<std::string> Patches_Data_Labels;

    std::vector<double> Functional_Response_Shape;
    std::vector<double> Half_Saturation;
    std::vector<double> y;
    std::vector<double> ar;
    std::vector<double> ax;
    std::vector<double> e;

    std::vector<double> M;
    std::vector<double> m;

    double Allometric_Scaled_H;
    double Allometric_Scaled_AMI;
    double Allometric_Scaled_Hc;

    double Community_Min_Consumer_Resource_Size_Ratio;
    double Community_Max_Consumer_Resource_Size_Ratio;
    double Community_Avg_Consumer_Resource_Size_Ratio;

    std::vector<double> Species_Min_Consumer_Resource_Size_Ratio;
    std::vector<double> Species_Max_Consumer_Resource_Size_Ratio;
    std::vector<double> Species_Avg_Consumer_Resource_Size_Ratio;

    int Food_Web_Unique;
    Food_Web_Metadata Food_Web_Data;

    int Spatial_Structure_Unique;
    Spatial_Structure_Metadata Spatial_Structure_Data;

    boost_matrix k;

    std::vector<double> Patch_Size;

    boost_matrix Initial_Abundance;

    Eqn_Parameters();

    Eqn_Parameters(Parameter_Generators Input_Generators);

    //void Set_to_Single_Community();

    int n_Community_Data() const;
    int n_Species_Data() const;
    int n_Patches_Data() const;

    double Get_Community_Data(int Index) const;
    std::vector<double> Get_Species_Data(int Index) const;
    std::vector<double> Get_Patches_Data(int Index) const;
};

//struct Simulation_Properties {

//    int nThreads;
//    int nIterations;
//    int nFood_Webs;
//    int nSpatial_Structures;
//    double Simulation_Resolution;
//    std::vector<double> Extinction_Threshold;

//    std::string Write_Path;
//    int Time_Series_Write;
//    int Limit_Cycle_Write;
//    int Species_Data_Write;
//    bool Write_Each_Set;
//    bool Test_against_no_Structure;
//    bool Write_Measure_Diagnostics;

//    double Minimum_Simulation_Duration;
//    double Maximum_Simulation_Duration;
//    double Equilibrium_Test_Interval_Duration;
//    double Time_between_Equilibrium_Tests;

//    double Phase_Lock_Threshold;

//    Cycle_Detection_Parameters Cycle_Params;
//};

struct Sim_Output {

    int ID;
    Eqn_Parameters Eqn_In;
    Dynamics Dyn_Out;

    Sim_Output();

    Sim_Output(int ID, Eqn_Parameters new_Eqn, Dynamics new_Dyn);

};

struct Sim_Set_Output {

    std::vector<Sim_Output> Sim_Out;

    std::vector<std::vector<std::string>> Set_Gen_Parameter_Labels;
    std::vector<std::vector<double>> Set_Gen_Parameter_Values;

    std::vector<int> Feasible_Simulations;
    std::vector<int> Non_Feasible_Simulations;
    std::vector<int> Single_Patch_Simulations;
    std::vector<int> Multi_Patch_Simulations;

    Sim_Set_Output(){}
    Sim_Set_Output(const Parameter_Generators &Set_Gens);
};

std::vector<Food_Web_Metadata> Generate_n_Unique_Food_Webs(int n_Food_Webs, const Parameter_Generators &Input_Gen);

std::vector<Spatial_Structure_Metadata> Generate_n_Unique_Spatial_Structures(int n_Spatial_Structures, const Parameter_Generators &Input_Gen);

class Meta_Allo_Simulation{

    public:

        //Store Parameters

        int counter;

        int iSet;

        int nSimulations;

        Parameter_Generators Par_In;

        Simulation_Properties Sim_In;

        Sim_Output Simulate (Eqn_Parameters x_In, const Simulation_Properties Sim_Props);

        Meta_Allo_Simulation(Parameter_Generators x_In, Simulation_Properties Simulation_Properties_Input);

        void Run_Careful_Simulation();
};


#endif // ALLOMETRIC_METACOMMUNITY

