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
#include "allometric_metacommunity.h"

// "parameter_generator_functions.h" from "allometric_metacommunity.h"
// <string> & "food_web_methods.h" from "parameter_generator_functions.h"
// <boost/multi_array.hpp> & <vector> from "food_web_methods.h"
// boost_matrix definition from "food_web_methods.h";

//#include "measure_time_series.h"

#include "lsoda_link.h"
#include "write_simulation_data_functions.h"

#include <iostream>
#include <numeric>
#include <fstream>
#include <random>
#include <cmath>
#include <math.h>
#include <ctime>
#define BOOST_THREAD_USE_LIB
#include <boost/thread.hpp>
#include <algorithm>

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <iterator>
#include <cstdlib>

#include <omp.h>

omp_lock_t counter_mutex;

namespace Pass_to_Fortran{

    boost::thread_specific_ptr<Eqn_Parameters> Thread_Local_Parameters_ptr;

    void Allo_Web_Eqn (int *NEQ, double *T, double Y[], double YDOT[])
    {
        //Uses (*Thread_Local_Parameters_ptr): make sure this is defined with the desired parameter values before calling
        std::vector<double> a_x_s(*NEQ,0);
        std::vector<double> Growth(*NEQ,0);
        std::vector<double> Dispersal(*NEQ,0);
        std::vector<double> Consumption(*NEQ,0);
        std::vector<double> Predation(*NEQ,0);
        std::vector<double> Mortality(*NEQ,0);

        for(int idx = 0; idx < (*Thread_Local_Parameters_ptr).nSpecies*(*Thread_Local_Parameters_ptr).nPatches; idx++)
            if(Y[idx]<(*Thread_Local_Parameters_ptr).Extinction_Threshold)
                Y[idx] = 0;


        for(int iP = 0; iP < (*Thread_Local_Parameters_ptr).nPatches; iP++)
        {
            for(int iS1 = 0; iS1 < (*Thread_Local_Parameters_ptr).nSpecies; iS1++)
            {
                Growth[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies] = (*Thread_Local_Parameters_ptr).ar[iS1]
                        *pow((*Thread_Local_Parameters_ptr).M[iS1],-.25)
                        *(*Thread_Local_Parameters_ptr).Food_Web_Data.is_Producer[iS1]
                        *(1-(Y[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies]
                          /((*Thread_Local_Parameters_ptr).k[iS1][iP]*(*Thread_Local_Parameters_ptr).Patch_Size[iP])));

                Mortality[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies] = (*Thread_Local_Parameters_ptr).ax[iS1]
                        *pow((*Thread_Local_Parameters_ptr).M[iS1],-.25)
                        *(*Thread_Local_Parameters_ptr).Food_Web_Data.is_Consumer[iS1];

                for(int iP2 = 0; iP2 <(*Thread_Local_Parameters_ptr).nPatches; iP2++)
                    Dispersal[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies] += (*Thread_Local_Parameters_ptr).m[iS1]
                            *(*Thread_Local_Parameters_ptr).Spatial_Structure_Data.Spatial_Structure[iP][iP2]
                            *Y[iS1+iP2*(*Thread_Local_Parameters_ptr).nSpecies];

                for(int iS2 = 0; iS2<(*Thread_Local_Parameters_ptr).nSpecies; iS2++)
                {
                    if(Y[iS2+iP*(*Thread_Local_Parameters_ptr).nSpecies]>0)
                        a_x_s[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies] += (*Thread_Local_Parameters_ptr).Food_Web_Data.Food_Web[iS2][iS1]
                                * pow(Y[iS2+iP*(*Thread_Local_Parameters_ptr).nSpecies],(*Thread_Local_Parameters_ptr).Functional_Response_Shape[iS1]);
                }

                Consumption[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies] =
                        (*Thread_Local_Parameters_ptr).ax[iS1]
                        *pow((*Thread_Local_Parameters_ptr).M[iS1],-.25)
                        *(*Thread_Local_Parameters_ptr).y[iS1]
                        *a_x_s[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies]
                        /(pow((*Thread_Local_Parameters_ptr).Half_Saturation[iS1]*(*Thread_Local_Parameters_ptr).Patch_Size[iP],(*Thread_Local_Parameters_ptr).Functional_Response_Shape[iS1])
                          +a_x_s[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies]);
            }

            for(int iS1 = 0; iS1 < (*Thread_Local_Parameters_ptr).nSpecies; iS1++)
                for(int iS2 = 0; iS2 < (*Thread_Local_Parameters_ptr).nSpecies; iS2++)
                    Predation[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies] +=
                        ((*Thread_Local_Parameters_ptr).ax[iS2]
                         *pow((*Thread_Local_Parameters_ptr).M[iS2],-.25)
                         *(*Thread_Local_Parameters_ptr).y[iS2]
                         *(*Thread_Local_Parameters_ptr).Food_Web_Data.Food_Web[iS1][iS2]
                         *Y[iS2+iP*(*Thread_Local_Parameters_ptr).nSpecies]
                            *pow(Y[iS1+iP*(*Thread_Local_Parameters_ptr).nSpecies],(*Thread_Local_Parameters_ptr).Functional_Response_Shape[iS1]-1))
                            /((*Thread_Local_Parameters_ptr).e[iS2]
                              *(pow((*Thread_Local_Parameters_ptr).Half_Saturation[iS1]*(*Thread_Local_Parameters_ptr).Patch_Size[iP],(*Thread_Local_Parameters_ptr).Functional_Response_Shape[iS1])
                                +a_x_s[iS2+iP*(*Thread_Local_Parameters_ptr).nSpecies]));
        }

        for(int idx = 0; idx < (*Thread_Local_Parameters_ptr).nSpecies*(*Thread_Local_Parameters_ptr).nPatches; idx++)
        {
            YDOT[idx] = Y[idx]*(Growth[idx]+Consumption[idx]-Predation[idx]-Mortality[idx])+Dispersal[idx];

            if(Y[idx]==0&&YDOT[idx]<(*Thread_Local_Parameters_ptr).Extinction_Threshold)
                YDOT[idx] = 0;
        }
    }

    void JDUM(){}

}

Sim_Output::Sim_Output(int new_ID, Eqn_Parameters new_Eqn, Dynamics new_Dyn):
    ID(new_ID),Eqn_In(new_Eqn),Dyn_Out(new_Dyn){}



Sim_Output Meta_Allo_Simulation::Simulate(Eqn_Parameters x_In, const Simulation_Properties Sim_Props)
{
    //Load into Simulation Function
    Pass_to_Fortran::Thread_Local_Parameters_ptr.reset(new Eqn_Parameters(x_In));

    omp_set_lock(&counter_mutex);
    int ID_temp = counter;
    counter++;
    omp_unset_lock(&counter_mutex);

    std::string Time_Series_Write_String = Sim_Props.Write_Path+std::to_string(static_cast<long long>(ID_temp))+"Timeseries.csv";

    if(Sim_Props.Time_Series_Write>0)
        Start_Time_series(Time_Series_Write_String,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches);

    std::vector<double> Initial_Abundances;

    for(int iP = 0; iP<(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches;iP++)
        for(int iS = 0;iS<(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies;iS++)
            Initial_Abundances.push_back(x_In.Initial_Abundance[iS][iP]);

    double T_Starting = 0;

    Dynamics x_Out;
    x_Out.Reset_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches);

    bool Force_Measure = false;

    for(double i=Sim_Props.Minimum_Simulation_Duration+Sim_Props.Equilibrium_Test_Interval_Duration;i<=Sim_Props.Maximum_Simulation_Duration;i+=Sim_Props.Time_between_Equilibrium_Tests+Sim_Props.Equilibrium_Test_Interval_Duration)
    {
        boost_matrix sim_ana;

        bool Feasible = LSODA_Solve(Initial_Abundances,
                                    &Pass_to_Fortran::Allo_Web_Eqn,
                                    &Pass_to_Fortran::JDUM,
                                    T_Starting,
                                    i-T_Starting,
                                    Sim_Props.Simulation_Resolution,
                                    x_In.Extinction_Threshold,
                                    (*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,
                                    sim_ana);

        if(Sim_Props.Time_Series_Write>0)
        {
            Write_or_Append_Time_series(Time_Series_Write_String,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,sim_ana);
        }

        if(i+Sim_Props.Time_between_Equilibrium_Tests+Sim_Props.Equilibrium_Test_Interval_Duration>Sim_Props.Maximum_Simulation_Duration) //Force exit if Time is max
        {
            Force_Measure = true;
        }

        //Will return true if nonfeasible, stationary, or Force_Measure is true
        bool Return_Test = Measure_Metacommunity_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,
                                                     T_Starting+sim_ana.shape()[0]*Sim_Props.Simulation_Resolution-Sim_Props.Equilibrium_Test_Interval_Duration,
                                                     std::pair<int,int>(sim_ana.shape()[0]-(int)(Sim_Props.Equilibrium_Test_Interval_Duration/Sim_Props.Simulation_Resolution),sim_ana.shape()[0]-1),
                                                     sim_ana,
                                                     Feasible,
                                                     x_In.Extinction_Threshold,
                                                     Sim_Props.Write_Path+std::to_string(static_cast<long long>(ID_temp)),
                                                     x_In.Spatial_Structure_Data.Spatial_Structure,
                                                     Sim_Props,
                                                     Force_Measure,
                                                     x_Out);


        if(Return_Test)
        {
            if(Sim_Props.Time_Series_Write>0)
            {
                if(Feasible)
                {
                    std::string Feas_Write_String = Sim_Props.Write_Path+"Feas_"+std::to_string(static_cast<long long>(ID_temp))+"Timeseries.csv";
                    rename(Time_Series_Write_String.data(),Feas_Write_String.data());
                }
                else if(Sim_Props.Time_Series_Write>1)
                {
                    std::string NonFeas_Write_String = Sim_Props.Write_Path+"NonFeas_"+std::to_string(static_cast<long long>(ID_temp))+"Timeseries.csv";
                    rename(Time_Series_Write_String.data(),NonFeas_Write_String.data());
                }
                else
                    remove(Time_Series_Write_String.data());
            }

            if(Feasible)
            {
                if(Sim_Props.Limit_Cycle_Write>0) // Change to Limit Cycle Write
                {
                    std::string Limit_Cycle_Write_String = Sim_Props.Write_Path+std::to_string(static_cast<long long>(ID_temp));

                    Write_Population_Cycle_Timeseries(Limit_Cycle_Write_String,Sim_Props.Simulation_Resolution,x_Out);
                }
            }


        Pass_to_Fortran::Thread_Local_Parameters_ptr.reset(NULL);
        return Sim_Output(ID_temp,x_In,x_Out); //x_Out -> Target dynamics output
        }

        //Starting abundances will be same as last abundances
        for(int j=1;j<sim_ana.shape()[1];j++)
            Initial_Abundances[j-1]=sim_ana[sim_ana.shape()[0]-1][j];

        //T_Starting will be same as last time
        T_Starting=i;
    }

    //Something went wrong if we got here
    //return statement with x_Out as an empty dynamics set

    x_Out.Reset_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches);

    Pass_to_Fortran::Thread_Local_Parameters_ptr.reset(NULL);
    return Sim_Output(ID_temp,x_In,x_Out);
}

Meta_Allo_Simulation::Meta_Allo_Simulation(Parameter_Generators x_In, Simulation_Properties Simulation_Properties_Input):
    Par_In(x_In),Sim_In(Simulation_Properties_Input)
{
    int nSets = Par_In.Food_Web->Get_nSets()
            *Par_In.Spatial_Structure->Get_nSets()
            *Par_In.Body_Mass->Get_nSets()
            *Par_In.Metabolism->Get_nSets()
            *Par_In.Functional_Response->Get_nSets()
            *Par_In.Dispersal->Get_nSets()
            *Par_In.Carrying_Capacity->Get_nSets()
            *Par_In.Initial_Abundance->Get_nSets()
            *Par_In.Extinction_Threshold->Get_nSets();

    int nPatches_Full = Par_In.Spatial_Structure->Get_Value().Spatial_Structure.shape()[0];

    if(Simulation_Properties_Input.Test_against_no_Structure)
        nSimulations = nSets*Sim_In.nIterations*(Sim_In.nFood_Webs*(nPatches_Full+Sim_In.nSpatial_Structures));
    else
        nSimulations = nSets*Sim_In.nIterations*Sim_In.nFood_Webs*Sim_In.nSpatial_Structures;

    counter = 0;
    iSet = 0;
}

std::vector<Food_Web_Metadata> Generate_n_Unique_Food_Webs(int n_Food_Webs, const Parameter_Generators &Input_Gen)
{
    Food_Web_Generator* Gen = Input_Gen.Food_Web->Clone();

    int Build_Iteration_Limit = 100;

    std::vector<Food_Web_Metadata> to_Return;

    int Build_Iters = 0;

    while(to_Return.size()<n_Food_Webs)
    {
        Build_Iters++;

        Gen->Generate_New();
        Food_Web_Metadata Test_Web(Gen->Get_Value());

        bool is_Unique = true;

        for(int i2=0;i2<to_Return.size();i2++)
            if(Test_Web.is_Equal(to_Return[i2]))
            {
                is_Unique = false;
                break;
            }

        if(is_Unique)
        {
            to_Return.push_back(Test_Web);
            Build_Iters = 0;
        }

        if(Build_Iters>Build_Iteration_Limit)
            break;
    }

    delete Gen;

    return to_Return;
}

std::vector<Spatial_Structure_Metadata> Generate_n_Unique_Spatial_Structures(int n_Spatial_Structures, const Parameter_Generators &Input_Gen)
{
    Spatial_Structure_Generator* Gen = Input_Gen.Spatial_Structure->Clone();

    int Build_Iteration_Limit = 100;

    std::vector<Spatial_Structure_Metadata> to_Return;

    int Build_Iters = 0;

    while(to_Return.size()<n_Spatial_Structures)
    {
        Build_Iters++;

        Gen->Generate_New();
        Spatial_Structure_Metadata Test_Web(Gen->Get_Value());

        bool is_Unique = true;

        for(int i2=0;i2<to_Return.size();i2++)
            if(Test_Web.is_Equal(to_Return[i2]))
            {
                is_Unique = false;
                break;
            }

        if(is_Unique)
        {
            to_Return.push_back(Test_Web);
            Build_Iters = 0;
        }

        if(Build_Iters>Build_Iteration_Limit)
            break;
    }

    delete Gen;

    return to_Return;
}


void Meta_Allo_Simulation::Run_Careful_Simulation()
{
    omp_lock_t write_mutex;

    omp_init_lock(&write_mutex);
    omp_init_lock(&counter_mutex);

    //std::vector<std::vector<Sim_Output>> Simulation_Results;

    std::vector<Sim_Set_Output> Simulation_Results;

    iSet = 0;

    while(true)
    {
        //Generate unique list of Food_Webs here

        std::vector<Food_Web_Metadata> Food_Webs(Generate_n_Unique_Food_Webs(Sim_In.nFood_Webs,Par_In));
        
        while(true)
        {
            //Generate unique list of Spatial_Structures here

            std::vector<Spatial_Structure_Metadata> Spatial_Structures(Generate_n_Unique_Spatial_Structures(Sim_In.nSpatial_Structures,Par_In));
            
            while(true)
            {while(true)
                {while(true)
                    {while(true)
                        {while(true)
                            {while(true)
                                {
                                        Simulation_Results.push_back(Sim_Set_Output(Par_In));

                                        for(int i_FW = 0;i_FW<Food_Webs.size();i_FW++)
                                        {
                                            Par_In.Food_Web->Set_Value(Food_Webs[i_FW]);

                                            if(Sim_In.Test_against_no_Structure)
                                            {

                                                int nPatches_Full = Par_In.Spatial_Structure->Get_Value().Spatial_Structure.shape()[0];

                                                int nPatches = 1;

                                                Par_In.Spatial_Structure->Set_Value(Spatial_Structure_Metadata(nPatches));

                                                while(true)
                                                {

#pragma omp parallel for schedule(dynamic) num_threads(Sim_In.nThreads)
                                                    for(int i=0;i<Sim_In.nIterations*nPatches_Full;i++)
                                                    {
                                                        Parameter_Generators Thread_Generators(Par_In);

                                                        //Don't generate new Food_Web or Spatial_Structure
                                                        for(int ix=2;ix<Thread_Generators.nGenerators();ix++)
                                                            Thread_Generators[ix]->Generate_New();

                                                        Sim_Output Iter_Out(Simulate(Eqn_Parameters(Thread_Generators),Sim_In));

                                                        Iter_Out.Eqn_In.Food_Web_Unique=i_FW;

                                                        omp_set_lock(&write_mutex);
                                                        Simulation_Results.back().Sim_Out.push_back(Iter_Out);

                                                        int Current_Index = Simulation_Results.back().Sim_Out.size()-1;

                                                        if(Iter_Out.Dyn_Out.Metacommunity_Feasible==1)
                                                            Simulation_Results.back().Feasible_Simulations.push_back(Current_Index);
                                                        else
                                                            Simulation_Results.back().Non_Feasible_Simulations.push_back(Current_Index);

                                                        Simulation_Results.back().Single_Patch_Simulations.push_back(Current_Index);

                                                        omp_unset_lock(&write_mutex);

                                                    }
                                                    if(!Par_In.Extinction_Threshold->Increment_Set())
                                                    {
                                                        Par_In.Extinction_Threshold->Reset();
                                                        break;
                                                    }
                                                }
                                            }

                                            for(int i_SS = 0;i_SS<Spatial_Structures.size();i_SS++)
                                            {
                                                Par_In.Spatial_Structure->Set_Value(Spatial_Structures[i_SS]);

                                                while(true)
                                                {

#pragma omp parallel for schedule(dynamic) num_threads(Sim_In.nThreads)
                                                    for(int i=0;i<Sim_In.nIterations;i++)
                                                    {
                                                        Parameter_Generators Thread_Generators(Par_In);

                                                        //Don't generate new Food_Web or Spatial_Structure
                                                        for(int ix=2;ix<Thread_Generators.nGenerators();ix++)
                                                            Thread_Generators[ix]->Generate_New();

                                                        Sim_Output Iter_Out(Simulate(Eqn_Parameters(Thread_Generators),Sim_In));

                                                        Iter_Out.Eqn_In.Food_Web_Unique=i_FW;
                                                        Iter_Out.Eqn_In.Spatial_Structure_Unique=i_SS;

                                                        omp_set_lock(&write_mutex);
                                                                Simulation_Results.back().Sim_Out.push_back(Iter_Out);

                                                                int Current_Index = Simulation_Results.back().Sim_Out.size()-1;

                                                                if(Iter_Out.Dyn_Out.Metacommunity_Feasible==1)
                                                                    Simulation_Results.back().Feasible_Simulations.push_back(Current_Index);
                                                                else
                                                                    Simulation_Results.back().Non_Feasible_Simulations.push_back(Current_Index);

                                                                Simulation_Results.back().Multi_Patch_Simulations.push_back(Current_Index);

                                                        omp_unset_lock(&write_mutex);

                                                    }
                                                    if(!Par_In.Extinction_Threshold->Increment_Set())
                                                    {
                                                        Par_In.Extinction_Threshold->Reset();
                                                        break;
                                                    }
                                                }
                                            }
                                        }

                                        iSet++;

                                        if(Sim_In.Write_Each_Set)
                                        {
                                            std::string Write_Path_Iterate = Sim_In.Write_Path + std::to_string(static_cast<long long>(iSet)) + "_";

                                            //Write_Population_Cycle_Timeseries(Sim_In.Write_Path, Sim_In.Simulation_Resolution, Simulation_Results);
                                            Write_Population_Output(Write_Path_Iterate,Simulation_Results);
                                            Write_Community_Output(Write_Path_Iterate, false, Simulation_Results);
                                            Write_Community_Parameters(Write_Path_Iterate, false, Simulation_Results,Sim_In.nFood_Webs,Sim_In.nSpatial_Structures);
                                            Write_Set_Summary_Output(Write_Path_Iterate,Simulation_Results);
                                            Write_Food_Web_Summary_Output(Write_Path_Iterate,Simulation_Results,Sim_In.nFood_Webs);
                                            Write_Spatial_Structure_Summary_Output(Write_Path_Iterate,Simulation_Results,Sim_In.nSpatial_Structures);
                                            Write_Patch_Output(Write_Path_Iterate,Simulation_Results);

                                            if(Sim_In.Species_Data_Write>0)
                                            {
                                                Write_Species_Output(Write_Path_Iterate,true,Simulation_Results);

                                                if(Sim_In.Species_Data_Write==2)
                                                {
                                                    Write_Species_Output(Write_Path_Iterate,false,Simulation_Results);
                                                }
                                            }
                                        }

                                    if(!Par_In.Initial_Abundance->Increment_Set())
                                    {
                                        Par_In.Initial_Abundance->Reset();
                                        break;
                                    }
                                }

                                if(!Par_In.Carrying_Capacity->Increment_Set())
                                {
                                    Par_In.Carrying_Capacity->Reset();
                                    break;
                                }
                            }
                            if(!Par_In.Dispersal->Increment_Set())
                            {
                                Par_In.Dispersal->Reset();
                                break;
                            }
                        }
                        if(!Par_In.Functional_Response->Increment_Set())
                        {
                            Par_In.Functional_Response->Reset();
                            break;
                        }
                    }
                    if(!Par_In.Metabolism->Increment_Set())
                    {
                        Par_In.Metabolism->Reset();
                        break;
                    }
                }
                if(!Par_In.Body_Mass->Increment_Set())
                {
                    Par_In.Body_Mass->Reset();
                    break;
                }
            }
            if(!Par_In.Spatial_Structure->Increment_Set())
            {
                Par_In.Spatial_Structure->Reset();
                break;
            }
        }
        if(!Par_In.Food_Web->Increment_Set())
        {
            Par_In.Food_Web->Reset();
            break;
        }
    }

    nSimulations = counter;

    omp_destroy_lock(&write_mutex);
    omp_destroy_lock(&counter_mutex);

    if(!Sim_In.Write_Each_Set)
    {
        //Write_Population_Cycle_Timeseries(Sim_In.Write_Path, Sim_In.Simulation_Resolution, Simulation_Results);
        Write_Population_Output(Sim_In.Write_Path,Simulation_Results);
        Write_Community_Output(Sim_In.Write_Path, false, Simulation_Results);
        Write_Community_Parameters(Sim_In.Write_Path, false, Simulation_Results,Sim_In.nFood_Webs,Sim_In.nSpatial_Structures);
        Write_Set_Summary_Output(Sim_In.Write_Path,Simulation_Results);
        Write_Food_Web_Summary_Output(Sim_In.Write_Path,Simulation_Results,Sim_In.nFood_Webs);
        Write_Spatial_Structure_Summary_Output(Sim_In.Write_Path,Simulation_Results,Sim_In.nSpatial_Structures);
        Write_Patch_Output(Sim_In.Write_Path,Simulation_Results);

        if(Sim_In.Species_Data_Write>0)
        {
            Write_Species_Output(Sim_In.Write_Path,true,Simulation_Results);

            if(Sim_In.Species_Data_Write==2)
            {
                Write_Species_Output(Sim_In.Write_Path,false,Simulation_Results);
            }
        }
    }
}

Sim_Set_Output::Sim_Set_Output(const Parameter_Generators &Set_Gens)
{
    for(int i=0;i<Set_Gens.nGenerators();i++)
    {
        Set_Gen_Parameter_Labels.push_back(Set_Gens[i]->Get_Generator_Method_Labels());
        Set_Gen_Parameter_Values.push_back(Set_Gens[i]->Get_Generator_Method_Data());
    }
}

Eqn_Parameters::Eqn_Parameters()
{
    Community_Data_Labels.push_back("Unique_Food_Web");
    Community_Data_Labels.push_back("H");
    Community_Data_Labels.push_back("AMI");
    Community_Data_Labels.push_back("Hc");
    Community_Data_Labels.push_back("Allo_H");
    Community_Data_Labels.push_back("Allo_AMI");
    Community_Data_Labels.push_back("Allo_Hc");
    Community_Data_Labels.push_back("Min_Consumer_Resource_Size_Ratio");
    Community_Data_Labels.push_back("Max_Consumer_Resource_Size_Ratio");
    Community_Data_Labels.push_back("Avg_Consumer_Resource_Size_Ratio");
    Community_Data_Labels.push_back("Unique_Spatial_Structure");
    Community_Data_Labels.push_back("Average_Degree");
    Community_Data_Labels.push_back("Degree_Skewness");
    Community_Data_Labels.push_back("Average_Path_Length");
    Community_Data_Labels.push_back("Clustering_Coefficient");
    Community_Data_Labels.push_back("Eigenratio");
    Community_Data_Labels.push_back("Extinction_Threshold");

    Species_Data_Labels.push_back("Min_Trophic_Position");
    Species_Data_Labels.push_back("Avg_Trophic_Position");
    Species_Data_Labels.push_back("Max_Trophic_Position");
    Species_Data_Labels.push_back("Min_Consumer_Resource_Size_Ratio");
    Species_Data_Labels.push_back("Max_Consumer_Resource_Size_Ratio");
    Species_Data_Labels.push_back("Avg_Consumer_Resource_Size_Ratio");
    Species_Data_Labels.push_back("Generality");
    Species_Data_Labels.push_back("Positional_Index");
    Species_Data_Labels.push_back("Mass");
    Species_Data_Labels.push_back("ar");
    Species_Data_Labels.push_back("ax");
    Species_Data_Labels.push_back("y");
    Species_Data_Labels.push_back("e");
    Species_Data_Labels.push_back("B0");
    Species_Data_Labels.push_back("h");
    Species_Data_Labels.push_back("Dispersal");

    Patches_Data_Labels.push_back("Patch_Size");

    n_Community_Parameters = 17;
    n_Species_Parameters = 16;
    n_Patches_Parameters = 1;
}

Eqn_Parameters::Eqn_Parameters(Parameter_Generators Input_Generators):
    Food_Web_Data(Input_Generators.Food_Web->Get_Value()),
    Spatial_Structure_Data(Input_Generators.Spatial_Structure->Get_Value()),
    k(Input_Generators.Carrying_Capacity->Get_K()),
    Patch_Size(Input_Generators.Carrying_Capacity->Get_Patch_Size()),
    Functional_Response_Shape(Input_Generators.Functional_Response->Get_h()),
    Half_Saturation(Input_Generators.Functional_Response->Get_B0()),
    y(Input_Generators.Metabolism->Get_y()),
    ar(Input_Generators.Metabolism->Get_ar()),
    ax(Input_Generators.Metabolism->Get_ax()),
    e(Input_Generators.Metabolism->Get_e()),
    M(Input_Generators.Body_Mass->Get_Value()),
    m(Input_Generators.Dispersal->Get_Value()),
    Initial_Abundance(Input_Generators.Initial_Abundance->Get_Value()),
    Extinction_Threshold(Input_Generators.Extinction_Threshold->Get_Value()),
    Food_Web_Unique(-1),
    Spatial_Structure_Unique(-1)
{
    nSpecies = Food_Web_Data.Food_Web.shape()[0];
    nPatches = Spatial_Structure_Data.Spatial_Structure.shape()[0];

    boost_matrix Complete_Food_Web(Food_Web_Data.Food_Web);

    for(int i1 =0;i1<Complete_Food_Web.shape()[0];i1++){
        for(int i2=0;i2<Complete_Food_Web.shape()[1];i2++){

        Complete_Food_Web[i2][i1] = ax[i1]*pow(M[i1],-.25)*y[i1]/e[i1]*Complete_Food_Web[i2][i1];

        }
    }

    Food_Web_Metadata Complete_Food_Web_Metadata(Complete_Food_Web);

    Allometric_Scaled_H = Complete_Food_Web_Metadata.H;
    Allometric_Scaled_AMI = Complete_Food_Web_Metadata.AMI;
    Allometric_Scaled_Hc = Complete_Food_Web_Metadata.Hc;

    std::vector<double> All_Interactions_Consumer_Resource_Size_Ratios;

    for(int i_Sp=0;i_Sp<Complete_Food_Web.shape()[0];i_Sp++)
    {
        std::vector<double> Consumer_Resource_Size_Ratios;

        for(int i_Prey=0;i_Prey<Complete_Food_Web.shape()[1];i_Prey++)
        {
            if(Complete_Food_Web[i_Prey][i_Sp]>0)
            {
                Consumer_Resource_Size_Ratios.push_back(M[i_Sp]/M[i_Prey]);
                All_Interactions_Consumer_Resource_Size_Ratios.push_back(M[i_Sp]/M[i_Prey]);
            }
        }

        if(Consumer_Resource_Size_Ratios.size()>0)
        {

            Species_Min_Consumer_Resource_Size_Ratio.push_back(*std::min_element(Consumer_Resource_Size_Ratios.begin(),Consumer_Resource_Size_Ratios.end()));
            Species_Max_Consumer_Resource_Size_Ratio.push_back(*std::max_element(Consumer_Resource_Size_Ratios.begin(),Consumer_Resource_Size_Ratios.end()));
            Species_Avg_Consumer_Resource_Size_Ratio.push_back(std::accumulate(Consumer_Resource_Size_Ratios.begin(),Consumer_Resource_Size_Ratios.end(),0.0)/Consumer_Resource_Size_Ratios.size());
        }
        else
        {
            Species_Min_Consumer_Resource_Size_Ratio.push_back(0);
            Species_Max_Consumer_Resource_Size_Ratio.push_back(0);
            Species_Avg_Consumer_Resource_Size_Ratio.push_back(0);
        }
    }

    if(All_Interactions_Consumer_Resource_Size_Ratios.size()>0)
    {

        Community_Min_Consumer_Resource_Size_Ratio = *std::min_element(All_Interactions_Consumer_Resource_Size_Ratios.begin(),All_Interactions_Consumer_Resource_Size_Ratios.end());
        Community_Max_Consumer_Resource_Size_Ratio = *std::max_element(All_Interactions_Consumer_Resource_Size_Ratios.begin(),All_Interactions_Consumer_Resource_Size_Ratios.end());
        Community_Avg_Consumer_Resource_Size_Ratio = std::accumulate(All_Interactions_Consumer_Resource_Size_Ratios.begin(),All_Interactions_Consumer_Resource_Size_Ratios.end(),0.0)/All_Interactions_Consumer_Resource_Size_Ratios.size();

    }
    else
    {
        Community_Min_Consumer_Resource_Size_Ratio = 0;
        Community_Max_Consumer_Resource_Size_Ratio = 0;
        Community_Avg_Consumer_Resource_Size_Ratio = 0;
    }

    Community_Data_Labels.push_back("Unique_Food_Web");
    Community_Data_Labels.push_back("H");
    Community_Data_Labels.push_back("AMI");
    Community_Data_Labels.push_back("Hc");
    Community_Data_Labels.push_back("Allo_H");
    Community_Data_Labels.push_back("Allo_AMI");
    Community_Data_Labels.push_back("Allo_Hc");
    Community_Data_Labels.push_back("Min_Consumer_Resource_Size_Ratio");
    Community_Data_Labels.push_back("Max_Consumer_Resource_Size_Ratio");
    Community_Data_Labels.push_back("Avg_Consumer_Resource_Size_Ratio");
    Community_Data_Labels.push_back("Unique_Spatial_Structure");
    Community_Data_Labels.push_back("Average_Degree");
    Community_Data_Labels.push_back("Degree_Skewness");
    Community_Data_Labels.push_back("Average_Path_Length");
    Community_Data_Labels.push_back("Clustering_Coefficient");
    Community_Data_Labels.push_back("Eigenratio");
    Community_Data_Labels.push_back("Extinction_Threshold");

    Species_Data_Labels.push_back("Min_Trophic_Position");
    Species_Data_Labels.push_back("Avg_Trophic_Position");
    Species_Data_Labels.push_back("Max_Trophic_Position");
    Species_Data_Labels.push_back("Min_Consumer_Resource_Size_Ratio");
    Species_Data_Labels.push_back("Max_Consumer_Resource_Size_Ratio");
    Species_Data_Labels.push_back("Avg_Consumer_Resource_Size_Ratio");
    Species_Data_Labels.push_back("Generality");
    Species_Data_Labels.push_back("Positional_Index");
    Species_Data_Labels.push_back("Mass");
    Species_Data_Labels.push_back("ar");
    Species_Data_Labels.push_back("ax");
    Species_Data_Labels.push_back("y");
    Species_Data_Labels.push_back("e");
    Species_Data_Labels.push_back("B0");
    Species_Data_Labels.push_back("h");
    Species_Data_Labels.push_back("Dispersal");

    Patches_Data_Labels.push_back("Patch_Size");

    n_Community_Parameters = 17;
    n_Species_Parameters = 16;
    n_Patches_Parameters = 1;
}

int Eqn_Parameters::n_Community_Data() const
{
    return n_Community_Parameters;
}

int Eqn_Parameters::n_Species_Data() const
{
    return n_Species_Parameters;
}

int Eqn_Parameters::n_Patches_Data() const
{
    return n_Patches_Parameters;
}

double Eqn_Parameters::Get_Community_Data(int Index) const
{
    switch(Index)
    {
    case(0):
    {
        return Food_Web_Unique;
    }
        break;
    case(1):
    {
        return Food_Web_Data.H;
    }
        break;
    case(2):
    {
        return Food_Web_Data.AMI;
    }
        break;
    case(3):
    {
        return Food_Web_Data.Hc;
    }
        break;
    case(4):
    {
        return Allometric_Scaled_H;
    }
        break;
    case(5):
    {
        return Allometric_Scaled_AMI;
    }
        break;
    case(6):
    {
        return Allometric_Scaled_Hc;
    }
        break;
    case(7):
    {
        return Community_Min_Consumer_Resource_Size_Ratio;
    }
        break;
    case(8):
    {
        return Community_Max_Consumer_Resource_Size_Ratio;
    }
        break;
    case(9):
    {
        return Community_Avg_Consumer_Resource_Size_Ratio;
    }
        break;
    case(10):
    {
        return Spatial_Structure_Unique;
    }
        break;
    case(11):
    {
        return Spatial_Structure_Data.Average_Degree;
    }
        break;
    case(12):
    {
        return Spatial_Structure_Data.Degree_Skewness;
    }
        break;
    case(13):
    {
        return Spatial_Structure_Data.Average_Path_Length;
    }
        break;
    case(14):
    {
        return Spatial_Structure_Data.Clustering_Coefficient;
    }
        break;
    case(15):
    {
        return Spatial_Structure_Data.Eigenratio;
    }
        break;
    case(16):
    {
        return Extinction_Threshold;
    }
        break;

    }

    return 0;
}

std::vector<double> Eqn_Parameters::Get_Species_Data(int Index) const
{    
    switch(Index)
    {
    case(0):
    {
        return Food_Web_Data.Min_Trophic_Position;
    }
    case(1):
    {
        return Food_Web_Data.Avg_Trophic_Position;
    }
    case(2):
    {
        return Food_Web_Data.Max_Trophic_Position;
    }
        break;
    case(3):
    {
        return Species_Min_Consumer_Resource_Size_Ratio;
    }
    case(4):
    {
        return Species_Max_Consumer_Resource_Size_Ratio;
    }
    case(5):
    {
        return Species_Avg_Consumer_Resource_Size_Ratio;
    }
        break;
    case(6):
    {
        return Food_Web_Data.Generality;
    }
        break;
    case(7):
    {
        return Food_Web_Data.Positional_Index;
    }
        break;
    case(8):
    {
        return M;
    }
        break;
    case(9):
    {
        return ar;
    }
        break;
    case(10):
    {
        return ax;
    }
        break;
    case(11):
    {
        return y;
    }
        break;
    case(12):
    {
        return e;
    }
        break;
    case(13):
    {
        return Half_Saturation;
    }
        break;
    case(14):
    {
        return Functional_Response_Shape;
    }
        break;
    case(15):
    {
        return m;
    }
        break;
    }

    return std::vector<double>(1,0);
}

std::vector<double> Eqn_Parameters::Get_Patches_Data(int Index) const
{
    switch(Index)
    {
    case(0):
    {
        return Patch_Size;
    }
        break;
    }

    return std::vector<double>(1,0);
}

