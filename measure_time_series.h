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
#ifndef MEASURE_TIME_SERIES
#define MEASURE_TIME_SERIES

#include "find_cycles.h"
//#include "allometric_metacommunity.h"

struct Simulation_Properties {

    int nThreads;
    int nIterations;
    int nFood_Webs;
    int nSpatial_Structures;
    double Simulation_Resolution;
    std::vector<double> Extinction_Threshold;

    std::string Write_Path;
    int Time_Series_Write;
    int Limit_Cycle_Write;
    int Species_Data_Write;
    bool Write_Each_Set;
    bool Test_against_no_Structure;
    bool Write_Measure_Diagnostics;

    double Minimum_Simulation_Duration;
    double Maximum_Simulation_Duration;
    double Equilibrium_Test_Interval_Duration;
    double Time_between_Equilibrium_Tests;

    double Phase_Lock_Threshold;

    Cycle_Detection_Parameters Cycle_Params;
};

struct Dynamics {

    //Metacommunity Metrics
    int Metacommunity_Feasible;
    int Metacommunity_Solver_Error;
    double Metacommunity_Stationary;
    double Metacommunity_Phase_Locked;
    double Metacommunity_Fixed_Point;

    double Metacommunity_Prop_Measured;
//    double Metacommunity_Avg_Extrema_Variability;
//    double Metacommunity_Avg_Extrema_Center_Fit;

    double Metacommunity_Transient_Length;
    double Metacommunity_Asynchronous;
    int Metacommunity_Eqb_Synch_Clusters;
    int Metacommunity_Eqb_Synch_Amp_Clusters;

    double Metacommunity_Prop_Eqb_Cycle;

    double Metacommunity_Eqb_Cycle_N_Extrema; //

    double Metacommunity_Eqb_Min_Max_Minima;
    double Metacommunity_Avg_Max_Amplitude;
    double Metacommunity_Avg_Eqb;

    double Metacommunity_Avg_CV_Max_Amplitude;
    double Metacommunity_Avg_CV_Eqb;

    double Metacommunity_Avg_Period;
    double Metacommunity_Avg_CV_Period;

    double Metacommunity_Prop_Common_Phase_Difference;

    double Metacommunity_Phase_Difference_Cycle_N_Diffs; //

    double Metacommunity_Avg_Common_Phase_Difference;
    double Metacommunity_Avg_SD_Common_Phase_Difference;

    double Metacommunity_Avg_Phase_Difference;

    double Metacommunity_Avg_Range_Common_Phase_Difference;
    double Metacommunity_Avg_SD_Avg_Phase_Difference;

    double Metacommunity_Quotient_Connec_In;
    double Metacommunity_Quotient_Connec_Out;
    double Metacommunity_Quotient_Skew_Connec_In;
    double Metacommunity_Quotient_Skew_Connec_Out;
    double Metacommunity_Quotient_Clustering;
    double Metacommunity_Quotient_Mean_Path_Length;
    double Metacommunity_Quotient_Eigenratio;

    //Species Metrics
    std::vector<int> Species_Extinct;
    std::vector<double> Species_Stationary;
    std::vector<double> Species_Phase_Locked;
    std::vector<double> Species_Fixed_Point;

    std::vector<double> Species_Prop_Measured;
    std::vector<double> Species_Asynchronous;
    std::vector<int> Species_Eqb_Synch_Clusters;
    std::vector<int> Species_Eqb_Synch_Amp_Clusters;

    std::vector<double> Species_Prop_Eqb_Cycle;

    std::vector<double> Species_Eqb_Cycle_N_Extrema; 

    std::vector<double> Species_Max_Eqb_Minima;
    std::vector<double> Species_Avg_Max_Amplitude;
    std::vector<double> Species_Avg_Eqb;

    std::vector<double> Species_CV_Max_Amplitude;
    std::vector<double> Species_CV_Eqb;

    std::vector<double> Species_Avg_Period;
    std::vector<double> Species_CV_Period;

    std::vector<double> Species_Prop_Common_Phase_Difference;

    std::vector<double> Species_Phase_Difference_Cycle_N_Diffs; 

    std::vector<double> Species_Avg_Common_Phase_Difference;
    std::vector<double> Species_Avg_SD_Common_Phase_Difference;

    std::vector<double> Species_Avg_Phase_Difference;

    std::vector<double> Species_Range_Common_Phase_Difference;
    std::vector<double> Species_SD_Avg_Phase_Difference;
    
    //New
    std::vector<double> Species_Total_Average;
    std::vector<double> Species_Total_CV;
    std::vector<double> Species_Total_Minimum;
    std::vector<double> Species_Total_Maximum;
    //End New

    //Community Metrics

    std::vector<int> Community_Eqb_Synch_Cluster_ID;
    std::vector<int> Community_Eqb_Synch_Amp_Cluster_ID;

    std::vector<double> Community_Prop_Measured;

    std::vector<double> Community_Prop_Eqb_Cycle;

    std::vector<double> Community_Eqb_N_Cycle_Extrema; //

    std::vector<double> Community_Eqb_Minima;
    std::vector<double> Community_Avg_Max_Amplitude;
    std::vector<double> Community_Avg_Eqb;

    std::vector<double> Community_Avg_Period;
    std::vector<double> Community_Prop_Common_Phase_Difference;

    std::vector<double> Community_Phase_Difference_Cycle_N_Diffs; //

    std::vector<double> Community_Avg_Common_Phase_Difference;
    std::vector<double> Community_Avg_SD_Common_Phase_Difference;

    std::vector<double> Community_Avg_Phase_Difference;

    std::vector<double> Community_Range_Common_Phase_Difference;
    std::vector<double> Community_SD_Avg_Phase_Difference;

    //Population Metrics

    std::vector<std::vector<Cycle_Details>> Population_Metrics; //[Species][Patches]
    std::vector<std::vector<int>> Population_Eqb_Synch_Cluster_ID;
    std::vector<std::vector<int>> Population_Eqb_Synch_Amp_Cluster_ID;
    
    //All New
    std::vector<std::vector<double>> Population_Time_Series_Average;
    std::vector<std::vector<double>> Population_Time_Series_CV;
    std::vector<std::vector<double>> Population_Time_Series_Skewness;
    std::vector<std::vector<double>> Population_Time_Series_Minimum;
    std::vector<std::vector<double>> Population_Time_Series_Maximum;
    //New End

//Data management
//Fix to new standard
    int n_Community_Metrics;
    int n_Species_Metrics;
    int n_Patch_Metrics;
    int n_Population_Metrics;

    std::vector<std::string> Community_Dynamics_Labels;
    std::vector<std::string> Species_Dynamics_Labels;
    std::vector<std::string> Patch_Dynamics_Labels;
    std::vector<std::string> Population_Dynamics_Labels;

    void Reset_Dynamics(int nSpecies, int nPatches);

    Dynamics();

    double Get_Community_Measures(int Index) const;
    std::vector<double> Get_Species_Measures(int Index) const;
    std::vector<double> Get_Patch_Measures(int Index) const;
    std::vector<std::vector<double>> Get_Population_Measures(int Index) const;
    std::vector<std::pair<double,double>> Get_Population_Time_Series(int Species, int Patch) const;

    int n_Community_Measures() const;
    int n_Species_Measures() const;
    int n_Patch_Measures() const;
    int n_Population_Measures() const;
};

int count_distinct_abs(std::vector<int> v);

std::vector<int> Which_Species_Extinct(int nSpecies, boost_matrix Simulation_Data, double Extinction_Threshold);

//New Seperate Dynamics gen for Nonfeas metacommunities
//Dynamics Measure_Metacommunity_Dynamics(const std::vector<Find_Cycles> &Cycles);

void Measure_Time_Series_Properties(const boost_matrix_2d_view &Species_Time_Series, int i_Species, Dynamics &Output);

bool Measure_Metacommunity_Dynamics(int nSpecies,
                                    double Transient_Test_Time,
                                    std::pair<int,int> Test_Interval,
                                    boost_matrix &Simulation_Data, bool Solver_Feas_Flag,
                                    double Extinction_Threshold,
                                    std::string Write_Path,
                                    const boost_matrix &Spatial_Structure,
                                    const Simulation_Properties &Sim_Props,
                                    bool Measure_if_Nonstationary,
                                    Dynamics &Output);

#endif // MEASURE_TIME_SERIES

