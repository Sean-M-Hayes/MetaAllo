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
#include "measure_time_series.h"
#include "spatial_structure_methods.h"
#include <fstream>

int count_distinct_abs(std::vector<int> v)
{
    sort(v.begin(), v.end()); // Average case O(n log n), worst case O(n^2) (usually implemented as quicksort.
    // To guarantee worst case O(n log n) replace with make_heap, then sort_heap.

    // Unique will take a sorted range, and move things around to get duplicated
    // items to the back and returns an iterator to the end of the unique section of the range
    auto unique_end = unique(v.begin(), v.end()); // Again n comparisons
    return distance(v.begin(),unique_end); // Constant time for random access iterators (like vector's)
}

std::vector< double > vSequence (double initial, double final, double by){

    std::vector< double > gen_Sequence;

    if(initial<=final)
    {
        for(double i = initial;i<=final;i=i+by)
            gen_Sequence.push_back(i);
    }
    else
    {
        for(double i = initial; i>final;i=i-by)
            gen_Sequence.push_back(i);
    }

    return gen_Sequence;
}

std::vector< int > vSequence (int initial, int final, int by){

    std::vector< int > gen_Sequence;

    if(initial<=final)
    {
        for(int i = initial;i<=final;i=i+by)
            gen_Sequence.push_back(i);
    }
    else
    {
        for(int i = initial; i>final;i=i-by)
            gen_Sequence.push_back(i);
    }

    return gen_Sequence;
}

std::vector<int> Which_Species_Extinct(int nSpecies, boost_matrix Simulation_Data, double Extinction_Threshold)
{
    std::vector<int> Extinct_Species(nSpecies,1);

    for(unsigned i = 0;i<Simulation_Data.shape()[1];i++)
    {
        if(Simulation_Data[Simulation_Data.shape()[0]-1][i]>Extinction_Threshold)
            Extinct_Species[i%nSpecies]=0;
    }

    return Extinct_Species;
}

void Measure_Nonfeasible_Metacommunity_Dynamics(int nSpecies, const boost_matrix &Simulation_Data, double Extinction_Threshold, Dynamics &Output)
{
    std::vector<int> Extinct_Species(nSpecies,1);

    for(unsigned i = 0;(i+1)<Simulation_Data.shape()[1];i++)
    {
        if(Simulation_Data[Simulation_Data.shape()[0]-1][i+1]>Extinction_Threshold)
            Extinct_Species[i%nSpecies]=0;
    }

    Output.Species_Extinct = Extinct_Species;
}

void Dynamics::Reset_Dynamics(int nSpecies, int nPatches)
{
    Metacommunity_Solver_Error = 0;

    Metacommunity_Stationary = 0;
    Metacommunity_Phase_Locked = 0;

    Metacommunity_Prop_Measured = 0;

    Metacommunity_Transient_Length = 0;

    Metacommunity_Feasible = 0;
    Metacommunity_Fixed_Point = 0;
    Metacommunity_Asynchronous = 0;
    Metacommunity_Eqb_Synch_Clusters = 0;
    Metacommunity_Eqb_Synch_Amp_Clusters = 0;

    Metacommunity_Prop_Eqb_Cycle = 0;

    Metacommunity_Eqb_Cycle_N_Extrema = 0;

    Metacommunity_Eqb_Min_Max_Minima = 0;
    Metacommunity_Avg_Max_Amplitude = 0;
    Metacommunity_Avg_Eqb = 0;

    Metacommunity_Avg_CV_Max_Amplitude = 0;
    Metacommunity_Avg_CV_Eqb = 0;

    Metacommunity_Avg_Period = 0;
    Metacommunity_Avg_CV_Period = 0;

    Metacommunity_Prop_Common_Phase_Difference = 0;

    Metacommunity_Phase_Difference_Cycle_N_Diffs = 0;

    Metacommunity_Avg_Common_Phase_Difference = 0;
    Metacommunity_Avg_SD_Common_Phase_Difference = 0;

    Metacommunity_Avg_Phase_Difference = 0;

    Metacommunity_Avg_Range_Common_Phase_Difference = 0;
    Metacommunity_Avg_SD_Avg_Phase_Difference = 0;

    //Metacommunity Quotient Network Measures
    Metacommunity_Quotient_Connec_In = 0;
    Metacommunity_Quotient_Connec_Out = 0;
    Metacommunity_Quotient_Skew_Connec_In = 0;
    Metacommunity_Quotient_Skew_Connec_Out = 0;
    Metacommunity_Quotient_Clustering = 0;
    Metacommunity_Quotient_Mean_Path_Length = 0;
    Metacommunity_Quotient_Eigenratio = 0;

    //Species Metrics

    Species_Prop_Measured.resize(nSpecies);

    Species_Extinct.resize(nSpecies);
    Species_Stationary.resize(nSpecies);
    Species_Phase_Locked.resize(nSpecies);
    Species_Fixed_Point.resize(nSpecies);
    Species_Asynchronous.resize(nSpecies);
    Species_Eqb_Synch_Clusters.resize(nSpecies);
    Species_Eqb_Synch_Amp_Clusters.resize(nSpecies);

    Species_Prop_Eqb_Cycle.resize(nSpecies);

    Species_Eqb_Cycle_N_Extrema.resize(nSpecies);

    Species_Max_Eqb_Minima.resize(nSpecies);
    Species_Avg_Max_Amplitude.resize(nSpecies);
    Species_Avg_Eqb.resize(nSpecies);

    Species_CV_Max_Amplitude.resize(nSpecies);
    Species_CV_Eqb.resize(nSpecies);

    Species_Avg_Period.resize(nSpecies);
    Species_CV_Period.resize(nSpecies);

    Species_Prop_Common_Phase_Difference.resize(nSpecies);

    Species_Phase_Difference_Cycle_N_Diffs.resize(nSpecies);

    Species_Avg_Common_Phase_Difference.resize(nSpecies);
    Species_Avg_SD_Common_Phase_Difference.resize(nSpecies);

    Species_Avg_Phase_Difference.resize(nSpecies);

    Species_Range_Common_Phase_Difference.resize(nSpecies);
    Species_SD_Avg_Phase_Difference.resize(nSpecies);
    
    Species_Total_Average.resize(nSpecies);
    Species_Total_CV.resize(nSpecies);
    Species_Total_Minimum.resize(nSpecies);
    Species_Total_Maximum.resize(nSpecies);

    std::fill(Species_Prop_Measured.begin(),Species_Prop_Measured.end(),0);

    std::fill(Species_Extinct.begin(),Species_Extinct.end(),0);
    std::fill(Species_Stationary.begin(),Species_Stationary.end(),0);
    std::fill(Species_Phase_Locked.begin(),Species_Phase_Locked.end(),0);
    std::fill(Species_Fixed_Point.begin(),Species_Fixed_Point.end(),0);
    std::fill(Species_Asynchronous.begin(),Species_Asynchronous.end(),0);
    std::fill(Species_Eqb_Synch_Clusters.begin(),Species_Eqb_Synch_Clusters.end(),0);
    std::fill(Species_Eqb_Synch_Amp_Clusters.begin(),Species_Eqb_Synch_Amp_Clusters.end(),0);

    std::fill(Species_Prop_Eqb_Cycle.begin(),Species_Prop_Eqb_Cycle.end(),0);

    std::fill(Species_Eqb_Cycle_N_Extrema.begin(),Species_Eqb_Cycle_N_Extrema.end(),0);

    std::fill(Species_Max_Eqb_Minima.begin(),Species_Max_Eqb_Minima.end(),0);
    std::fill(Species_Avg_Max_Amplitude.begin(),Species_Avg_Max_Amplitude.end(),0);
    std::fill(Species_Avg_Eqb.begin(),Species_Avg_Eqb.end(),0);

    std::fill(Species_CV_Max_Amplitude.begin(),Species_CV_Max_Amplitude.end(),0);
    std::fill(Species_CV_Eqb.begin(),Species_CV_Eqb.end(),0);

    std::fill(Species_Avg_Period.begin(),Species_Avg_Period.end(),0);
    std::fill(Species_CV_Period.begin(),Species_CV_Period.end(),0);

    std::fill(Species_Prop_Common_Phase_Difference.begin(),Species_Prop_Common_Phase_Difference.end(),0);

    std::fill(Species_Phase_Difference_Cycle_N_Diffs.begin(),Species_Phase_Difference_Cycle_N_Diffs.end(),0);

    std::fill(Species_Avg_Common_Phase_Difference.begin(),Species_Avg_Common_Phase_Difference.end(),0);
    std::fill(Species_Avg_SD_Common_Phase_Difference.begin(),Species_Avg_SD_Common_Phase_Difference.end(),0);

    std::fill(Species_Avg_Phase_Difference.begin(),Species_Avg_Phase_Difference.end(),0);

    std::fill(Species_Range_Common_Phase_Difference.begin(),Species_Range_Common_Phase_Difference.end(),0);
    std::fill(Species_SD_Avg_Phase_Difference.begin(),Species_SD_Avg_Phase_Difference.end(),0);
    
    std::fill(Species_Total_Average.begin(),Species_Total_Average.end(),0);
    std::fill(Species_Total_CV.begin(),Species_Total_CV.end(),0);
    std::fill(Species_Total_Minimum.begin(),Species_Total_Minimum.end(),0);
    std::fill(Species_Total_Maximum.begin(),Species_Total_Maximum.end(),0);

    //Community Metrics

    Community_Prop_Measured.resize(nPatches);

    Community_Eqb_Synch_Cluster_ID.resize(nPatches);
    Community_Eqb_Synch_Amp_Cluster_ID.resize(nPatches);

    Community_Prop_Eqb_Cycle.resize(nPatches);

    Community_Eqb_N_Cycle_Extrema.resize(nPatches);

    Community_Eqb_Minima.resize(nPatches);
    Community_Avg_Max_Amplitude.resize(nPatches);
    Community_Avg_Eqb.resize(nPatches);

    Community_Avg_Period.resize(nPatches);
    Community_Prop_Common_Phase_Difference.resize(nPatches);

    Community_Phase_Difference_Cycle_N_Diffs.resize(nPatches);

    Community_Avg_Common_Phase_Difference.resize(nPatches);
    Community_Avg_SD_Common_Phase_Difference.resize(nPatches);

    Community_Avg_Phase_Difference.resize(nPatches);

    Community_Range_Common_Phase_Difference.resize(nPatches);
    Community_SD_Avg_Phase_Difference.resize(nPatches);

    std::fill(Community_Prop_Measured.begin(),Community_Prop_Measured.end(),0);

    std::fill(Community_Eqb_Synch_Cluster_ID.begin(),Community_Eqb_Synch_Cluster_ID.end(),0);
    std::fill(Community_Eqb_Synch_Amp_Cluster_ID.begin(),Community_Eqb_Synch_Amp_Cluster_ID.end(),0);

    std::fill(Community_Prop_Eqb_Cycle.begin(),Community_Prop_Eqb_Cycle.end(),0);

    std::fill(Community_Eqb_N_Cycle_Extrema.begin(),Community_Eqb_N_Cycle_Extrema.end(),0);

    std::fill(Community_Eqb_Minima.begin(),Community_Eqb_Minima.end(),0);
    std::fill(Community_Avg_Max_Amplitude.begin(),Community_Avg_Max_Amplitude.end(),0);
    std::fill(Community_Avg_Eqb.begin(),Community_Avg_Eqb.end(),0);

    std::fill(Community_Avg_Period.begin(),Community_Avg_Period.end(),0);
    std::fill(Community_Prop_Common_Phase_Difference.begin(),Community_Prop_Common_Phase_Difference.end(),0);

    std::fill(Community_Phase_Difference_Cycle_N_Diffs.begin(),Community_Phase_Difference_Cycle_N_Diffs.end(),0);

    std::fill(Community_Avg_Common_Phase_Difference.begin(),Community_Avg_Common_Phase_Difference.end(),0);
    std::fill(Community_Avg_SD_Common_Phase_Difference.begin(),Community_Avg_SD_Common_Phase_Difference.end(),0);

    std::fill(Community_Avg_Phase_Difference.begin(),Community_Avg_Phase_Difference.end(),0);

    std::fill(Community_Range_Common_Phase_Difference.begin(),Community_Range_Common_Phase_Difference.end(),0);
    std::fill(Community_SD_Avg_Phase_Difference.begin(),Community_SD_Avg_Phase_Difference.end(),0);

    //Population Metrics

    Population_Metrics.clear();
    Population_Eqb_Synch_Cluster_ID.clear();
    Population_Eqb_Synch_Amp_Cluster_ID.clear();
    
    for(int i=0;i<nSpecies;i++)
    {
        Population_Metrics.push_back(std::vector<Cycle_Details>(nPatches));
        Population_Eqb_Synch_Cluster_ID.push_back(std::vector<int>(nPatches,0));
        Population_Eqb_Synch_Amp_Cluster_ID.push_back(std::vector<int>(nPatches,0));
    }

    Population_Time_Series_Average.resize(nSpecies);
    Population_Time_Series_CV.resize(nSpecies);
    Population_Time_Series_Skewness.resize(nSpecies);
    Population_Time_Series_Minimum.resize(nSpecies);
    Population_Time_Series_Maximum.resize(nSpecies);

    std::fill(Population_Time_Series_Average.begin(),Population_Time_Series_Average.end(),std::vector<double>(nPatches,0));
    std::fill(Population_Time_Series_CV.begin(),Population_Time_Series_CV.end(),std::vector<double>(nPatches,0));
    std::fill(Population_Time_Series_Skewness.begin(),Population_Time_Series_Skewness.end(),std::vector<double>(nPatches,0));
    std::fill(Population_Time_Series_Minimum.begin(),Population_Time_Series_Minimum.end(),std::vector<double>(nPatches,0));
    std::fill(Population_Time_Series_Maximum.begin(),Population_Time_Series_Maximum.end(),std::vector<double>(nPatches,0));

}

void Write_Measure_Diagnostics(std::string Write_Path, const Find_Cycles &Output)
{
    std::ofstream Diagnostic_Output(Write_Path+"_Diagnostic.csv");

    int max_n_Rows = 0;

    std::vector<std::vector<std::vector<int>>> Element_Sequence;

    for(int i=0;i<Output.Patterns_to_Compare.size();i++)
    {
        Diagnostic_Output<<"Matched_"<<i<<",Element_"<<i<<",Dist_"<<i<<",First_Max_"<<i<<",Min_"<<i<<",Second_Max_"<<i<<",";

        if(Output.Sequences_to_Compare[i].Get_Sequence().size()>max_n_Rows)
            max_n_Rows = Output.Sequences_to_Compare[i].Get_Sequence().size();

        Element_Sequence.push_back(Output.Sequences_to_Compare[i].Get_Sequence());
    }

    Diagnostic_Output<<std::endl;

    std::vector<std::vector<int>> Matched;

    //std::vector<std::vector<std::vector<int>>> Element_Sequence = Output.Sequences_to_Compare[i].Get_Sequence();

    for(int i=0;i<Output.Patterns_to_Compare.size();i++)
    {
        Matched.push_back(std::vector<int>());

        std::vector<int> Sequence_Map = Output.Patterns_to_Compare[i].Get_Max_Sequence_Map();
        int Sequence_Length = Output.Patterns_to_Compare[i].Get_Max_Sequence_Key().size();

        for(int j=0;j<Sequence_Map.size();j++)
        {
            while(Sequence_Map[j]>Matched[i].size())
                Matched[i].push_back(-1);

            for(int k=0;k<Sequence_Length;k++)
                Matched[i].push_back(j);
        }

        while(Matched[i].size()<Element_Sequence[i].size())
            Matched[i].push_back(-1);
    }

    std::vector<boost_matrix> Cycle_Stats = Output.Extrema_to_Compare.Get_Cycle_Stats();

    for(int i=0;i<max_n_Rows;i++)
    {
        for(int j=0;j<Output.Patterns_to_Compare.size();j++)
        {
            if(i>=Element_Sequence[j].size())
            {
                Diagnostic_Output<<"-1,-1,-1,-1,-1,-1,";
            }
            else
            {
                Diagnostic_Output<<Matched[j][i]<<","<<Element_Sequence[j][i][0]<<","<<Cycle_Stats[j][i][0]<<","<<Cycle_Stats[j][i][1]<<","<<Cycle_Stats[j][i][2]<<","<<Cycle_Stats[j][i][3]<<",";
            }
        }

        Diagnostic_Output<<std::endl;
    }

    Diagnostic_Output.close();
}

bool Measure_Metacommunity_Dynamics(int nSpecies,
                                    double Transient_Test_Time,
                                    std::pair<int,int> Test_Interval,
                                    boost_matrix &Simulation_Data,
                                    bool Solver_Feas_Flag,
                                    double Extinction_Threshold,
                                    std::string Write_Path,
                                    const boost_matrix &Spatial_Structure,
                                    const Simulation_Properties &Sim_Props,
                                    bool Measure_if_Nonstationary,
                                    Dynamics &Output)
{
    int nPatches = (Simulation_Data.shape()[1]-1)/nSpecies;

    Output.Reset_Dynamics(nSpecies,nPatches);

    //Test for feasibility
    Output.Species_Extinct = std::vector<int>(nSpecies,1);

    for(unsigned i = 0;(i+1)<Simulation_Data.shape()[1];i++)
    {
        if(Simulation_Data[Simulation_Data.shape()[0]-1][i+1]>Extinction_Threshold)
            Output.Species_Extinct[i%nSpecies]=0;
    }

    if(*std::max_element(Output.Species_Extinct.begin(),Output.Species_Extinct.end())>0)
        Output.Metacommunity_Feasible = 0;
    else
        Output.Metacommunity_Feasible = 1;

    if(Output.Metacommunity_Feasible==0) //Exit if non-feasible
        return true;

    if(!Solver_Feas_Flag)
    {
        Output.Metacommunity_Solver_Error = 1;
        return true;
    }

    //Test for everything else
    //std::vector< std::vector<int> > Eqb_Synch_Cluster_ID_Species;
    //std::vector< std::vector<int> > Eqb_Synch_Amp_Cluster_ID_Species;

    for(int i=0;i<nPatches;i++) //Start Eqb_Minima with final observations of first species from Simulation_Data
        Output.Community_Eqb_Minima[i] = Simulation_Data[Test_Interval.second][i*nSpecies+1];

    std::vector< std::vector<double> > Patch_Common_Phase_Diff(nPatches,std::vector<double>(nSpecies,0));
    std::vector< std::vector<double> > Patch_Common_Phase_Diff_SD(nPatches,std::vector<double>(nSpecies,0));

    std::vector< std::vector<double> > Patch_Avg_Phase_Diff(nPatches,std::vector<double>(nSpecies,0));

    std::vector< std::vector<double> > Patch_Species_measured(nPatches,std::vector<double>(nSpecies,0));

    for(int i=0;i<nSpecies;i++)
    {
        std::pair<boost_matrix_1d_view,boost_matrix_2d_view> Species_Time_series = Subset_Time_Series(Simulation_Data,Test_Interval,i,nSpecies);
        
        Measure_Time_Series_Properties(Species_Time_series.second,i,Output);

        Find_Cycles Species_Cycles(Species_Time_series,
                                   Sim_Props.Simulation_Resolution,
                                   Sim_Props.Cycle_Params,
                                   Measure_if_Nonstationary);

        //Output.Species_Extrema_CV[i] = Species_Cycles.Extrema_CV;
        //Output.Species_Extrema_Center_Fit[i] = Species_Cycles.Center_Fit;

        //Output.Metacommunity_Avg_Extrema_CV += Species_Cycles.Extrema_CV;
        //Output.Metacommunity_Avg_Extrema_Center_Fit += Species_Cycles.Center_Fit;

        //Output.Species_Stationary[i] = (int)Species_Cycles.Stationary;
        //Output.Metacommunity_Stationary += (double)Species_Cycles.Stationary/(double)nSpecies;

        if(!Species_Cycles.Stationary&&!Measure_if_Nonstationary) //Exit if non-stationary and you're not forcing measures
            return false;

        if(nPatches>1)
        {
            for(int j=0;j<nPatches;j++)
            {
                //Check if Population is Phase Locked
                if(Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Frequency>=Sim_Props.Phase_Lock_Threshold)
                    Output.Species_Phase_Locked[i]+= 1/((double)nPatches-1);
            }

            if(Output.Species_Phase_Locked[i]<Sim_Props.Cycle_Params.Proportion_Series_for_Stationary&&!Measure_if_Nonstationary) //Exit if not phase locked and you're not forcing measures
                return false;
        }
        else
        {
            Output.Species_Phase_Locked[i]=1;
        }

        if(Sim_Props.Write_Measure_Diagnostics)
        {
            Write_Measure_Diagnostics(Write_Path+"_Sp"+std::to_string(static_cast<long long>(i)),Species_Cycles);
        }

        std::vector<double> Common_Phase_Diff;
        std::vector<double> Common_Phase_Diff_SD;

        for(int j=0;j<nPatches;j++)
        {
            if(Species_Cycles.Max_Cycle_Details[j].Relative_Duration>=Sim_Props.Cycle_Params.Stationary_Cycle_Frequency_Threshold)
                Output.Species_Stationary[i]+= 1/(double)nPatches;

            Output.Species_Fixed_Point[i] += (double)Species_Cycles.Fixed_Point[j]/(double)nPatches;

            Output.Population_Metrics[i][j] = Species_Cycles.Max_Cycle_Details[j];

            if(Species_Cycles.Max_Cycle_Details[j].Maxima.size()>0)
            {
                Output.Species_Prop_Measured[i]++;

                Output.Species_Prop_Eqb_Cycle[i] += Species_Cycles.Max_Cycle_Details[j].Relative_Duration;

                Output.Species_Eqb_Cycle_N_Extrema[i] += Species_Cycles.Max_Cycle_Details[j].Maxima.size();

                Output.Species_Avg_Max_Amplitude[i] += Species_Cycles.Max_Cycle_Details[j].Max_Amplitude;
                Output.Species_Avg_Eqb[i] += Species_Cycles.Max_Cycle_Details[j].Average;

                Output.Species_Avg_Period[i] += Species_Cycles.Max_Cycle_Details[j].Period;

                for(int k=0;k<Species_Cycles.Max_Cycle_Details[j].Minima.size();k++)
                {
                    if(Species_Cycles.Max_Cycle_Details[j].Minima[k]>Output.Species_Max_Eqb_Minima[i])
                        Output.Species_Max_Eqb_Minima[i]=Species_Cycles.Max_Cycle_Details[j].Minima[k];

                    if(Species_Cycles.Max_Cycle_Details[j].Minima[k]<Output.Community_Eqb_Minima[j])
                        Output.Community_Eqb_Minima[j]=Species_Cycles.Max_Cycle_Details[j].Minima[k];
                }

                Output.Species_Prop_Common_Phase_Difference[i] +=  Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Frequency;

                //if(Species_Cycles.Base_Key==j)
                    //Output.Species_Prop_Common_Phase_Difference[i] += 1; //0 usually for Base_Key, made 1 to ease comparison

                Output.Species_Phase_Difference_Cycle_N_Diffs[i] += Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff.size();

                Output.Species_Avg_Common_Phase_Difference[i] +=  Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Average;
                Output.Species_Avg_SD_Common_Phase_Difference[i] +=  Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_SD;


                Output.Species_Avg_Phase_Difference[i] +=  Species_Cycles.Max_Cycle_Details[j].Average_Phase_Diff;

                Common_Phase_Diff.push_back(Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Average);
                Common_Phase_Diff_SD.push_back(Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_SD);

                //Community-sorted stats

                 Output.Community_Prop_Measured[j]++;

                Output.Community_Prop_Eqb_Cycle[j] += Species_Cycles.Max_Cycle_Details[j].Relative_Duration;

                Output.Community_Eqb_N_Cycle_Extrema[j] += Species_Cycles.Max_Cycle_Details[j].Maxima.size();

                Output.Community_Avg_Max_Amplitude[j] += Species_Cycles.Max_Cycle_Details[j].Max_Amplitude;
                Output.Community_Avg_Eqb[j] += Species_Cycles.Max_Cycle_Details[j].Average;

                Output.Community_Avg_Period[j] += Species_Cycles.Max_Cycle_Details[j].Period;
                Output.Community_Prop_Common_Phase_Difference[j] += Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Frequency;

                if(Species_Cycles.Base_Key==j)
                    Output.Community_Prop_Common_Phase_Difference[j] += 1; //0 usually for Base_Key, made 1 to ease comparison

                Output.Community_Phase_Difference_Cycle_N_Diffs[j] += Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff.size();

                Output.Community_Avg_Common_Phase_Difference[j] += (Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Average-Species_Cycles.Max_Cycle_Details[0].Common_Phase_Diff_Average);
                Output.Community_Avg_SD_Common_Phase_Difference[j] += (Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_SD-Species_Cycles.Max_Cycle_Details[0].Common_Phase_Diff_SD);
                Output.Community_Avg_Phase_Difference[j] += (Species_Cycles.Max_Cycle_Details[j].Average_Phase_Diff-Species_Cycles.Max_Cycle_Details[0].Average_Phase_Diff);

                Patch_Common_Phase_Diff[j][i] = Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Average-Species_Cycles.Max_Cycle_Details[0].Common_Phase_Diff_Average;
                Patch_Common_Phase_Diff_SD[j][i] = Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_SD-Species_Cycles.Max_Cycle_Details[0].Common_Phase_Diff_SD;

                Patch_Avg_Phase_Diff[j][i] = Species_Cycles.Max_Cycle_Details[j].Average_Phase_Diff-Species_Cycles.Max_Cycle_Details[0].Average_Phase_Diff;

                Patch_Species_measured[j][i] = 1;
            }

            Output.Population_Eqb_Synch_Cluster_ID[i][j]=j;
            Output.Population_Eqb_Synch_Amp_Cluster_ID[i][j]=j;
        }

        if(Output.Species_Prop_Measured[i]>0)
        {

            Output.Species_Prop_Eqb_Cycle[i] /=  Output.Species_Prop_Measured[i];

            Output.Species_Eqb_Cycle_N_Extrema[i] /= Output.Species_Prop_Measured[i];

            Output.Species_Avg_Max_Amplitude[i] /=  Output.Species_Prop_Measured[i];
            Output.Species_Avg_Eqb[i] /=  Output.Species_Prop_Measured[i];

            Output.Species_Avg_Period[i] /=  Output.Species_Prop_Measured[i];

            if(Output.Species_Prop_Measured[i]>1)
            {
                //Phase differences are all compared to Base key, so total is n-1
                //All metrics will be 0 regardless if n<=1

                Output.Species_Prop_Common_Phase_Difference[i] /=  Output.Species_Prop_Measured[i]-1;

                Output.Species_Phase_Difference_Cycle_N_Diffs[i] /= Output.Species_Prop_Measured[i]-1;

                Output.Species_Avg_Common_Phase_Difference[i] /=  Output.Species_Prop_Measured[i]-1;
                Output.Species_Avg_SD_Common_Phase_Difference[i] /=  Output.Species_Prop_Measured[i]-1;

                Output.Species_Avg_Phase_Difference[i] /=  Output.Species_Prop_Measured[i]-1;
            }

            Output.Species_Range_Common_Phase_Difference[i] = *std::max_element(Common_Phase_Diff.begin(),Common_Phase_Diff.end())-*std::min_element(Common_Phase_Diff.begin(),Common_Phase_Diff.end());

            for(int j=0;j<nPatches;j++)
            {
                if(Patch_Species_measured[j][i]>0)
                {
                    Output.Species_CV_Max_Amplitude[i] += pow(Species_Cycles.Max_Cycle_Details[j].Max_Amplitude-Output.Species_Avg_Max_Amplitude[i],2);
                    Output.Species_CV_Eqb[i] += pow(Species_Cycles.Max_Cycle_Details[j].Average-Output.Species_Avg_Eqb[i],2);

                    Output.Species_CV_Period[i] += pow(Species_Cycles.Max_Cycle_Details[j].Period-Output.Species_Avg_Period[i],2);
                    Output.Species_SD_Avg_Phase_Difference[i] += pow(Species_Cycles.Max_Cycle_Details[j].Average_Phase_Diff-Output.Species_Avg_Phase_Difference[i],2);
                }
            }

            Output.Species_CV_Max_Amplitude[i] /=  Output.Species_Prop_Measured[i];
            Output.Species_CV_Max_Amplitude[i]=pow(Output.Species_CV_Max_Amplitude[i],.5);
            Output.Species_CV_Max_Amplitude[i] /= Output.Species_Avg_Max_Amplitude[i];

            Output.Species_CV_Eqb[i] /=  Output.Species_Prop_Measured[i];
            Output.Species_CV_Eqb[i]=pow(Output.Species_CV_Eqb[i],.5);
            Output.Species_CV_Eqb[i] /= Output.Species_Avg_Eqb[i];

            Output.Species_CV_Period[i] /=  Output.Species_Prop_Measured[i];
            Output.Species_CV_Period[i]=pow(Output.Species_CV_Period[i],.5);
            Output.Species_CV_Period[i] /= Output.Species_Avg_Period[i];

            if(Output.Species_Prop_Measured[i]>1) //As above, Phase diff metrics will be 0 if n<=1, only meaningful if n>1
                Output.Species_SD_Avg_Phase_Difference[i] /=  Output.Species_Prop_Measured[i]-1;

            Output.Species_SD_Avg_Phase_Difference[i]=pow(Output.Species_SD_Avg_Phase_Difference[i],.5);

        }

        //Synch Clusters calc
        for(int j=0;j<nPatches;j++)
        {
            for(int k = j+1;k<nPatches;k++)
            {
                bool Exact_Time_Series_Match = true;

                for(int l=0;l<20;l++) //match 20 time steps, super magic-numbery but whatever I'm tired
                {
                    if(Species_Time_series.first.shape()[0]-1-l>=0)
                    {
                        if(std::abs(Species_Time_series.second[Species_Time_series.second.shape()[0]-1-l][j]-Species_Time_series.second[Species_Time_series.second.shape()[0]-1-l][k])>=Sim_Props.Cycle_Params.Difference_Error_Tolerance)
                        {
                            Exact_Time_Series_Match = false;
                            break;
                        }
                    }
                    else
                        break;
                }

                if(Exact_Time_Series_Match)
                    Output.Population_Eqb_Synch_Amp_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Amp_Cluster_ID[i][j];


                if(Species_Cycles.Fixed_Point[j]==1&&Species_Cycles.Fixed_Point[k]==1)
                {
                    if(std::abs(Species_Time_series.second[Species_Time_series.second.shape()[0]-1][j]-Species_Time_series.second[Species_Time_series.second.shape()[0]-1][k])<Sim_Props.Cycle_Params.Difference_Error_Tolerance)
                    {
                        Output.Population_Eqb_Synch_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Cluster_ID[i][j];
                        //Output.Population_Eqb_Synch_Amp_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Amp_Cluster_ID[i][j];
                    }

                }else{
                    if(Species_Cycles.Max_Cycle_Details[j].Relative_Duration>=Sim_Props.Cycle_Params.Stationary_Cycle_Frequency_Threshold&&
                            Species_Cycles.Max_Cycle_Details[k].Relative_Duration>=Sim_Props.Cycle_Params.Stationary_Cycle_Frequency_Threshold&&
                            std::abs(Species_Cycles.Max_Cycle_Details[j].Period-Species_Cycles.Max_Cycle_Details[k].Period)<Sim_Props.Cycle_Params.Difference_Error_Tolerance)
                    {
                        if(std::abs(Species_Cycles.Max_Cycle_Details[j].Common_Phase_Diff_Average-Species_Cycles.Max_Cycle_Details[k].Common_Phase_Diff_Average)<Sim_Props.Cycle_Params.Difference_Error_Tolerance)
                        {
                            Output.Population_Eqb_Synch_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Cluster_ID[i][j];

//                            if(std::abs(Species_Cycles.Max_Cycle_Details[j].Average_Phase_Diff-Species_Cycles.Max_Cycle_Details[k].Average_Phase_Diff)<Sim_Props.Cycle_Params.Difference_Error_Tolerance&&
//                               std::abs(Species_Cycles.Max_Cycle_Details[j].SD_Phase_Diff-Species_Cycles.Max_Cycle_Details[k].SD_Phase_Diff)<Sim_Props.Cycle_Params.Difference_Error_Tolerance)
//                            {
//                                Output.Population_Eqb_Synch_Amp_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Amp_Cluster_ID[i][j];
//                            }
                        }
                    }
                    else
                    {
                        double sum_x = 0;
                        double sum_y = 0;

                        for(int l=0;l<Species_Time_series.second.shape()[0];l++)
                        {
                            sum_x = sum_x + Species_Time_series.second[l][j];
                            sum_y = sum_y + Species_Time_series.second[l][k];
                        }

                        double nPairs = Species_Time_series.second.shape()[0]-1;

                        double mean_x = sum_x/nPairs;
                        double mean_y = sum_y/nPairs;

                        double sum_x_2_dev = 0;
                        double sum_y_2_dev = 0;
                        double sum_xy_dev = 0;

                        for(int l=0;l<Species_Time_series.second.shape()[0];l++)
                        {
                            sum_x_2_dev = sum_x_2_dev + pow(Species_Time_series.second[l][j]-mean_x,2.0);
                            sum_y_2_dev = sum_y_2_dev + pow(Species_Time_series.second[l][k]-mean_y,2.0);
                            sum_xy_dev = sum_xy_dev + (Species_Time_series.second[l][j]-mean_x)*(Species_Time_series.second[l][k]-mean_y);
                        }

                        double Clust_test = sum_xy_dev/sqrt(sum_x_2_dev*sum_y_2_dev);

                        if(Clust_test>.999)
                        {
                            Output.Population_Eqb_Synch_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Cluster_ID[i][j];
                            //Output.Population_Eqb_Synch_Amp_Cluster_ID[i][k]=Output.Population_Eqb_Synch_Amp_Cluster_ID[i][j];
                        }
                    }
                }
            }
        }

        Output.Species_Eqb_Synch_Clusters[i] = count_distinct_abs(Output.Population_Eqb_Synch_Cluster_ID[i]);
        Output.Species_Eqb_Synch_Amp_Clusters[i] = count_distinct_abs(Output.Population_Eqb_Synch_Amp_Cluster_ID[i]);

        if(Output.Species_Eqb_Synch_Amp_Clusters[i]>1)
            Output.Species_Asynchronous[i] = 1;
        //End Synch Clusters calc

        Output.Metacommunity_Stationary += Output.Species_Stationary[i];
        Output.Metacommunity_Phase_Locked += Output.Species_Phase_Locked[i];

        Output.Metacommunity_Fixed_Point += Output.Species_Fixed_Point[i];

        if(Output.Species_Asynchronous[i]>0)
            Output.Metacommunity_Asynchronous = 1;

        Output.Metacommunity_Prop_Eqb_Cycle += Output.Species_Prop_Eqb_Cycle[i];

        Output.Metacommunity_Eqb_Cycle_N_Extrema += Output.Species_Eqb_Cycle_N_Extrema[i];

        Output.Metacommunity_Avg_Max_Amplitude += Output.Species_Avg_Max_Amplitude[i];
        Output.Metacommunity_Avg_Eqb += Output.Species_Avg_Eqb[i];

        Output.Metacommunity_Avg_CV_Max_Amplitude += Output.Species_CV_Max_Amplitude[i];
        Output.Metacommunity_Avg_CV_Eqb += Output.Species_CV_Eqb[i];

        Output.Metacommunity_Avg_Period += Output.Species_Avg_Period[i];
        Output.Metacommunity_Avg_CV_Period += Output.Species_CV_Period[i];

        Output.Metacommunity_Prop_Common_Phase_Difference += Output.Species_Prop_Common_Phase_Difference[i];

        Output.Metacommunity_Phase_Difference_Cycle_N_Diffs += Output.Species_Phase_Difference_Cycle_N_Diffs[i];

        Output.Metacommunity_Avg_Common_Phase_Difference += Output.Species_Avg_Common_Phase_Difference[i];
        Output.Metacommunity_Avg_SD_Common_Phase_Difference += Output.Species_Avg_SD_Common_Phase_Difference[i];

        Output.Metacommunity_Avg_Phase_Difference += Output.Species_Avg_Phase_Difference[i];

        Output.Metacommunity_Avg_SD_Avg_Phase_Difference += Output.Species_SD_Avg_Phase_Difference[i];
        Output.Metacommunity_Avg_Range_Common_Phase_Difference += Output.Species_Range_Common_Phase_Difference[i];

        Output.Species_Prop_Measured[i] /= (double)nPatches;

        Output.Metacommunity_Prop_Measured += Output.Species_Prop_Measured[i];
    }

    if(Output.Metacommunity_Stationary>=Sim_Props.Cycle_Params.Proportion_Series_for_Stationary&&Output.Metacommunity_Phase_Locked>=Sim_Props.Cycle_Params.Proportion_Series_for_Stationary)
        Output.Metacommunity_Transient_Length = Transient_Test_Time;

    Output.Metacommunity_Stationary /= (double)nSpecies;
    Output.Metacommunity_Phase_Locked /= (double)nSpecies;

    Output.Metacommunity_Fixed_Point /= (double)nSpecies;

    Output.Metacommunity_Prop_Eqb_Cycle /= (double)nSpecies;

    Output.Metacommunity_Eqb_Cycle_N_Extrema /= (double)nSpecies;

    Output.Metacommunity_Avg_Max_Amplitude /= (double)nSpecies;
    Output.Metacommunity_Avg_Eqb /= (double)nSpecies;

    Output.Metacommunity_Avg_CV_Max_Amplitude /= (double)nSpecies;
    Output.Metacommunity_Avg_CV_Eqb /= (double)nSpecies;

    Output.Metacommunity_Avg_Period /= (double)nSpecies;
    Output.Metacommunity_Avg_CV_Period /= (double)nSpecies;

    Output.Metacommunity_Prop_Common_Phase_Difference /= (double)nSpecies;

    Output.Metacommunity_Phase_Difference_Cycle_N_Diffs /= (double)nSpecies;

    Output.Metacommunity_Avg_Common_Phase_Difference /= (double)nSpecies;
    Output.Metacommunity_Avg_SD_Common_Phase_Difference /= (double)nSpecies;

    Output.Metacommunity_Avg_Phase_Difference /= (double)nSpecies;

    Output.Metacommunity_Avg_SD_Avg_Phase_Difference /= (double)nSpecies;
    Output.Metacommunity_Avg_Range_Common_Phase_Difference /= (double)nSpecies;

    Output.Metacommunity_Prop_Measured /= (double)nSpecies;

    //Output.Metacommunity_Avg_Extrema_CV /= (double)nSpecies;
    //Output.Metacommunity_Avg_Extrema_Center_Fit /= (double)nSpecies;

    std::vector<int>::iterator Max_Synch_Clusters = std::max_element(Output.Species_Eqb_Synch_Clusters.begin(),Output.Species_Eqb_Synch_Clusters.end());
    auto Most_Synch_Clusters = std::distance(Output.Species_Eqb_Synch_Clusters.begin(),Max_Synch_Clusters);
    Output.Metacommunity_Eqb_Synch_Clusters = *Max_Synch_Clusters;
    Output.Community_Eqb_Synch_Cluster_ID = Output.Population_Eqb_Synch_Cluster_ID[(int)Most_Synch_Clusters];

    std::vector<int>::iterator Max_Synch_Amp_Clusters = std::max_element(Output.Species_Eqb_Synch_Amp_Clusters.begin(),Output.Species_Eqb_Synch_Amp_Clusters.end());
    auto Most_Synch_Amp_Clusters = std::distance(Output.Species_Eqb_Synch_Amp_Clusters.begin(),Max_Synch_Amp_Clusters);
    Output.Metacommunity_Eqb_Synch_Amp_Clusters = *Max_Synch_Amp_Clusters;
    Output.Community_Eqb_Synch_Amp_Cluster_ID = Output.Population_Eqb_Synch_Amp_Cluster_ID[(int)Most_Synch_Amp_Clusters];

    Output.Metacommunity_Eqb_Min_Max_Minima = *std::min_element(Output.Species_Max_Eqb_Minima.begin(),Output.Species_Max_Eqb_Minima.end());

    for(int i=0;i<nPatches;i++)
    {
        if( Output.Community_Prop_Measured[i]>0)
        {

            Output.Community_Prop_Eqb_Cycle[i] /=  Output.Community_Prop_Measured[i];

            Output.Community_Eqb_N_Cycle_Extrema[i] /= Output.Community_Prop_Measured[i];

            Output.Community_Avg_Max_Amplitude[i] /=  Output.Community_Prop_Measured[i];
            Output.Community_Avg_Eqb[i] /=  Output.Community_Prop_Measured[i];

            Output.Community_Avg_Period[i] /=  Output.Community_Prop_Measured[i];

            Output.Community_Prop_Common_Phase_Difference[i] /=  Output.Community_Prop_Measured[i];

            Output.Community_Phase_Difference_Cycle_N_Diffs[i] /= Output.Community_Phase_Difference_Cycle_N_Diffs[i];

            Output.Community_Avg_Common_Phase_Difference[i] /=  Output.Community_Prop_Measured[i];
            Output.Community_Avg_SD_Common_Phase_Difference[i] /=  Output.Community_Prop_Measured[i];

            Output.Community_Avg_Phase_Difference[i] /=  Output.Community_Prop_Measured[i];

            std::vector<double> Patch_Measured_Common_Phase_Diff;

            for(int j=0;j<nSpecies;j++)
            {
                if(Patch_Species_measured[i][j]>0)
                {
                    Output.Community_SD_Avg_Phase_Difference[i] += pow(Patch_Avg_Phase_Diff[i][j]-Output.Community_Avg_Phase_Difference[i],2);

                    Patch_Measured_Common_Phase_Diff.push_back(Patch_Common_Phase_Diff[i][j]);
                }
            }

            Output.Community_SD_Avg_Phase_Difference[i] /=  Output.Community_Prop_Measured[i];
            Output.Community_SD_Avg_Phase_Difference[i] = pow(Output.Community_SD_Avg_Phase_Difference[i],.5);

            if(Patch_Measured_Common_Phase_Diff.size()>0)
            {
                Output.Community_Range_Common_Phase_Difference[i] = *std::max_element(Patch_Measured_Common_Phase_Diff.begin(),Patch_Measured_Common_Phase_Diff.end())
                        -*std::min_element(Patch_Measured_Common_Phase_Diff.begin(),Patch_Measured_Common_Phase_Diff.end());
            }
        }

        Output.Community_Prop_Measured[i] /= (double)nSpecies;
    }

    if(Output.Metacommunity_Eqb_Synch_Amp_Clusters>1)
    {
        Quotient_Structure_Metadata Quotient_Structure_Measures = Calculate_Quotient_Structure_Metadata(Spatial_Structure,Output.Community_Eqb_Synch_Amp_Cluster_ID);

        Output.Metacommunity_Quotient_Connec_In = Quotient_Structure_Measures.Average_In_Degree;
        Output.Metacommunity_Quotient_Connec_Out = Quotient_Structure_Measures.Average_Out_Degree;
        Output.Metacommunity_Quotient_Skew_Connec_In = Quotient_Structure_Measures.In_Degree_Skewness;
        Output.Metacommunity_Quotient_Skew_Connec_Out = Quotient_Structure_Measures.Out_Degree_Skewness;

        Output.Metacommunity_Quotient_Clustering = Quotient_Structure_Measures.Clustering_Coefficient;
        Output.Metacommunity_Quotient_Mean_Path_Length = Quotient_Structure_Measures.Average_Path_Length;
        Output.Metacommunity_Quotient_Eigenratio = Quotient_Structure_Measures.Eigenratio;
    }

    return true;
}

//Output elements for i_Species must be set first
void Measure_Time_Series_Properties(const boost_matrix_2d_view &Species_Time_Series, int i_Species, Dynamics &Output)
{    
    Output.Population_Time_Series_Average[i_Species].resize(Species_Time_Series.shape()[1]);
    Output.Population_Time_Series_CV[i_Species].resize(Species_Time_Series.shape()[1]);
    Output.Population_Time_Series_Skewness[i_Species].resize(Species_Time_Series.shape()[1]);
    Output.Population_Time_Series_Minimum[i_Species].resize(Species_Time_Series.shape()[1]);
    Output.Population_Time_Series_Maximum[i_Species].resize(Species_Time_Series.shape()[1]);
    
    std::fill(Output.Population_Time_Series_Average[i_Species].begin(),Output.Population_Time_Series_Average[i_Species].end(),0);
    std::fill(Output.Population_Time_Series_CV[i_Species].begin(),Output.Population_Time_Series_CV[i_Species].end(),0);
    std::fill(Output.Population_Time_Series_Skewness[i_Species].begin(),Output.Population_Time_Series_Skewness[i_Species].end(),0);
    std::fill(Output.Population_Time_Series_Minimum[i_Species].begin(),Output.Population_Time_Series_Minimum[i_Species].end(),0);
    std::fill(Output.Population_Time_Series_Maximum[i_Species].begin(),Output.Population_Time_Series_Maximum[i_Species].end(),0);
    
    for(int j=0;j<Species_Time_Series.shape()[1];j++)
    {
        Output.Population_Time_Series_Minimum[i_Species][j] = Species_Time_Series[0][j];
        Output.Population_Time_Series_Maximum[i_Species][j] = Species_Time_Series[0][j];
        
    }
    
    std::vector<double> Species_Total_Time_Series;

    Output.Species_Total_Average[i_Species] = 0;
    Output.Species_Total_CV[i_Species] = 0;
    
    for(int i=0;i<Species_Time_Series.shape()[0];i++)
    {
        Species_Total_Time_Series.push_back(0);
        
        for(int j=0;j<Species_Time_Series.shape()[1];j++)
        {
            Species_Total_Time_Series[i] += Species_Time_Series[i][j];
            
            Output.Population_Time_Series_Average[i_Species][j] += Species_Time_Series[i][j];
        
            if(Species_Time_Series[i][j]<Output.Population_Time_Series_Minimum[i_Species][j])
                Output.Population_Time_Series_Minimum[i_Species][j]=Species_Time_Series[i][j];
        
            if(Species_Time_Series[i][j]>Output.Population_Time_Series_Maximum[i_Species][j])
                Output.Population_Time_Series_Maximum[i_Species][j]=Species_Time_Series[i][j];
        }
        
        if(i==0)
        {
            Output.Species_Total_Minimum[i_Species] = Species_Total_Time_Series[0];
            Output.Species_Total_Maximum[i_Species] = Species_Total_Time_Series[0];
        }

        Output.Species_Total_Average[i_Species]+=Species_Total_Time_Series[i];
        
        if(Species_Total_Time_Series[i]<Output.Species_Total_Minimum[i_Species])
            Output.Species_Total_Minimum[i_Species]=Species_Total_Time_Series[i];
        if(Species_Total_Time_Series[i]>Output.Species_Total_Maximum[i_Species])
            Output.Species_Total_Maximum[i_Species]=Species_Total_Time_Series[i];
        
    }
    
    for(int j=0;j<Species_Time_Series.shape()[1];j++)
    {
        Output.Population_Time_Series_Average[i_Species][j] /= (double)Species_Time_Series.shape()[0];
    }
    
    Output.Species_Total_Average[i_Species] /= (double)Species_Time_Series.shape()[0];
    
    for(int i=0;i<Species_Time_Series.shape()[0];i++)
    {
        for(int j=0;j<Species_Time_Series.shape()[1];j++)
        {
            Output.Population_Time_Series_CV[i_Species][j] += pow(Species_Time_Series[i][j]-Output.Population_Time_Series_Average[i_Species][j],2);
            Output.Population_Time_Series_Skewness[i_Species][j] += pow(Species_Time_Series[i][j]-Output.Population_Time_Series_Average[i_Species][j],3);
        }
        
        Output.Species_Total_CV[i_Species] += pow(Species_Total_Time_Series[i]-Output.Species_Total_Average[i_Species],2);
    }
    
    for(int j=0;j<Species_Time_Series.shape()[1];j++)
    {
        Output.Population_Time_Series_CV[i_Species][j] /= (double)Species_Time_Series.shape()[0];
        Output.Population_Time_Series_CV[i_Species][j] = pow(Output.Population_Time_Series_CV[i_Species][j],.5); //this is SD at this point
        
        Output.Population_Time_Series_Skewness[i_Species][j] /= (double)Species_Time_Series.shape()[0];
        Output.Population_Time_Series_Skewness[i_Species][j] /= pow(Output.Population_Time_Series_CV[i_Species][j],3); //skewness / SD^3
        
        Output.Population_Time_Series_CV[i_Species][j] /= Output.Population_Time_Series_Average[i_Species][j]; //turn SD into CV
    }
    
    Output.Species_Total_CV[i_Species] /= (double)Species_Time_Series.shape()[0];
    Output.Species_Total_CV[i_Species] = pow(Output.Species_Total_CV[i_Species],.5);
   
}


Dynamics::Dynamics()
{
    Community_Dynamics_Labels.push_back("Feasible");
    Community_Dynamics_Labels.push_back("Solver_Error");
    Community_Dynamics_Labels.push_back("Stationary");
    Community_Dynamics_Labels.push_back("Phase_Locked");
    Community_Dynamics_Labels.push_back("Fixed_Point");
    Community_Dynamics_Labels.push_back("Prop_Measured");
    Community_Dynamics_Labels.push_back("Transient_Length");
    Community_Dynamics_Labels.push_back("Asynchronous");
    Community_Dynamics_Labels.push_back("Eqb_Synch_Clusters");
    Community_Dynamics_Labels.push_back("Eqb_Synch_Amp_Clusters");
    Community_Dynamics_Labels.push_back("Prop_Eqb_Cycle");
    Community_Dynamics_Labels.push_back("Eqb_Cycle_N_Extrema");
    Community_Dynamics_Labels.push_back("Eqb_Min_Max_Minima");
    Community_Dynamics_Labels.push_back("Avg_Max_Amplitude");
    Community_Dynamics_Labels.push_back("Avg_Eqb");
    Community_Dynamics_Labels.push_back("Avg_CV_Max_Amplitude");
    Community_Dynamics_Labels.push_back("Avg_CV_Eqb");
    Community_Dynamics_Labels.push_back("Avg_Period");
    Community_Dynamics_Labels.push_back("Avg_CV_Period");
    Community_Dynamics_Labels.push_back("Prop_Common_Phase_Difference");
    Community_Dynamics_Labels.push_back("Phase_Difference_Cycle_N_Diffs");
    Community_Dynamics_Labels.push_back("Avg_Common_Phase_Difference");
    Community_Dynamics_Labels.push_back("Avg_SD_Common_Phase_Difference");
    Community_Dynamics_Labels.push_back("Avg_Phase_Difference");
    Community_Dynamics_Labels.push_back("Avg_Range_Common_Phase_Difference");
    Community_Dynamics_Labels.push_back("Avg_SD_Avg_Phase_Difference");
    Community_Dynamics_Labels.push_back("Quotient_Connec_In");
    Community_Dynamics_Labels.push_back("Quotient_Connec_Out");
    Community_Dynamics_Labels.push_back("Quotient_Skew_Connec_In");
    Community_Dynamics_Labels.push_back("Quotient_Skew_Connec_Out");
    Community_Dynamics_Labels.push_back("Quotient_Clustering");
    Community_Dynamics_Labels.push_back("Quotient_Mean_Path_Length");
    Community_Dynamics_Labels.push_back("Quotient_Eigenratio");

    Species_Dynamics_Labels.push_back("Extinct");
    Species_Dynamics_Labels.push_back("Stationary");
    Species_Dynamics_Labels.push_back("Phase_Locked");
    Species_Dynamics_Labels.push_back("Fixed_Point");
    Species_Dynamics_Labels.push_back("Prop_Measured");
    Species_Dynamics_Labels.push_back("Asynchronous");
    Species_Dynamics_Labels.push_back("Eqb_Synch_Clusters");
    Species_Dynamics_Labels.push_back("Eqb_Synch_Amp_Clusters");
    Species_Dynamics_Labels.push_back("Prop_Eqb_Cycle");
    Species_Dynamics_Labels.push_back("Eqb_Cycle_N_Extrema");
    Species_Dynamics_Labels.push_back("Eqb_Max_Minima");
    Species_Dynamics_Labels.push_back("Avg_Max_Amplitude");
    Species_Dynamics_Labels.push_back("Avg_Eqb");
    Species_Dynamics_Labels.push_back("CV_Max_Amplitude");
    Species_Dynamics_Labels.push_back("CV_Eqb");
    Species_Dynamics_Labels.push_back("Avg_Period");
    Species_Dynamics_Labels.push_back("CV_Period");
    Species_Dynamics_Labels.push_back("Prop_Common_Phase_Difference");
    Species_Dynamics_Labels.push_back("Phase_Difference_Cycle_N_Diffs");
    Species_Dynamics_Labels.push_back("Avg_Common_Phase_Difference");
    Species_Dynamics_Labels.push_back("Avg_SD_Common_Phase_Difference");
    Species_Dynamics_Labels.push_back("Avg_Phase_Difference");
    Species_Dynamics_Labels.push_back("Range_Common_Phase_Difference");
    Species_Dynamics_Labels.push_back("SD_Avg_Phase_Difference");
    Species_Dynamics_Labels.push_back("Total_Time_Series_Average");
    Species_Dynamics_Labels.push_back("Total_Time_Series_CV");
    Species_Dynamics_Labels.push_back("Total_Time_Series_Minimum");
    Species_Dynamics_Labels.push_back("Total_Time_Series_Maximum");

    Patch_Dynamics_Labels.push_back("Eqb_Synch_Cluster_ID");
    Patch_Dynamics_Labels.push_back("Eqb_Synch_Amp_Cluster_ID");
    Patch_Dynamics_Labels.push_back("Community_Prop_Measured");
    Patch_Dynamics_Labels.push_back("Prop_Eqb_Cycle");
    Patch_Dynamics_Labels.push_back("Eqb_Cycle_N_Extrema");
    Patch_Dynamics_Labels.push_back("Eqb_Minima");
    Patch_Dynamics_Labels.push_back("Avg_Max_Amplitude");
    Patch_Dynamics_Labels.push_back("Avg_Eqb");
    Patch_Dynamics_Labels.push_back("Avg_Period");
    Patch_Dynamics_Labels.push_back("Prop_Common_Phase_Difference");
    Patch_Dynamics_Labels.push_back("Phase_Difference_Cycle_N_Diffs");
    Patch_Dynamics_Labels.push_back("Avg_Common_Phase_Difference");
    Patch_Dynamics_Labels.push_back("Avg_SD_Common_Phase_Difference");
    Patch_Dynamics_Labels.push_back("Avg_Phase_Difference");
    Patch_Dynamics_Labels.push_back("Range_Common_Phase_Difference");
    Patch_Dynamics_Labels.push_back("SD_Avg_Phase_Difference");

    Population_Dynamics_Labels.push_back("Extrema_Variability");
    Population_Dynamics_Labels.push_back("Center_Fit");
    Population_Dynamics_Labels.push_back("Relative_Duration");
    Population_Dynamics_Labels.push_back("Period");
    Population_Dynamics_Labels.push_back("Eqb_Cycle_N_Extrema");
    Population_Dynamics_Labels.push_back("Max_Amplitude");
    Population_Dynamics_Labels.push_back("Eqb_Average");
    Population_Dynamics_Labels.push_back("Eqb_SD");
    Population_Dynamics_Labels.push_back("Eqb_Skewness");
    Population_Dynamics_Labels.push_back("Avg_Common_Phase_Diff");
    Population_Dynamics_Labels.push_back("SD_Common_Phase_Diff");
    Population_Dynamics_Labels.push_back("Common_Phase_Diff_Frequency");
    Population_Dynamics_Labels.push_back("Phase_Difference_Cycle_N_Diffs");
    Population_Dynamics_Labels.push_back("Average_Phase_Diff");
    Population_Dynamics_Labels.push_back("SD_Phase_Diff");
    Population_Dynamics_Labels.push_back("Eqb_Synch_Cluster_ID");
    Population_Dynamics_Labels.push_back("Eqb_Synch_Amp_Cluster_ID");
    Population_Dynamics_Labels.push_back("Time_Series_Average");
    Population_Dynamics_Labels.push_back("Time_Series_CV");
    Population_Dynamics_Labels.push_back("Time_Series_Skewness");
    Population_Dynamics_Labels.push_back("Time_Series_Minimum");
    Population_Dynamics_Labels.push_back("Time_Series_Maximum");


    n_Community_Metrics = 33;
    n_Species_Metrics = 28;
    n_Patch_Metrics = 16;
    n_Population_Metrics = 22;
}

double Dynamics::Get_Community_Measures(int Index) const
{
    switch(Index)
    {
    case 0:
    {
        return (double)Metacommunity_Feasible;
    }
        break;
    case 1:
    {
        return Metacommunity_Solver_Error;
    }
        break;
    case 2:
    {
        return Metacommunity_Stationary;
    }
        break;
    case 3:
    {
        return Metacommunity_Phase_Locked;
    }
        break;
    case 4:
    {
        return Metacommunity_Fixed_Point;
    }
        break;
    case 5:
    {
        return Metacommunity_Prop_Measured;
    }
        break;
    case 6:
    {
        return Metacommunity_Transient_Length;
    }
        break;
    case 7:
    {
        return Metacommunity_Asynchronous;
    }
        break;
    case 8:
    {
        return (double) Metacommunity_Eqb_Synch_Clusters;
    }
        break;
    case 9:
    {
        return (double) Metacommunity_Eqb_Synch_Amp_Clusters;
    }
        break;
    case 10:
    {
        return Metacommunity_Prop_Eqb_Cycle;
    }
        break;
    case 11:
    {
        return Metacommunity_Eqb_Cycle_N_Extrema;
    }
        break;
    case 12:
    {
        return Metacommunity_Eqb_Min_Max_Minima;
    }
        break;
    case 13:
    {
        return Metacommunity_Avg_Max_Amplitude;
    }
        break;
    case 14:
    {
        return Metacommunity_Avg_Eqb;
    }
        break;
    case 15:
    {
        return Metacommunity_Avg_CV_Max_Amplitude;
    }
        break;
    case 16:
    {
        return Metacommunity_Avg_CV_Eqb;
    }
        break;
    case 17:
    {
        return Metacommunity_Avg_Period;
    }
        break;
    case 18:
    {
        return Metacommunity_Avg_CV_Period;
    }
        break;
    case 19:
    {
        return Metacommunity_Prop_Common_Phase_Difference;
    }
        break;
    case 20:
    {
        return Metacommunity_Phase_Difference_Cycle_N_Diffs;
    }
        break;
    case 21:
    {
        return Metacommunity_Avg_Common_Phase_Difference;
    }
        break;
    case 22:
    {
        return Metacommunity_Avg_SD_Common_Phase_Difference;
    }
        break;
    case 23:
    {
        return Metacommunity_Avg_Phase_Difference;
    }
        break;
    case 24:
    {
        return Metacommunity_Avg_Range_Common_Phase_Difference;
    }
        break;
    case 25:
    {
        return Metacommunity_Avg_SD_Avg_Phase_Difference;
    }
        break;
    case 26:
    {
        return Metacommunity_Quotient_Connec_In;
    }
        break;
    case 27:
    {
        return Metacommunity_Quotient_Connec_Out;
    }
        break;
    case 28:
    {
        return Metacommunity_Quotient_Skew_Connec_In;
    }
        break;
    case 29:
    {
        return Metacommunity_Quotient_Skew_Connec_Out;
    }
        break;
    case 30:
    {
        return Metacommunity_Quotient_Clustering;
    }
        break;
    case 31:
    {
        return Metacommunity_Quotient_Mean_Path_Length;
    }
        break;
    case 32:
    {
        return Metacommunity_Quotient_Eigenratio;
    }
        break;
    }
    return 0;
}

std::vector<double> Dynamics::Get_Species_Measures(int Index) const
{
    switch(Index)
    {
    case 0:
    {
        return std::vector<double>(Species_Extinct.begin(),Species_Extinct.end());
    }
        break;
    case 1:
    {
        return Species_Stationary;
    }
        break;
    case 2:
    {
        return Species_Phase_Locked;
    }
        break;
    case 3:
    {
        return Species_Fixed_Point;
    }
        break;
    case 4:
    {
        return Species_Prop_Measured;
    }
        break;
    case 5:
    {
        return Species_Asynchronous;
    }
        break;
    case 6:
    {
        return std::vector<double>(Species_Eqb_Synch_Clusters.begin(),Species_Eqb_Synch_Clusters.end());
    }
        break;
    case 7:
    {
        return std::vector<double>(Species_Eqb_Synch_Amp_Clusters.begin(),Species_Eqb_Synch_Amp_Clusters.end());
    }
        break;
    case 8:
    {
        return Species_Prop_Eqb_Cycle;
    }
        break;
    case 9:
    {
        return Species_Eqb_Cycle_N_Extrema;
    }
        break;
    case 10:
    {
        return Species_Max_Eqb_Minima;
    }
        break;
    case 11:
    {
        return Species_Avg_Max_Amplitude;
    }
        break;
    case 12:
    {
        return Species_Avg_Eqb;
    }
        break;
    case 13:
    {
        return Species_CV_Max_Amplitude;
    }
        break;
    case 14:
    {
        return Species_CV_Eqb;
    }
        break;
    case 15:
    {
        return Species_Avg_Period;
    }
        break;
    case 16:
    {
        return Species_CV_Period;
    }
        break;
    case 17:
    {
        return Species_Prop_Common_Phase_Difference;
    }
        break;
    case 18:
    {
        return Species_Phase_Difference_Cycle_N_Diffs;
    }
        break;
    case 19:
    {
        return Species_Avg_Common_Phase_Difference;
    }
        break;
    case 20:
    {
        return Species_Avg_SD_Common_Phase_Difference;
    }
        break;
    case 21:
    {
        return Species_Avg_Phase_Difference;
    }
        break;
    case 22:
    {
        return Species_Range_Common_Phase_Difference;
    }
        break;
    case 23:
    {
        return Species_SD_Avg_Phase_Difference;
    }
        break;
    case 24:
    {
        return Species_Total_Average;
    }
        break;
    case 25:
    {
        return Species_Total_CV;
    }
        break;
    case 26:
    {
        return Species_Total_Minimum;
    }
        break;
    case 27:
    {
        return Species_Total_Maximum;
    }
        break;
    }

    return std::vector<double>(1,0);
}

std::vector<double> Dynamics::Get_Patch_Measures(int Index) const
{
    switch(Index)
    {
    case 0:
    {
        return std::vector<double>(Community_Eqb_Synch_Cluster_ID.begin(),Community_Eqb_Synch_Cluster_ID.end());
    }
        break;
    case 1:
    {
        return std::vector<double>(Community_Eqb_Synch_Amp_Cluster_ID.begin(),Community_Eqb_Synch_Amp_Cluster_ID.end());
    }
        break;
    case 2:
    {
        return Community_Prop_Measured;
    }
        break;
    case 3:
    {
        return Community_Prop_Eqb_Cycle;
    }
        break;
    case 4:
    {
        return Community_Eqb_N_Cycle_Extrema;
    }
        break;
    case 5:
    {
        return Community_Eqb_Minima;
    }
        break;
    case 6:
    {
        return Community_Avg_Max_Amplitude;
    }
        break;
    case 7:
    {
        return Community_Avg_Eqb;
    }
        break;
    case 8:
    {
        return Community_Avg_Period;
    }
        break;
    case 9:
    {
        return Community_Prop_Common_Phase_Difference;
    }
        break;
    case 10:
    {
        return Community_Phase_Difference_Cycle_N_Diffs;
    }
        break;
    case 11:
    {
        return Community_Avg_Common_Phase_Difference;
    }
        break;
    case 12:
    {
        return Community_Avg_SD_Common_Phase_Difference;
    }
        break;
    case 13:
    {
        return Community_Avg_Phase_Difference;
    }
        break;
    case 14:
    {
        return Community_Range_Common_Phase_Difference;
    }
        break;
    case 15:
    {
        return Community_SD_Avg_Phase_Difference;
    }
        break;
    }
    return std::vector<double>(1,0);
}

std::vector<std::vector<double>> Dynamics::Get_Population_Measures(int Index) const
{
    switch(Index)
    {
    case 0:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Extrema_Variability);
        }

        return to_Return;
    }
    case 1:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Center_Fit);
        }

        return to_Return;
    }
    case 2:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Relative_Duration);
        }

        return to_Return;
    }
        break;
    case 3:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Period);
        }

        return to_Return;
    }
        break;
    case 4:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Maxima.size());
        }

        return to_Return;
    }
        break;
    case 5:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Max_Amplitude);
        }

        return to_Return;
    }
        break;
    case 6:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Average);
        }

        return to_Return;
    }
        break;
    case 7:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].SD);
        }

        return to_Return;
    }
        break;
    case 8:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Skewness);
        }

        return to_Return;
    }
        break;
    case 9:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Common_Phase_Diff_Average);
        }

        return to_Return;
    }
    case 10:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Common_Phase_Diff_SD);
        }

        return to_Return;
    }
        break;
    case 11:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Common_Phase_Diff_Frequency);
        }

        return to_Return;
    }
        break;
    case 12:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Common_Phase_Diff.size());
        }

        return to_Return;
    }
        break;
    case 13:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].Average_Phase_Diff);
        }

        return to_Return;
    }
        break;
    case 14:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Metrics.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Metrics[i].size();j++)
                    to_Return[i].push_back(Population_Metrics[i][j].SD_Phase_Diff);
        }

        return to_Return;
    }
        break;
    case 15:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Eqb_Synch_Cluster_ID.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Eqb_Synch_Cluster_ID[i].size();j++)
                    to_Return[i].push_back((double)Population_Eqb_Synch_Cluster_ID[i][j]);
        }

        return to_Return;
    }
        break;
    case 16:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Eqb_Synch_Amp_Cluster_ID.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Eqb_Synch_Amp_Cluster_ID[i].size();j++)
                    to_Return[i].push_back((double)Population_Eqb_Synch_Amp_Cluster_ID[i][j]);
        }

        return to_Return;
    }
        break;
    case 17:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Time_Series_Average.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Time_Series_Average[i].size();j++)
                    to_Return[i].push_back((double)Population_Time_Series_Average[i][j]);
        }

        return to_Return;
    }
        break;
    case 18:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Time_Series_CV.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Time_Series_CV[i].size();j++)
                    to_Return[i].push_back((double)Population_Time_Series_CV[i][j]);
        }

        return to_Return;
    }
        break;
    case 19:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Time_Series_Skewness.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Time_Series_Skewness[i].size();j++)
                    to_Return[i].push_back((double)Population_Time_Series_Skewness[i][j]);
        }

        return to_Return;
    }
        break;
    case 20:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Time_Series_Minimum.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Time_Series_Minimum[i].size();j++)
                    to_Return[i].push_back((double)Population_Time_Series_Minimum[i][j]);
        }

        return to_Return;
    }
        break;
    case 21:
    {
        std::vector<std::vector<double>> to_Return;

        for(int i=0;i<Population_Time_Series_Average.size();i++)
        {
            to_Return.push_back(std::vector<double>());
            for(int j=0;j<Population_Time_Series_Maximum[i].size();j++)
                    to_Return[i].push_back((double)Population_Time_Series_Maximum[i][j]);
        }

        return to_Return;
    }
        break;

    }
    return std::vector<std::vector<double>>(1,std::vector<double>(1,0));
}

std::vector<std::pair<double,double>> Dynamics::Get_Population_Time_Series(int Species, int Patch) const
{
    return Population_Metrics[Species][Patch].Example_Time_series;
}

int Dynamics::n_Community_Measures() const
{
    return n_Community_Metrics;
}

int Dynamics::n_Species_Measures() const
{
    return n_Species_Metrics;
}

int Dynamics::n_Patch_Measures() const
{
    return n_Patch_Metrics;
}

int Dynamics::n_Population_Measures() const
{
    return n_Population_Metrics;
}
