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
#include "write_simulation_data_functions.h"
//#include <boost/filesystem.hpp>
#include <fstream>

void Write_Timeseries(std::string Write_Path, int nSpecies, const boost_matrix &Time_Series)
{
    std::ofstream Timeseries_Output( Write_Path+"Timeseries.csv");

    Timeseries_Output<<"t,";

    for(unsigned i_col = 0; i_col<Time_Series.shape()[1]/nSpecies;i_col++)
       for(unsigned i_sp_col = 0;i_sp_col<nSpecies;i_sp_col++)
            Timeseries_Output<<"s"<<i_sp_col<<"_p"<<i_col<<",";

    Timeseries_Output<<std::endl;

    for(unsigned i_row = 0; i_row < Time_Series.shape()[0];i_row++)
    {
        for(unsigned i_col = 0; i_col <Time_Series.shape()[1];i_col++)
            Timeseries_Output<<Time_Series[i_row][i_col]<<",";

        Timeseries_Output<<std::endl;
    }

    Timeseries_Output.close();
}

//For New_Simulate

void Start_Time_series(std::string Write_File, int nSpecies, int nPatches)
{
    std::ofstream Timeseries_Output(Write_File);

    Timeseries_Output<<"t,";

    for(unsigned i_col = 0; i_col<nPatches;i_col++)
        for(unsigned i_sp_col = 0;i_sp_col<nSpecies;i_sp_col++)
            Timeseries_Output<<"s"<<i_sp_col<<"_p"<<i_col<<",";

    Timeseries_Output<<std::endl;

    Timeseries_Output.close();
}

void Write_or_Append_Time_series(std::string Write_File, int nSpecies, const boost_matrix &Time_Series)
{
    std::ofstream Timeseries_Output(Write_File,std::ios_base::app);

    for(unsigned i_row = 0; i_row < Time_Series.shape()[0];i_row++)
    {
        for(unsigned i_col = 0; i_col <Time_Series.shape()[1];i_col++)
            Timeseries_Output<<Time_Series[i_row][i_col]<<",";

        Timeseries_Output<<std::endl;
    }

    Timeseries_Output.close();

//    boost::filesystem::path Check_File(Write_File);

//    if(boost::filesystem::exists(Check_File))
//    {
//        std::ofstream Timeseries_Output(Write_File,std::ios_base::app);

//        for(unsigned i_row = 0; i_row < Time_Series.shape()[0];i_row++)
//        {
//            for(unsigned i_col = 0; i_col <Time_Series.shape()[1];i_col++)
//                Timeseries_Output<<Time_Series[i_row][i_col]<<",";

//            Timeseries_Output<<std::endl;
//        }

//        Timeseries_Output.close();
//    }
//    else
//    {
//        std::ofstream Timeseries_Output(Write_File);

//        Timeseries_Output<<"t,";

//        for(unsigned i_col = 0; i_col<Time_Series.shape()[1]/nSpecies;i_col++)
//            for(unsigned i_sp_col = 0;i_sp_col<nSpecies;i_sp_col++)
//                Timeseries_Output<<"s"<<i_sp_col<<"_p"<<i_col<<",";

//        Timeseries_Output<<std::endl;

//        for(unsigned i_row = 0; i_row < Time_Series.shape()[0];i_row++)
//        {
//            for(unsigned i_col = 0; i_col <Time_Series.shape()[1];i_col++)
//                Timeseries_Output<<Time_Series[i_row][i_col]<<",";

//            Timeseries_Output<<std::endl;
//        }

//        Timeseries_Output.close();
//    }
}

void Write_Food_Web_Summary_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output, int n_Food_Webs)
{
    std::ofstream Write_FW_Summary(Write_Path+"Food_Web_Summary.csv");

    Write_FW_Summary<<"iSet,FW_ID,SS,";

    Write_FW_Summary<<"H,AMI,Hc,N_Producers,Min_Trophic_Levels,Min_Consumer_Resource_Size_Ratio,Max_Consumer_Resource_Size_Ratio,Avg_Consumer_Resource_Size_Ratio,";

    Write_FW_Summary<<"Feas_Frequency,";

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels.size();i++)
    {
        Write_FW_Summary<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<",";
        Write_FW_Summary<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<".SD,";
    }

    Write_FW_Summary<<std::endl;

    for(int i_Set=0;i_Set<Output.size();i_Set++)
    {
        std::vector<std::vector<double>> Feas_Freq(2,std::vector<double>(n_Food_Webs,0));
        std::vector<std::vector<double>> n_Sims(2,std::vector<double>(n_Food_Webs,0));

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Sim_Out.size();i_Set_2++)
        {
            if(Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique==-1)
            {
                Feas_Freq[0][Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique] += Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.Metacommunity_Feasible;
                n_Sims[0][Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique]++;
            }
            else
            {
                Feas_Freq[1][Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique] += Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.Metacommunity_Feasible;
                n_Sims[1][Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique]++;
            }
        }

        for(int i=0;i<Feas_Freq.size();i++)
            for(int i2=0;i2<Feas_Freq[i].size();i2++)
            {
                Feas_Freq[i][i2] /= n_Sims[i][i2];
            }

        std::vector<double> H_list(n_Food_Webs,0);
        std::vector<double> AMI_list(n_Food_Webs,0);
        std::vector<double> Hc_list(n_Food_Webs,0);

        std::vector<double> N_Producers_list(n_Food_Webs,0);
        std::vector<double> N_Trophic_Levels_list(n_Food_Webs,0);

        std::vector<double> Min_Consumer_Resource_Size_Ratio(n_Food_Webs,0);
        std::vector<double> Max_Consumer_Resource_Size_Ratio(n_Food_Webs,0);
        std::vector<double> Avg_Consumer_Resource_Size_Ratio(n_Food_Webs,0);

        std::vector<std::vector<double>> Means(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Food_Webs,0));
        std::vector<std::vector<double>> n_Measures(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Food_Webs,0));

        std::vector<std::vector<double>> SDs(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Food_Webs,0));

        std::vector<double> H_list_SS(n_Food_Webs,0);
        std::vector<double> AMI_list_SS(n_Food_Webs,0);
        std::vector<double> Hc_list_SS(n_Food_Webs,0);

        std::vector<double> N_Producers_list_SS(n_Food_Webs,0);
        std::vector<double> N_Trophic_Levels_list_SS(n_Food_Webs,0);

        std::vector<double> Min_Consumer_Resource_Size_Ratio_SS(n_Food_Webs,0);
        std::vector<double> Max_Consumer_Resource_Size_Ratio_SS(n_Food_Webs,0);
        std::vector<double> Avg_Consumer_Resource_Size_Ratio_SS(n_Food_Webs,0);

        std::vector<std::vector<double>> Means_SS(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Food_Webs,0));
        std::vector<std::vector<double>> n_Measures_SS(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Food_Webs,0));

        std::vector<std::vector<double>> SDs_SS(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Food_Webs,0));

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
        {
            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
            {
                if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique==-1)
                {
                    Means[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] += Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i);
                    n_Measures[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique]++;

                    H_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.H;
                    AMI_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.AMI;
                    Hc_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Hc;

                    N_Producers_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = std::accumulate(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.is_Producer.begin(),
                                                                                                                                                    Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.is_Producer.end(),
                                                                                                                                                    (double)0);

                    N_Trophic_Levels_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = *std::max_element(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Min_Trophic_Position.begin(),
                                                                                                                                                         Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Min_Trophic_Position.end());

                    Min_Consumer_Resource_Size_Ratio[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Community_Min_Consumer_Resource_Size_Ratio;
                    Max_Consumer_Resource_Size_Ratio[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Community_Max_Consumer_Resource_Size_Ratio;
                    Avg_Consumer_Resource_Size_Ratio[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Community_Avg_Consumer_Resource_Size_Ratio;
                }
                else
                {
                    Means_SS[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] += Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i);
                    n_Measures_SS[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique]++;

                    H_list_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.H;
                    AMI_list_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.AMI;
                    Hc_list_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Hc;

                    N_Producers_list_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = std::accumulate(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.is_Producer.begin(),
                                                                                                                                                     Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.is_Producer.end(),
                                                                                                                                                     (double)0);

                    N_Trophic_Levels_list_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = *std::max_element(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Min_Trophic_Position.begin(),
                                                                                                                                                            Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Min_Trophic_Position.end());

                    Min_Consumer_Resource_Size_Ratio_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Community_Min_Consumer_Resource_Size_Ratio;
                    Max_Consumer_Resource_Size_Ratio_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Community_Max_Consumer_Resource_Size_Ratio;
                    Avg_Consumer_Resource_Size_Ratio_SS[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Community_Avg_Consumer_Resource_Size_Ratio;
                }
            }
        }

        for(int i=0;i<Means.size();i++)
            for(int i2=0;i2<Means[i].size();i2++)
            {
                Means[i][i2] /= n_Measures[i][i2];
                Means_SS[i][i2] /= n_Measures_SS[i][i2];
            }

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
        {
            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
            {
                if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique==-1)
                    SDs[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] += pow(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i)-Means[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique],2);
                else
                    SDs_SS[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique] += pow(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i)-Means_SS[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique],2);
            }
        }

        for(int i=0;i<SDs.size();i++)
            for(int i2=0;i2<SDs[i].size();i2++)
            {
                SDs[i][i2] = sqrt(SDs[i][i2]);
                SDs_SS[i][i2] = sqrt(SDs_SS[i][i2]);
            }

        for(int i_FW = 0;i_FW<n_Food_Webs;i_FW++)
        {
            if(Feas_Freq[0][i_FW]>0||Feas_Freq[1][i_FW]>0)
            {
                //No SS Row

                if(Output[i_Set].Single_Patch_Simulations.size()>0)
                {
                    Write_FW_Summary<<i_Set<<","<<i_FW<<","<<0<<",";

                    Write_FW_Summary<<H_list[i_FW]<<","<<AMI_list[i_FW]<<","<<Hc_list[i_FW]<<","<<N_Producers_list[i_FW]<<","<<N_Trophic_Levels_list[i_FW]<<","<<Min_Consumer_Resource_Size_Ratio[i_FW]<<","<<Max_Consumer_Resource_Size_Ratio[i_FW]<<","<<Avg_Consumer_Resource_Size_Ratio[i_FW]<<",";

                    Write_FW_Summary<<Feas_Freq[0][i_FW]<<",";

                    for(int i=0;i<Means.size();i++)
                        Write_FW_Summary<<Means[i][i_FW]<<","<<SDs[i][i_FW]<<",";

                    Write_FW_Summary<<std::endl;
                }

                //SS Row
                Write_FW_Summary<<i_Set<<","<<i_FW<<","<<1<<",";

                Write_FW_Summary<<H_list_SS[i_FW]<<","<<AMI_list_SS[i_FW]<<","<<Hc_list_SS[i_FW]<<","<<N_Producers_list_SS[i_FW]<<","<<N_Trophic_Levels_list_SS[i_FW]<<","<<Min_Consumer_Resource_Size_Ratio_SS[i_FW]<<","<<Max_Consumer_Resource_Size_Ratio_SS[i_FW]<<","<<Avg_Consumer_Resource_Size_Ratio_SS[i_FW]<<",";

                Write_FW_Summary<<Feas_Freq[1][i_FW]<<",";

                for(int i=0;i<Means_SS.size();i++)
                    Write_FW_Summary<<Means_SS[i][i_FW]<<","<<SDs_SS[i][i_FW]<<",";

                Write_FW_Summary<<std::endl;
            }
        }
    }
}

void Write_Spatial_Structure_Summary_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output, int n_Spatial_Structures)
{
    std::ofstream Write_SS_Summary(Write_Path+"Spatial_Structure_Summary.csv");

    Write_SS_Summary<<"iSet,SS_ID,";

    Write_SS_Summary<<"Feas_Frequency,";

    Write_SS_Summary<<"Average_Degree,Degree_Skewness,Average_Path_Length,Clustering_Coefficient,Eigenratio,";

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels.size();i++)
    {
        Write_SS_Summary<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<",";
        Write_SS_Summary<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<".SD,";
    }

    Write_SS_Summary<<std::endl;

    for(int i_Set=0;i_Set<Output.size();i_Set++)
    {
        std::vector<double> Feas_Freq(n_Spatial_Structures+1,0);
        std::vector<double> n_Sims(n_Spatial_Structures+1,0);

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Sim_Out.size();i_Set_2++)
        {
            if(Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique>=0)
            {
                Feas_Freq[Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique] += Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.Metacommunity_Feasible;
                n_Sims[Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique]++;
            }
            else if(Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique==-1)
            {
                Feas_Freq[n_Spatial_Structures] += Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.Metacommunity_Feasible;
                n_Sims[n_Spatial_Structures]++;
            }

        }

        for(int i=0;i<Feas_Freq.size();i++)
                Feas_Freq[i] /= n_Sims[i];

        std::vector<double> Average_Degree_list(n_Spatial_Structures+1,0);
        std::vector<double> Degree_Skewness_list(n_Spatial_Structures+1,0);
        std::vector<double> Average_Path_Length_list(n_Spatial_Structures+1,0);
        std::vector<double> Clustering_Coefficient_list(n_Spatial_Structures+1,0);
        std::vector<double> Eigenratio_list(n_Spatial_Structures+1,0);

        std::vector<std::vector<double>> Means(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Spatial_Structures+1,0));
        std::vector<std::vector<double>> n_Measures(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Spatial_Structures+1,0));

        std::vector<std::vector<double>> SDs(Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),std::vector<double>(n_Spatial_Structures+1,0));

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
        {
            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
            {
                if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique>=0)
                {
                    Means[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] += Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i);
                    n_Measures[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique]++;

                    Average_Degree_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Average_Degree;
                    Degree_Skewness_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Degree_Skewness;
                    Average_Path_Length_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Average_Path_Length;
                    Clustering_Coefficient_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Clustering_Coefficient;
                    Eigenratio_list[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Eigenratio;
                }
                else if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique==-1)
                {
                    Means[i][n_Spatial_Structures] += Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i);
                    n_Measures[i][n_Spatial_Structures]++;
                }
            }
        }

        for(int i=0;i<Means.size();i++)
            for(int i2=0;i2<Means[i].size();i2++)
                Means[i][i2] /= n_Measures[i][i2];

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
                if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique>=0)
                    SDs[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique] += pow(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i)-Means[i][Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique],2);

        for(int i=0;i<SDs.size();i++)
            for(int i2=0;i2<SDs[i].size();i2++)
                SDs[i][i2] = sqrt(SDs[i][i2]);

        if(n_Sims[n_Spatial_Structures]>0)
        {
        //Write for isolated metas
            Write_SS_Summary<<i_Set<<","<<-1<<",";

            Write_SS_Summary<<Feas_Freq[n_Spatial_Structures]<<",";

            Write_SS_Summary<<Average_Degree_list[n_Spatial_Structures]<<","<<Degree_Skewness_list[n_Spatial_Structures]<<","<<Average_Path_Length_list[n_Spatial_Structures]<<","<<Clustering_Coefficient_list[n_Spatial_Structures]<<","<<Eigenratio_list[n_Spatial_Structures]<<",";

            for(int i=0;i<Means.size();i++)
                Write_SS_Summary<<Means[i][n_Spatial_Structures]<<","<<SDs[i][n_Spatial_Structures]<<",";

            Write_SS_Summary<<std::endl;
        }

        for(int i_SS = 0;i_SS<n_Spatial_Structures;i_SS++)
        {
            Write_SS_Summary<<i_Set<<","<<i_SS<<",";

            Write_SS_Summary<<Feas_Freq[i_SS]<<",";

            Write_SS_Summary<<Average_Degree_list[i_SS]<<","<<Degree_Skewness_list[i_SS]<<","<<Average_Path_Length_list[i_SS]<<","<<Clustering_Coefficient_list[i_SS]<<","<<Eigenratio_list[i_SS]<<",";

            for(int i=0;i<Means.size();i++)
                Write_SS_Summary<<Means[i][i_SS]<<","<<SDs[i][i_SS]<<",";

            Write_SS_Summary<<std::endl;
        }
    }
}

void Write_Population_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output)
{
    std::ofstream Write_Population(Write_Path+"Population_Data.csv");

    Write_Population<<"Index,iSet,unique_FW_ID,unique_SS_ID,Species,Patch,";

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Population_Dynamics_Labels.size();i++)
        Write_Population<<Output[0].Sim_Out[0].Dyn_Out.Population_Dynamics_Labels[i]<<",";

    Write_Population<<std::endl;

        for(int i_Set=0;i_Set<Output.size();i_Set++)
            for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
            {
                int nSpecies = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Food_Web.shape()[0];
                int nPatches = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0];

                for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                {
                    for(int i_Pa = 0;i_Pa<nPatches;i_Pa++)
                    {
                        Write_Population<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID<<","
                                                                                                           <<i_Set<<","
                                                                                                          <<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique<<","
                                                                                                         <<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique<<","
                                                                                                         <<i_Sp<<","
                                                                                                        <<i_Pa<<",";

                        for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Population_Measures();i++)
                            Write_Population<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Population_Measures(i)[i_Sp][i_Pa]<<",";

                        Write_Population<<std::endl;
                    }
                }
            }

    Write_Population.close();
}

void Write_Population_Cycle_Timeseries(std::string Write_Path, double Simulation_Resolution, const Dynamics &Output)
{
    int nSpecies = Output.Species_Avg_Eqb.size(); // Not best access solution
    int nPatches = Output.Community_Avg_Eqb.size(); // Not best access solution

    std::ofstream Write_Limit_Cycle(Write_Path+"_Limit_Cycle.csv");

    std::vector<std::vector<double>> Time_Series;

    Write_Limit_Cycle<<"Times,";

    for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
        for(int i_Pa=0;i_Pa<nPatches;i_Pa++)
        {
            Time_Series.push_back(std::vector<double>());

            std::vector<std::pair<double,double>> Pop_Time_Series = Output.Get_Population_Time_Series(i_Sp,i_Pa);

            for(int i_Ts = 0;i_Ts<Pop_Time_Series.size();i_Ts++)
                Time_Series.back().push_back(Pop_Time_Series[i_Ts].second);

            Write_Limit_Cycle<<"s"<<i_Sp<<"_p"<<i_Pa<<",";

        }

    Write_Limit_Cycle<<std::endl;

    int i_row = 0;
    bool rows_to_write = true;
    while(rows_to_write)
    {
        rows_to_write = false;
        std::vector<double> Temp_Row;

        for(int i=0;i<Time_Series.size();i++)
        {
            if(Time_Series[i].size()>i_row)
            {
                Temp_Row.push_back(Time_Series[i][i_row]);
                rows_to_write = true;
            }
            else
            {
                Temp_Row.push_back(-1);
            }
        }

        if(rows_to_write)
        {
            Write_Limit_Cycle<<((double)i_row)*Simulation_Resolution<<",";

            for(int i=0;i<Temp_Row.size();i++)
            {
                Write_Limit_Cycle<<Temp_Row[i]<<",";
            }

            Write_Limit_Cycle<<std::endl;
        }

        i_row++;
    }
    Write_Limit_Cycle.close();
}

void Write_Population_Cycle_Timeseries(std::string Write_Path, double Simulation_Resolution, const std::vector<Sim_Set_Output> &Output)
{
    for(int i_Set=0;i_Set<Output.size();i_Set++)
        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
        {
            std::ofstream Write_Timeseries(Write_Path+"_"+std::to_string(static_cast<long long>(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID))+"_"+"Population_Cycle_Timeseries.csv");

            Write_Timeseries<<"Times,";

            //Write_Timeseries<<"Species,Patch,Type"<<std::endl;

            int nSpecies = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Food_Web.shape()[0];
            int nPatches = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0];

            std::vector<std::vector<double>> Time_Series;

            for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                for(int i_Pa=0;i_Pa<nPatches;i_Pa++)
                {
                    Time_Series.push_back(std::vector<double>());

                    std::vector<std::pair<double,double>> Pop_Time_Series = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Population_Time_Series(i_Sp,i_Pa);

                    for(int i_Ts = 0;i_Ts<Pop_Time_Series.size();i_Ts++)
                        Time_Series.back().push_back(Pop_Time_Series[i_Ts].second);

                    Write_Timeseries<<"s"<<i_Sp<<"_p"<<i_Pa<<",";

                }

            Write_Timeseries<<std::endl;

            int i_row = 0;
            bool rows_to_write = true;
            while(rows_to_write)
            {
                rows_to_write = false;
                std::vector<double> Temp_Row;

                for(int i=0;i<Time_Series.size();i++)
                {
                    if(Time_Series[i].size()>i_row)
                    {
                        Temp_Row.push_back(Time_Series[i][i_row]);
                        rows_to_write = true;
                    }
                    else
                    {
                        Temp_Row.push_back(-1);
                    }
                }

                if(rows_to_write)
                {
                    Write_Timeseries<<((double)i_row)*Simulation_Resolution<<",";

                    for(int i=0;i<Temp_Row.size();i++)
                    {
                        Write_Timeseries<<Temp_Row[i]<<",";
                    }

                    Write_Timeseries<<std::endl;
                }

                i_row++;
            }
            Write_Timeseries.close();
        }
}

void Write_Species_Output(std::string Write_Path, bool Feasible_Only, const std::vector<Sim_Set_Output> &Output)
{
    if(Feasible_Only)
        Write_Path = Write_Path + "Feas_";
    else
        Write_Path = Write_Path + "NonFeas_";

    std::ofstream Write_Species(Write_Path+"Species_Data.csv");

    Write_Species<<"Index,iSet,unique_FW_ID,unique_SS_ID,Species,";

    for(int i=0;i<Output[0].Sim_Out[0].Eqn_In.Species_Data_Labels.size();i++)
        Write_Species<<Output[0].Sim_Out[0].Eqn_In.Species_Data_Labels[i]<<",";

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Species_Dynamics_Labels.size();i++)
        Write_Species<<Output[0].Sim_Out[0].Dyn_Out.Species_Dynamics_Labels[i]<<",";

    Write_Species<<std::endl;

    if(Feasible_Only)
        {
            for(int i_Set=0;i_Set<Output.size();i_Set++)
                for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
                {
                    int nSpecies = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Food_Web.shape()[0];

                    for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                    {
                        Write_Species<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID<<","
                                                                                                           <<i_Set<<","
                                                                                                          <<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique<<","
                                                                                                         <<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique<<","
                                                                                                         <<i_Sp<<",";

                        for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Species_Data();i++)
                            Write_Species<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Get_Species_Data(i)[i_Sp]<<",";

                        for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Species_Measures();i++)
                            Write_Species<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Species_Measures(i)[i_Sp]<<",";

                        Write_Species<<std::endl;

                    }
                }
        }
    else
    {
        for(int i_Set=0;i_Set<Output.size();i_Set++)
            for(int i_Set_2=0;i_Set_2<Output[i_Set].Non_Feasible_Simulations.size();i_Set_2++)
            {
                int nSpecies = Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Food_Web.shape()[0];

                for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                {
                    Write_Species<<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].ID<<","
                                                                                                       <<i_Set<<","
                                                                                                      <<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique<<","
                                                                                                     <<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique<<","
                                                                                                     <<i_Sp<<",";

                    for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Eqn_In.n_Species_Data();i++)
                        Write_Species<<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Eqn_In.Get_Species_Data(i)[i_Sp]<<",";

                    for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Dyn_Out.n_Species_Measures();i++)
                        Write_Species<<Output[i_Set].Sim_Out[Output[i_Set].Non_Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Species_Measures(i)[i_Sp]<<",";

                    Write_Species<<std::endl;

                }
            }
    }

    Write_Species.close();

}

void Write_Community_Output(std::string Write_Path, bool Feasible_Only, const std::vector<Sim_Set_Output> &Output)
{    
    if(Feasible_Only)
        Write_Path = Write_Path + "Feas_";

    std::ofstream Write_Summary(Write_Path+"Community_Data.csv");

    //Column Names

    Write_Summary<<"Index,iSet,";

//    for(int i_Gen = 0;i_Gen<Output[0].Set_Gen_Parameter_Labels.size();i_Gen++)
//        for(int i = 0; i<Output[0].Set_Gen_Parameter_Labels[i_Gen].size();i++)
//            Write_Summary<<Output[0].Set_Gen_Parameter_Labels[i_Gen][i]<<",";

    for(int i=0;i<Output[0].Sim_Out[0].Eqn_In.Community_Data_Labels.size();i++)
        Write_Summary<<Output[0].Sim_Out[0].Eqn_In.Community_Data_Labels[i]<<",";

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels.size();i++)
        Write_Summary<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<",";

    Write_Summary<<std::endl;

    //Write Data

    if(Feasible_Only)
    {
        for(int i_Set=0;i_Set<Output.size();i_Set++)
            for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
            {
                //Sim_Output Output_Temp(Output[i_Set].Sim_Out[i_Set_2]);

                Write_Summary<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID<<","<<i_Set<<",";

//                for(int i_Gen = 0;i_Gen<Output[i_Set].Set_Gen_Parameter_Values.size();i_Gen++)
//                    for(int i = 0; i<Output[i_Set].Set_Gen_Parameter_Values[i_Gen].size();i++)
//                        Write_Summary<<Output[i_Set].Set_Gen_Parameter_Values[i_Gen][i]<<",";

                for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Community_Data();i++)
                    Write_Summary<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Get_Community_Data(i)<<",";

                for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
                    Write_Summary<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i)<<",";

                Write_Summary<<std::endl;
            }
    }
    else
    {
        for(int i_Set=0;i_Set<Output.size();i_Set++)
            for(int i_Set_2=0;i_Set_2<Output[i_Set].Sim_Out.size();i_Set_2++)
            {
                //Sim_Output Output_Temp(Output[i_Set].Sim_Out[i_Set_2]);

                Write_Summary<<Output[i_Set].Sim_Out[i_Set_2].ID<<","<<i_Set<<",";

//                for(int i_Gen = 0;i_Gen<Output[i_Set].Set_Gen_Parameter_Values.size();i_Gen++)
//                    for(int i = 0; i<Output[i_Set].Set_Gen_Parameter_Values[i_Gen].size();i++)
//                        Write_Summary<<Output[i_Set].Set_Gen_Parameter_Values[i_Gen][i]<<",";

                for(int i=0;i<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.n_Community_Data();i++)
                    Write_Summary<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Get_Community_Data(i)<<",";

                for(int i=0;i<Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.n_Community_Measures();i++)
                    Write_Summary<<Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.Get_Community_Measures(i)<<",";

                Write_Summary<<std::endl;
            }
    }

    Write_Summary.close();
}

void Write_Community_Parameters(std::string Write_Path, bool Feasible_Only, const std::vector<Sim_Set_Output> &Output, int n_Food_Webs, int n_Spatial_Structures)
{
    if(Feasible_Only)
        Write_Path = Write_Path + "Feas_";

    std::ofstream Write_Food_Web(Write_Path+"Food_Web_Matricies.csv");
    std::ofstream Write_Spatial_Structure(Write_Path+"Spatial_Structure_Matricies.csv");
    std::ofstream Write_Carrying_Capacity(Write_Path+"Carrying_Capacities.csv");
    std::ofstream Write_Initial_Abundances(Write_Path+"Initial_Abundances.csv");

    int max_nSpecies = 0;
    int max_nPatches = 0;

    for(int i_Set=0;i_Set<Output.size();i_Set++)
        for(int i_Sim=0;i_Sim<Output[i_Set].Sim_Out.size();i_Sim++)
        {
            if(Output[i_Set].Sim_Out[i_Sim].Eqn_In.Food_Web_Data.Food_Web.shape()[0]>max_nSpecies)
            {
                max_nSpecies = Output[i_Set].Sim_Out[i_Sim].Eqn_In.Food_Web_Data.Food_Web.shape()[0];
            }

            if(Output[i_Set].Sim_Out[i_Sim].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0]>max_nPatches)
            {
                max_nPatches = Output[i_Set].Sim_Out[i_Sim].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0];
            }
        }

    Write_Food_Web<<"iSet,unique_ID,";

    for(int i_FW = 0; i_FW<max_nSpecies;i_FW++)
    {
        Write_Food_Web<<"S_"<<i_FW<<",";
    }

    Write_Food_Web<<std::endl;

    Write_Spatial_Structure<<"iSet,unique_ID,";
    Write_Carrying_Capacity<<"i,";
    Write_Initial_Abundances<<"i,";

    for(int i_SS = 0; i_SS<max_nPatches;i_SS++)
    {
        Write_Spatial_Structure<<"P_"<<i_SS<<",";
        Write_Carrying_Capacity<<"P_"<<i_SS<<",";
        Write_Initial_Abundances<<"P_"<<i_SS<<",";
    }

    Write_Spatial_Structure<<std::endl;
    Write_Carrying_Capacity<<std::endl;
    Write_Initial_Abundances<<std::endl;

    if(Feasible_Only)
    {
        int count_i = 0;

        for(int i_Set=0;i_Set<Output.size();i_Set++)
        {
            std::vector<bool> Food_Web_IDs(n_Food_Webs,true);
            std::vector<bool> Spatial_Structure_IDs(n_Spatial_Structures,true);

            for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
            {
                int nSpecies = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Food_Web.shape()[0];
                int nPatches = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0];

                if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique>=0)
                    if(Food_Web_IDs[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique])
                    {
                        for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                        {

                            Write_Food_Web<<i_Set<<","<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique<<",";

                            for(int i_Sp_2 = 0;i_Sp_2<nSpecies;i_Sp_2++)
                            {
                                Write_Food_Web<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Data.Food_Web[i_Sp][i_Sp_2]<<",";
                            }

                            Write_Food_Web<<std::endl;


                        }

                        Food_Web_IDs[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique]=false;
                    }

                for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                {
                    Write_Carrying_Capacity<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID<<",";
                    Write_Initial_Abundances<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID<<",";

                    for(int i_Pa = 0;i_Pa<nPatches;i_Pa++)
                    {
                        Write_Carrying_Capacity<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.k[i_Sp][i_Pa]<<",";
                        Write_Initial_Abundances<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Initial_Abundance[i_Sp][i_Pa]<<",";
                    }

                    Write_Carrying_Capacity<<std::endl;
                    Write_Initial_Abundances<<std::endl;
                }

                if(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique>=0)
                    if(Spatial_Structure_IDs[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique])
                    {
                        for(int i_Pa = 0;i_Pa<nPatches;i_Pa++)
                        {
                            Write_Spatial_Structure<<i_Set<<","<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique<<",";

                            for(int i_Pa_2 = 0;i_Pa_2<nPatches;i_Pa_2++)
                            {
                                Write_Spatial_Structure<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Spatial_Structure[i_Pa][i_Pa_2]<<",";
                            }

                            Write_Spatial_Structure<<std::endl;

                        }

                        Spatial_Structure_IDs[Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique]=false;
                    }
                    count_i++;
            }
        }
    }
    else
    {
        int count_i = 0;

        for(int i_Set=0;i_Set<Output.size();i_Set++)
        {
            std::vector<bool> Food_Web_IDs(n_Food_Webs,true);
            std::vector<bool> Spatial_Structure_IDs(n_Spatial_Structures,true);

            for(int i_Set_2=0;i_Set_2<Output[i_Set].Sim_Out.size();i_Set_2++)
            {
                int nSpecies = Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Data.Food_Web.shape()[0];
                int nPatches = Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0];

                if(Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique>=0)
                    if(Food_Web_IDs[Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique])
                    {

                        for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                        {


                            Write_Food_Web<<i_Set<<","<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique<<",";


                            for(int i_Sp_2 = 0;i_Sp_2<nSpecies;i_Sp_2++)
                            {
                                Write_Food_Web<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Data.Food_Web[i_Sp][i_Sp_2]<<",";
                            }

                            Write_Food_Web<<std::endl;


                        }

                        Food_Web_IDs[Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Food_Web_Unique]=false;

                    }

                for(int i_Sp = 0;i_Sp<nSpecies;i_Sp++)
                {                    
                    Write_Carrying_Capacity<<Output[i_Set].Sim_Out[i_Set_2].ID<<",";
                    Write_Initial_Abundances<<Output[i_Set].Sim_Out[i_Set_2].ID<<",";

                    for(int i_Pa = 0;i_Pa<nPatches;i_Pa++)
                    {
                        Write_Carrying_Capacity<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.k[i_Sp][i_Pa]<<",";
                        Write_Initial_Abundances<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Initial_Abundance[i_Sp][i_Pa]<<",";
                    }

                    Write_Carrying_Capacity<<std::endl;
                    Write_Initial_Abundances<<std::endl;

                }

                if(Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique>=0)
                    if(Spatial_Structure_IDs[Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique])
                    {

                        for(int i_Pa = 0;i_Pa<nPatches;i_Pa++)
                        {
                            Write_Spatial_Structure<<i_Set<<","<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique<<",";

                            for(int i_Pa_2 = 0;i_Pa_2<nPatches;i_Pa_2++)
                            {
                                Write_Spatial_Structure<<Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Data.Spatial_Structure[i_Pa][i_Pa_2]<<",";
                            }

                            Write_Spatial_Structure<<std::endl;

                        }
                        Spatial_Structure_IDs[Output[i_Set].Sim_Out[i_Set_2].Eqn_In.Spatial_Structure_Unique]=false;

                    }
                    count_i++;
            }
        }
    }

    Write_Food_Web.close();
    Write_Spatial_Structure.close();
    Write_Carrying_Capacity.close();
    Write_Initial_Abundances.close();
}

void Write_Set_Summary_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output)
{
    std::ofstream Write_Set(Write_Path+"Set_Summary.csv");

    //Column Names

    Write_Set<<"iSet,Feas_Freq,";

        for(int i_Gen = 0;i_Gen<Output[0].Set_Gen_Parameter_Labels.size();i_Gen++)
            for(int i = 0; i<Output[0].Set_Gen_Parameter_Labels[i_Gen].size();i++)
                Write_Set<<Output[0].Set_Gen_Parameter_Labels[i_Gen][i]<<",";

    for(int i=0;i<Output[0].Sim_Out[0].Eqn_In.Community_Data_Labels.size();i++)
    {
        Write_Set<<Output[0].Sim_Out[0].Eqn_In.Community_Data_Labels[i]<<",";
        Write_Set<<Output[0].Sim_Out[0].Eqn_In.Community_Data_Labels[i]<<".SD,";
    }

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels.size();i++)
    {
        Write_Set<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<",";
        Write_Set<<Output[0].Sim_Out[0].Dyn_Out.Community_Dynamics_Labels[i]<<".SD,";
    }

    Write_Set<<std::endl;

    //Write Data

    for(int i_Set=0;i_Set<Output.size();i_Set++)
    {
        double Feasible_Frequency = 0;

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Sim_Out.size();i_Set_2++)
        {
            Feasible_Frequency+=Output[i_Set].Sim_Out[i_Set_2].Dyn_Out.Metacommunity_Feasible;
        }

        Feasible_Frequency /= Output[i_Set].Sim_Out.size();

        std::vector<double> Means(Output[i_Set].Sim_Out[0].Eqn_In.n_Community_Data()+Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),0);
        std::vector<double> SDs(Output[i_Set].Sim_Out[0].Eqn_In.n_Community_Data()+Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures(),0);

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
        {
            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Community_Data();i++)
                Means[i] += Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Get_Community_Data(i);

            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
                Means[i+Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Community_Data()] += Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i);
        }

        for(int i=0;i<Output[i_Set].Sim_Out[0].Eqn_In.n_Community_Data();i++)
            Means[i] /= Output[i_Set].Feasible_Simulations.size();

        for(int i=0;i<Output[i_Set].Sim_Out[0].Dyn_Out.n_Community_Measures();i++)
            Means[i+Output[i_Set].Sim_Out[0].Eqn_In.n_Community_Data()] /= Output[i_Set].Feasible_Simulations.size();

        for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
        {
            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Community_Data();i++)
                SDs[i] += pow(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Get_Community_Data(i)-Means[i],2);

            for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Community_Measures();i++)
                SDs[i+Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Community_Data()] += pow(Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Community_Measures(i)-Means[i+Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Community_Data()],2);
        }

        for(int i=0;i<SDs.size();i++)
            SDs[i] = sqrt(SDs[i]);

        Write_Set<<i_Set<<",";

        Write_Set<<Feasible_Frequency<<",";

        for(int i_Gen = 0;i_Gen<Output[i_Set].Set_Gen_Parameter_Values.size();i_Gen++)
            for(int i = 0; i<Output[i_Set].Set_Gen_Parameter_Values[i_Gen].size();i++)
                Write_Set<<Output[i_Set].Set_Gen_Parameter_Values[i_Gen][i]<<",";

        for(int i=0;i<Means.size();i++)
            Write_Set<<Means[i]<<","<<SDs[i]<<",";

        Write_Set<<std::endl;
    }
}

void Write_Patch_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output)
{
    std::ofstream Write_Patches(Write_Path+"Patch_Data.csv");

    Write_Patches<<"Index,iSet,unique_FW_ID,unique_SS_ID,Species,";

    for(int i=0;i<Output[0].Sim_Out[0].Eqn_In.Patches_Data_Labels.size();i++)
        Write_Patches<<Output[0].Sim_Out[0].Eqn_In.Patches_Data_Labels[i]<<",";

    for(int i=0;i<Output[0].Sim_Out[0].Dyn_Out.Patch_Dynamics_Labels.size();i++)
        Write_Patches<<Output[0].Sim_Out[0].Dyn_Out.Patch_Dynamics_Labels[i]<<",";

    Write_Patches<<std::endl;

            for(int i_Set=0;i_Set<Output.size();i_Set++)
                for(int i_Set_2=0;i_Set_2<Output[i_Set].Feasible_Simulations.size();i_Set_2++)
                {
                    int nPatches = Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Data.Spatial_Structure.shape()[0];

                    for(int i_Pa = 0;i_Pa<nPatches;i_Pa++)
                    {
                        Write_Patches<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].ID<<","
                                                                                                           <<i_Set<<","
                                                                                                          <<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Food_Web_Unique<<","
                                                                                                         <<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Spatial_Structure_Unique<<","
                                                                                                         <<i_Pa<<",";

                        for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.n_Patches_Data();i++)
                            Write_Patches<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Eqn_In.Get_Patches_Data(i)[i_Pa]<<",";

                        for(int i=0;i<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.n_Patch_Measures();i++)
                            Write_Patches<<Output[i_Set].Sim_Out[Output[i_Set].Feasible_Simulations[i_Set_2]].Dyn_Out.Get_Patch_Measures(i)[i_Pa]<<",";

                        Write_Patches<<std::endl;

                    }
                }
}
