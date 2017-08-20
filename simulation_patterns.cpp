//#include "simulation_patterns.h"

//void Simulate_Random(Parameter_Generators x_In, double Simulation_Length, double Simulation_Resolution, std::string Write_Path, int Time_Series_Print, double Extinction_Threshold)
//{

//    //Generate all relevant parameters for simulation from x_In (requires mutex)
//    //omp_set_lock(&gen_mutex);

//    std::vector<std::vector<double>> Generator_Method_Data;

//    for(int i=0;i<x_In.nGenerators();i++)
//    {
//        x_In[i]->Generate_New();
//        Generator_Method_Data.push_back(x_In[i]->Get_Generator_Method_Data());
//    }

//    Generator_Method_Data.push_back(x_In.Get_Generator_Set_Method_Data());

//    //Set parameters for use in the simulation function

//    //Eqn_Parameters Parameter_Set(x_In);

//    Pass_to_Fortran::Thread_Local_Parameters_ptr.reset(new Eqn_Parameters(x_In));

//    //Pass_to_Fortran::Thread_Local_Parameters_ptr.reset(&Parameter_Set);

//    boost_matrix Input_Initial_Abundances(x_In.Initial_Abundance->Get_Value());

//    //Put initial abundances into format for simulation function
//    std::vector<double> Initial_Abundances;

//    for(int iP = 0; iP<(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches;iP++)
//        for(int iS = 0;iS<(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies;iS++)
//            Initial_Abundances.push_back(Input_Initial_Abundances[iS][iP]);

//    //Call solver
//    {
//        boost_matrix sim_ana = LSODA_Solve(Initial_Abundances,
//                                           &Pass_to_Fortran::Allo_Web_Eqn,
//                                           &Pass_to_Fortran::JDUM,
//                                           Simulation_Length,
//                                           Simulation_Resolution,
//                                           Extinction_Threshold,
//                                           (*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies);

//        //Save results

//        if(sim_ana.size()==vSequence(1.0,Simulation_Length,Simulation_Resolution).size())
//        {
//            Dynamics Measured = Measure_Metacommunity_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches, sim_ana);

//            omp_set_lock(&write_mutex);

//            iSet_Index.push_back(iSet);

//            All_Community_Dynamics[iSet].Add_Dynamics(Measured);
//            Feasible_Community_Dynamics[iSet].Add_Dynamics(Measured);

//            Food_Web_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).Food_Web);
//            Spatial_Structure_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).a_m);

//            Carrying_Capacity_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).k);
//            Initial_Abundances_list.push_back(Input_Initial_Abundances);

//            //Species Level Data

//            Trophic_Position.push_back(x_In.Food_Web->Get_Metadata().Trophic_Position);
//            Generality.push_back(x_In.Food_Web->Get_Metadata().Generality);

//            Mass_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).M);

//            ar_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).ar);
//            ax_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).ax);
//            y_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).y);
//            e_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).e);
//            B0_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).Half_Saturation);
//            h_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).Functional_Response_Shape);
//            Dispersal_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).m);

//            std::vector<int> Dummy_dynamics ((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,0);

//            Extinct.push_back(Dummy_dynamics);
//            Transient_Length.push_back(Measured.Transient_Length_Species);
//            Eqb_Minima.push_back(Measured.Eqb_Max_Minima_Species);
//            Transient_Minima.push_back(Measured.Transient_Max_Minima_Species);
//            Eqb_N_Clusters.push_back(Measured.Eqb_N_Clusters_Species);
//            Amplitude.push_back(Measured.Max_Amplitude_Species);

//            for(unsigned i=0;i<Generator_Cache.size();i++)
//                for(unsigned i2=0;i2<Generator_Method_Data[i].size();i2++)
//                {
//                    Generator_Cache[i].Set_Data[set_counter][i2] = Generator_Method_Data[i][i2];
//                    Generator_Cache[i].Generation_Data[counter][i2] = Generator_Method_Data[i][i2];
//                }

//            if(Time_Series_Print>0)
//            {
//                std::ofstream Timeseries_Output (Write_Path+"Feas_Timeseries_"+std::to_string(static_cast<long long>(counter))+".csv");

//                //Timeseries_Output<<"Trans_Length"<<","<<Trans_Length_list[counter]<<","<<"Eqb_Clusters"<<","<<Clusters_list[counter]<<","<<"Fixed_Point"<<","<<Fixed_Point_list[counter]<<","<<"Amplitude"<<","<<Amplitude_list[counter]<<endl;

//                for(unsigned i_col = 0; i_col<sim_ana.shape()[1];i_col++)
//                {
//                    Timeseries_Output<<i_col<<",";
//                }

//                Timeseries_Output<<std::endl;

//                for(unsigned i_row = 0; i_row < sim_ana.shape()[0];i_row++)
//                {
//                    for(unsigned i_col = 0; i_col <sim_ana.shape()[1];i_col++)
//                    {
//                        Timeseries_Output<<sim_ana[i_row][i_col]<<",";
//                    }

//                    Timeseries_Output<<std::endl;
//                }

//                Timeseries_Output.close();
//            }

//            counter++;
//            set_counter++;

//            omp_unset_lock(&write_mutex);

//        }
//        else
//        {
//            Dynamics Measured;

//            Measured.Empty_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches);

//            omp_set_lock(&write_mutex);

//            iSet_Index.push_back(iSet);

//            All_Community_Dynamics[iSet].Add_Dynamics(Measured);

//            Food_Web_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).Food_Web);
//            Spatial_Structure_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).a_m);

//            Carrying_Capacity_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).k);
//            Initial_Abundances_list.push_back(Input_Initial_Abundances);

//            //Species Level Data

//            //Species_Level_Data Add_to_Cache;

//            Trophic_Position.push_back(x_In.Food_Web->Get_Metadata().Trophic_Position);
//            Generality.push_back(x_In.Food_Web->Get_Metadata().Generality);

//            Mass_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).M);

//            ar_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).ar);
//            ax_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).ax);
//            y_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).y);
//            e_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).e);
//            B0_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).Half_Saturation);
//            h_list.push_back((*Pass_to_Fortran::Thread_Local_Parameters_ptr).Functional_Response_Shape);
//            Dispersal_list.push_back( (*Pass_to_Fortran::Thread_Local_Parameters_ptr).m);

//            std::vector<int> Dummy_dynamics ((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,0);
//            std::vector<double> Dummy_dynamics_too ((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,0);

//            Extinct.push_back(Which_Species_Extinct((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,sim_ana,Extinction_Threshold));
//            Transient_Length.push_back(Dummy_dynamics_too);
//            Eqb_Minima.push_back(Dummy_dynamics_too);
//            Transient_Minima.push_back(Dummy_dynamics_too);
//            Eqb_N_Clusters.push_back(Dummy_dynamics);
//            Amplitude.push_back(Dummy_dynamics_too);

//            for(unsigned i=0;i<Generator_Cache.size();i++)
//                for(unsigned i2=0;i2<Generator_Method_Data[i].size();i2++)
//                {
//                    Generator_Cache[i].Set_Data[set_counter][i2] = Generator_Method_Data[i][i2];
//                    Generator_Cache[i].Generation_Data[counter][i2] = Generator_Method_Data[i][i2];
//                }

//            if(Time_Series_Print>1)
//            {
//                std::ofstream Timeseries_Output (Write_Path+"NonFeas_Timeseries_"+std::to_string(static_cast<long long>(counter))+".csv");

//                //Timeseries_Output<<"Trans_Length"<<","<<Trans_Length_list[counter]<<","<<"Eqb_Clusters"<<","<<Clusters_list[counter]<<","<<"Fixed_Point"<<","<<Fixed_Point_list[counter]<<","<<"Amplitude"<<","<<Amplitude_list[counter]<<endl;

//                for(unsigned i_col = 0; i_col<sim_ana.shape()[1];i_col++)
//                {
//                    Timeseries_Output<<i_col<<",";
//                }

//                Timeseries_Output<<std::endl;

//                for(unsigned i_row = 0; i_row < sim_ana.shape()[0];i_row++)
//                {
//                    for(unsigned i_col = 0; i_col <sim_ana.shape()[1];i_col++)
//                    {
//                        Timeseries_Output<<sim_ana[i_row][i_col]<<",";
//                    }

//                    Timeseries_Output<<std::endl;
//                }

//                Timeseries_Output.close();
//            }

//            counter++;
//            set_counter++;

//            omp_unset_lock(&write_mutex);

//        }
//    }

//    //Generate for isolated community (patch 0 if dissimilar)

//    Pass_to_Fortran::Thread_Local_Parameters_ptr->Set_to_Single_Community();

//    //Put initial abundances into format for simulation function
//    Initial_Abundances.erase(Initial_Abundances.begin(),Initial_Abundances.end());

//    for(int iS = 0;iS<(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies;iS++)
//        Initial_Abundances.push_back(Input_Initial_Abundances[iS][0]);

//    //Call solver

//    {
//        boost_matrix sim_ana_iso = LSODA_Solve(Initial_Abundances,
//                                               &Pass_to_Fortran::Allo_Web_Eqn,
//                                               &Pass_to_Fortran::JDUM,
//                                               Simulation_Length,
//                                               Simulation_Resolution,
//                                               Extinction_Threshold,
//                                               (*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies);

//        //Save results

//        if(sim_ana_iso.size()==vSequence(1.0,Simulation_Length,Simulation_Resolution).size())
//        {
//            Dynamics Measured = Measure_Metacommunity_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches, sim_ana_iso);

//            omp_set_lock(&write_mutex);

//            Individual_Community_Dynamics[iSet].Add_Dynamics(Measured);

//            //        std::vector<int> Dummy_dynamics ((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,0);

//            //        Extinct.push_back(Dummy_dynamics);
//            //        Transient_Length.push_back(Measured.Transient_Length_Species);
//            //        Eqb_Minima.push_back(Measured.Eqb_Max_Minima_Species);
//            //        Transient_Minima.push_back(Measured.Transient_Max_Minima_Species);
//            //        Eqb_N_Clusters.push_back(Measured.Eqb_N_Clusters_Species);
//            //        Amplitude.push_back(Measured.Max_Amplitude_Species);

//            omp_unset_lock(&write_mutex);
//        }
//        else
//        {
//            Dynamics Measured;

//            Measured.Empty_Dynamics((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,(*Pass_to_Fortran::Thread_Local_Parameters_ptr).nPatches);

//            omp_set_lock(&write_mutex);

//            Individual_Community_Dynamics[iSet].Add_Dynamics(Measured);

//            //        std::vector<int> Dummy_dynamics ((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,0);
//            //        std::vector<double> Dummy_dynamics_too ((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,0);

//            //        Extinct.push_back(Which_Species_Extinct((*Pass_to_Fortran::Thread_Local_Parameters_ptr).nSpecies,sim_ana,Extinction_Threshold));
//            //        Transient_Length.push_back(Dummy_dynamics_too);
//            //        Eqb_Minima.push_back(Dummy_dynamics_too);
//            //        Transient_Minima.push_back(Dummy_dynamics_too);
//            //        Eqb_N_Clusters.push_back(Dummy_dynamics);
//            //        Amplitude.push_back(Dummy_dynamics_too);

//            omp_unset_lock(&write_mutex);

//        }
//    }

//    Pass_to_Fortran::Thread_Local_Parameters_ptr.reset(NULL);


//}
