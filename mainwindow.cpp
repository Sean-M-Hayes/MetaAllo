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
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "allometric_metacommunity.h"

// "parameter_generator_functions.h" from "mainwindow.h" and "allometric_metacommunity.h"
// <string> & "food_web_methods.h" from "parameter_generator_functions.h"
// <boost/multi_array.hpp> & <vector> from "food_web_methods.h"
// boost_matrix definition from "forward_declarations.h";

#include "table_validator_delegate.h"
#include "input_reading_functions.h"

#include "xml_methods.h"
// <QXmlReader> & <QTableWidget> from "xml_methods.h"

#include "read_csv.h"


#define BOOST_THREAD_USE_LIB
#include <boost/thread.hpp>
#include <time.h>
#include <iostream>

#include <boost/ref.hpp>

#include <functional>

#include <QFile>
#include <QMessageBox>

#include <QFileDialog>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    Multi_Food_Web.push_back(boost_matrix(boost::extents[1][1]));
    Multi_Food_Web[0][0][0] = 0;

    Multi_Spatial_Structure.push_back(boost_matrix(boost::extents[1][1]));
    Multi_Spatial_Structure[0][0][0] = 0;

//    QRegExp comma_delineated_int_regex("^\\d+(,\\d+)*$");
//    QRegExpValidator *comma_delineated_int_validator = new QRegExpValidator(comma_delineated_int_regex);

    QRegExp comma_delineated_double_regex("^\\d*\\.?\\d+(,\\d*\\.?\\d+)*$");
    QRegExpValidator *comma_delineated_double_validator = new QRegExpValidator(comma_delineated_double_regex);

    ui->Connectance_Input->setValidator(comma_delineated_double_validator);

    ui->Connectance_Spatial_Input->setValidator(comma_delineated_double_validator);

    ui->Branching_Input->setValidator(comma_delineated_double_validator);

    ui->nBranches_Input->setValidator(comma_delineated_double_validator);
    ui->Neighbor_Distance_Input->setValidator(comma_delineated_double_validator);

    ui->SS_N_Rewire_Input->setValidator(comma_delineated_double_validator);

    ui->ar_Input->setValidator(comma_delineated_double_validator);
    ui->ax_Input->setValidator(comma_delineated_double_validator);
    ui->e_Input->setValidator(comma_delineated_double_validator);
    ui->y_Input->setValidator(comma_delineated_double_validator);
    ui->B0_Input->setValidator(comma_delineated_double_validator);
    ui->h_Input->setValidator(comma_delineated_double_validator);
    ui->m_Input->setValidator(comma_delineated_double_validator);

    QRegExp comma_delineated_double_negative_regex("^\\-?\\d*\\.?\\d+(,\\-?\\d*\\.?\\d+)*$");
    QRegExpValidator *comma_delineated_double_negative_validator = new QRegExpValidator(comma_delineated_double_negative_regex);

    ui->Mass_trophic_Input->setValidator(comma_delineated_double_negative_validator);
    ui->m_Trophic_Input->setValidator(comma_delineated_double_negative_validator);
    ui->m_Competitive_Input->setValidator(comma_delineated_double_negative_validator);
    ui->m_Random_Input->setValidator(comma_delineated_double_negative_validator);

    QRegExp nThreads_regex("^[1-9]$");
    QRegExpValidator *nThreads_validator = new QRegExpValidator(nThreads_regex);

    ui->nThreads_Input->setValidator(nThreads_validator);

    QRegExp single_int_regex("^\\d+$");
    QRegExpValidator *single_int_validator = new QRegExpValidator(single_int_regex);

    ui->nIterations_Input->setValidator(single_int_validator);

    ui->MaxSimLength_Input->setValidator(single_int_validator);
    ui->MinSimLength_Input->setValidator(single_int_validator);
    ui->TestInt_Input->setValidator(single_int_validator);
    ui->LenBetweenTest_Input->setValidator(single_int_validator);

    ui->MinPatternObs_Input->setValidator(single_int_validator);
    ui->SeqRecenter_Input->setValidator(single_int_validator);
    ui->SeqNewCenter_Input->setValidator(single_int_validator);
    ui->SeqMinCenters_Input->setValidator(single_int_validator);
    ui->SeqMaxCenters_Input->setValidator(single_int_validator);

    QRegExp single_double_regex("^\\d*\\.?\\d+$");
    QRegExpValidator *single_double_validator = new QRegExpValidator(single_double_regex);

    ui->SimRes_Input->setValidator(single_double_validator);
    ui->nFood_Webs_Input->setValidator(single_double_validator);
    ui->nSpatial_Structures_Input->setValidator(single_double_validator);

    ui->CommonPatternFreq_Input->setValidator(single_double_validator);
    ui->CycleSwitch_Input->setValidator(single_double_validator);
    ui->DiffErrTol_Input->setValidator(single_double_validator);
    ui->MinExtremaPropDiff_Input->setValidator(single_double_validator);
    ui->SeqCenterCV_Input->setValidator(single_double_validator);
    ui->SeqThreshIncCenters_Input->setValidator(single_double_validator);
    ui->SeqConvergence_Input->setValidator(single_double_validator);
    ui->SeqRecenter_Input->setValidator(single_double_validator);
    ui->SeqNewCenter_Input->setValidator(single_double_validator);
    ui->SeqMinCenters_Input->setValidator(single_double_validator);
    ui->SeqMaxCenters_Input->setValidator(single_double_validator);

    ui->StationaryCycleFreq_Input->setValidator(single_double_validator);
    ui->StationaryCycleProp_Input->setValidator(single_double_validator);

    ui->Extinction_Input->setValidator(comma_delineated_double_validator);

    //Table Validators

    ui->Food_Web_Multi_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Food_Web_Multi_Input));
    ui->Spatial_Structure_Multi_Input->setItemDelegate(new TableValidatorDelegate_Int(ui->Spatial_Structure_Multi_Input));
    ui->Body_Mass_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Body_Mass_Input));
    ui->Metabolic_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Metabolic_Input));
    ui->Functional_Response_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Functional_Response_Input));
    ui->Dispersal_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Dispersal_Input));
    ui->Carrying_Capacity_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Carrying_Capacity_Input));
    ui->Initial_Abundance_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Initial_Abundance_Input));
    ui->Patch_Size_Input->setItemDelegate(new TableValidatorDelegate_Double(ui->Patch_Size_Input));

    QMessageBox::information(this,"Copyright Information",
                          QString(
                             "Copyright 2017 Sean M. Hayes \n\nThis program is free software: you can redistribute it and/or modify \nit under the terms of the GNU General Public License as published by \nthe Free Software Foundation, either version 3 of the License, or \n(at your option) any later version. \n\nThis program is distributed in the hope that it will be useful, \nbut WITHOUT ANY WARRANTY; without even the implied warranty of \nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \nGNU General Public License for more details. \n\nYou should have received a copy of the GNU General Public License \nalong with this program.  If not, see <http://www.gnu.org/licenses/>.\n\nFull source code is available at <https://github.com/SMHayes/MetaAllo>"),
                          QMessageBox::Ok);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::Food_Web_Method()
{   
    Food_Web_Generator * Food_Web_Parameter_Method = NULL;

    switch(ui->Food_Web_Method->currentIndex())
    {
    case 0:
    {
        Food_Web_Parameter_Method = new Food_Web_Read_Generator(&Input_Cache,Multi_Food_Web);
    }
        break;
    case 1:
    {
        int nSpecies;
        std::vector<double> Connectance;

        if(Collect_Input(ui->Connectance_Input,Connectance))
            if(Collect_Input(ui->Species_Input,nSpecies))
                Food_Web_Parameter_Method = new Food_Web_Niche_Model_Generator(&Input_Cache,nSpecies,Connectance);
    }
        break;
    }

    if(!(Input_Cache.Food_Web==NULL))
        delete Input_Cache.Food_Web;

    Input_Cache.Food_Web = Food_Web_Parameter_Method;
}

void MainWindow::Spatial_Structure_Method()
{
    Spatial_Structure_Generator * Spatial_Structure_Parameter_Method = NULL;

    switch(ui->Spatial_Structure_Method->currentIndex())
    {
    case 0:
    {
        Spatial_Structure_Parameter_Method = new Spatial_Structure_Read_Generator(&Input_Cache,Multi_Spatial_Structure);
    }
        break;
    case 1:
    {
        int nPatches;
        std::vector<double> Connectance;

        if(Collect_Input(ui->Connectance_Spatial_Input,Connectance))
            if(Collect_Input(ui->Patches_Input,nPatches))
                Spatial_Structure_Parameter_Method = new Spatial_Structure_Undirected_Erdos_Renyi_Generator(&Input_Cache,nPatches,Connectance);
    }
        break;
    case 2:
    {
        int nPatches;
        std::vector<double> Branching_Probability;

        if(Collect_Input(ui->Branching_Input,Branching_Probability))
            if(Collect_Input(ui->Patches_Input,nPatches))
                Spatial_Structure_Parameter_Method = new Spatial_Structure_Undirected_Dendritic_Generator(&Input_Cache,nPatches,Branching_Probability);
    }
        break;
    case 3:
    {
        int nPatches;
        std::vector<int> nBranches;
        std::vector<int> Neighbor_Distance;

        if(Collect_Input(ui->nBranches_Input,nBranches))
            if(Collect_Input(ui->Neighbor_Distance_Input,Neighbor_Distance))
                if(Collect_Input(ui->Patches_Input,nPatches))
                    Spatial_Structure_Parameter_Method = new Spatial_Structure_Kurt_Comparison_Generator(&Input_Cache,nPatches,nBranches,Neighbor_Distance);
    }
        break;
    case 4:
    {
        int nPatches;
        std::vector<int> nRewire;

        if(Collect_Input(ui->SS_N_Rewire_Input,nRewire))
            if(Collect_Input(ui->Patches_Input,nPatches))
                if(nPatches==18) //Shipley-Skinner only works with 19 patches
                    Spatial_Structure_Parameter_Method = new Spatial_Structure_Shipley_Skinner_Generator(&Input_Cache,nRewire);
    }
        break;
    }

    if(!(Input_Cache.Spatial_Structure==NULL))
        delete Input_Cache.Spatial_Structure;

    Input_Cache.Spatial_Structure = Spatial_Structure_Parameter_Method;
}

void MainWindow::Body_Mass_Method()
{
    Body_Mass_Generator * Body_Mass_Parameter_Method = NULL;

    switch(ui->Body_Mass_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> Body_Mass;

        if(Collect_Input_Table_Column(0,ui->Body_Mass_Input,Body_Mass))
            Body_Mass_Parameter_Method = new Body_Mass_Generator(&Input_Cache,Body_Mass);
    }
        break;
    case 1:
    {
        std::vector<double> Input_Base_Mass;
        std::vector<double> Input_Trophic_Scaling;
        std::vector<double> Input_Mass_Variation;

        if(!(Input_Cache.Food_Web==NULL))
            if(Collect_Input(ui->Base_mass_Input,Input_Base_Mass))
                if(Collect_Input(ui->Mass_trophic_Input,Input_Trophic_Scaling))
                    if(Collect_Input(ui->Mass_variation_Input,Input_Mass_Variation))
                        Body_Mass_Parameter_Method = new Body_Mass_Trophic_Scaling_Generator(&Input_Cache,ui->Body_Mass_Scaling_Input->currentIndex(),Input_Base_Mass,Input_Trophic_Scaling,Input_Mass_Variation);

    }
        break;
    }

    if(!(Input_Cache.Body_Mass==NULL))
        delete Input_Cache.Body_Mass;

    Input_Cache.Body_Mass = Body_Mass_Parameter_Method;
}


void MainWindow::Metabolism_Method()
{
    Metabolism_Generator * Metabolism_Parameter_Method = NULL;

    switch(ui->Metabolism_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> ar;
        std::vector<double> ax;
        std::vector<double> y;
        std::vector<double> e;

        if(Collect_Input_Table_Column(0,ui->Metabolic_Input,ar))
            if(Collect_Input_Table_Column(1,ui->Metabolic_Input,ax))
                if(Collect_Input_Table_Column(2,ui->Metabolic_Input,y))
                    if(Collect_Input_Table_Column(3,ui->Metabolic_Input,e))
                        Metabolism_Parameter_Method = new Metabolism_Generator(&Input_Cache,ar,ax,y,e);
    }
        break;
    case 1:
    {
        int nSpecies;

        std::vector<double> ar_fixed;
        std::vector<double> ax_fixed;
        std::vector<double> y_fixed;
        std::vector<double> e_fixed;

        if(Collect_Input(ui->Species_Input,nSpecies))
            if(Collect_Input(ui->ar_Input,ar_fixed))
                if(Collect_Input(ui->ax_Input,ax_fixed))
                    if(Collect_Input(ui->y_Input,y_fixed))
                        if(Collect_Input(ui->e_Input,e_fixed))
                            Metabolism_Parameter_Method = new Metabolism_Fixed_Species_Generator(&Input_Cache,nSpecies,ar_fixed,ax_fixed,y_fixed,e_fixed);

    }
        break;
    }

    if(!(Input_Cache.Metabolism==NULL))
        delete Input_Cache.Metabolism;

    Input_Cache.Metabolism = Metabolism_Parameter_Method;
}

void MainWindow::Functional_Response_Method()
{
    Functional_Response_Generator * Functional_Response_Parameter_Method = NULL;

    switch(ui->Functional_Response_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> B0;
        std::vector<double> h;

        if(Collect_Input_Table_Column(0,ui->Functional_Response_Input,B0))
            if(Collect_Input_Table_Column(1,ui->Functional_Response_Input,h))
                Functional_Response_Parameter_Method = new Functional_Response_Generator(&Input_Cache,B0,h);
    }
        break;
    case 1:
    {
        int nSpecies;

        std::vector<double> B0_fixed;
        std::vector<double> h_fixed;

        if(Collect_Input(ui->Species_Input,nSpecies))
            if(Collect_Input(ui->B0_Input,B0_fixed))
                if(Collect_Input(ui->h_Input,h_fixed))
                    Functional_Response_Parameter_Method = new Functional_Response_Fixed_Species_Generator(&Input_Cache,nSpecies,B0_fixed,h_fixed);
    }
        break;
    }

    if(!(Input_Cache.Functional_Response==NULL))
        delete Input_Cache.Functional_Response;

    Input_Cache.Functional_Response = Functional_Response_Parameter_Method;

}

void MainWindow::Dispersal_Method()
{
    Dispersal_Generator * Dispersal_Parameter_Method = NULL;

    switch(ui->Dispersal_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> Dispersal;

        if(Collect_Input_Table_Column(0,ui->Dispersal_Input,Dispersal))
            Dispersal_Parameter_Method = new Dispersal_Generator(&Input_Cache,Dispersal);
    }
        break;
    case 1:
    {
        std::vector<double> Dispersal_Base;
        std::vector<double> Trophic_Scaling;
        std::vector<double> Competitive_Scaling;
        std::vector<double> Random_Scaling;

        if(!(Input_Cache.Food_Web==NULL))
            if(Collect_Input(ui->m_Input,Dispersal_Base))
                if(Collect_Input(ui->m_Trophic_Input,Trophic_Scaling))
                    if(Collect_Input(ui->m_Competitive_Input,Competitive_Scaling))
                        if(Collect_Input(ui->m_Random_Input,Random_Scaling))
                            Dispersal_Parameter_Method = new Dispersal_Allometric_Linear_Scaling_Generator(&Input_Cache,Dispersal_Base,Trophic_Scaling,Competitive_Scaling,Random_Scaling);
    }
        break;
    case 2:
    {
        std::vector<double> Dispersal_Base;
        std::vector<double> Trophic_Scaling;
        std::vector<double> Competitive_Scaling;
        std::vector<double> Random_Scaling;

        if(!(Input_Cache.Food_Web==NULL))
            if(Collect_Input(ui->m_Input_2,Dispersal_Base))
                if(Collect_Input(ui->m_Trophic_Input_2,Trophic_Scaling))
                    if(Collect_Input(ui->m_Competitive_Input_2,Competitive_Scaling))
                        if(Collect_Input(ui->m_Random_Input_2,Random_Scaling))
                            Dispersal_Parameter_Method = new Dispersal_Allometric_Exponential_Scaling_Generator(&Input_Cache,Dispersal_Base,Trophic_Scaling,Competitive_Scaling,Random_Scaling);
    }
        break;
    }

    if(!(Input_Cache.Dispersal==NULL))
        delete Input_Cache.Dispersal;

    Input_Cache.Dispersal = Dispersal_Parameter_Method;
}

void MainWindow::Carrying_Capacity_Method()
{
    Carrying_Capacity_Generator * Carrying_Capacity_Parameter_Method = NULL;

    switch(ui->Carrying_Capacity_Method->currentIndex())
    {
    case 0:
    {
        boost_matrix Carrying_Capacity;
        std::vector<double> Patch_Size;

        if(Collect_Input_Table(ui->Carrying_Capacity_Input,Carrying_Capacity))
            if(Collect_Input_Table_Row(0,ui->Patch_Size_Input,Patch_Size))
                Carrying_Capacity_Parameter_Method = new Carrying_Capacity_Generator(&Input_Cache,Carrying_Capacity,Patch_Size);

    }
        break;
    case 1:
    {
        int nSpecies;
        int nPatches;

        std::vector<double> Carrying_Capacity_Fixed_Community;

        if(Collect_Input(ui->Species_Input,nSpecies))
            if(Collect_Input(ui->Patches_Input,nPatches))
                if(Collect_Input(ui->Carrying_Capacity_Fixed_Community_Input,Carrying_Capacity_Fixed_Community))
                    Carrying_Capacity_Parameter_Method = new Carrying_Capacity_Fixed_Community_Generator(&Input_Cache,nSpecies,nPatches,Carrying_Capacity_Fixed_Community);
    }
        break;
    case 2:
    {
        int nSpecies;

        boost_matrix Patch_Size;

        if(Collect_Input(ui->Species_Input,nSpecies))
            if(Collect_Input_Table(ui->Patch_Size_Input_2,Patch_Size))
                Carrying_Capacity_Parameter_Method = new Carrying_Capacity_Scale_By_Size_Generator(&Input_Cache,nSpecies,Patch_Size);
    }
        break;
    }

    if(!(Input_Cache.Carrying_Capacity==NULL))
        delete Input_Cache.Carrying_Capacity;

    Input_Cache.Carrying_Capacity = Carrying_Capacity_Parameter_Method;
}

void MainWindow::Initial_Abundance_Method()
{
    Initial_Abundance_Generator * Initial_Abundance_Parameter_Method = NULL;

    switch(ui->Initial_Abundance_Method->currentIndex())
    {
    case 0:
    {
        boost_matrix Initial_Abundance;

        if(Collect_Input_Table(ui->Initial_Abundance_Input,Initial_Abundance))
            Initial_Abundance_Parameter_Method = new Initial_Abundance_Generator(&Input_Cache,Initial_Abundance);
    }
        break;
    case 1:
    {
        int nPatches;

        std::vector<double> Initial_Abundance_Mean;
        std::vector<double> Initial_Abundance_Dev;

        if(Collect_Input(ui->Patches_Input,nPatches))
            if(Collect_Input_Table_Column(0,ui->Initial_Abundance_Random_Mean_Input,Initial_Abundance_Mean))
                if(Collect_Input_Table_Column(0,ui->Initial_Abundance_Random_Dev_Input,Initial_Abundance_Dev))
                    Initial_Abundance_Parameter_Method = new Initial_Abundance_Random_Generator(&Input_Cache,nPatches,Initial_Abundance_Mean,Initial_Abundance_Dev);
    }
        break;
    }

    if(!(Input_Cache.Initial_Abundance==NULL))
        delete Input_Cache.Initial_Abundance;

    Input_Cache.Initial_Abundance = Initial_Abundance_Parameter_Method;
}

void MainWindow::Extinction_Threshold_Method()
{
    Extinction_Threshold_Generator * Extinction_Threshold_Method = NULL;

    //Only does sets, default is unnecessary

    std::vector<double> Extinction_Threshold_Set;

    if(Collect_Input(ui->Extinction_Input,Extinction_Threshold_Set))
        Extinction_Threshold_Method = new Extinction_Threshold_Set_Generator(&Input_Cache,Extinction_Threshold_Set);

    if(!(Input_Cache.Extinction_Threshold==NULL))
        delete Input_Cache.Extinction_Threshold;

    Input_Cache.Extinction_Threshold = Extinction_Threshold_Method;
}

void MainWindow::on_Simulate_clicked()
{
    ui->Simulate->setEnabled(false);

    Simulation_Properties Sim_Parms;

    bool Sim_Methods_Check = false;

    if(Collect_Input(ui->nThreads_Input,Sim_Parms.nThreads))
        if(Collect_Input(ui->nIterations_Input,Sim_Parms.nIterations))
            if(Collect_Input(ui->nFood_Webs_Input,Sim_Parms.nFood_Webs))
                if(Collect_Input(ui->nSpatial_Structures_Input,Sim_Parms.nSpatial_Structures))
                    if(Collect_Input(ui->MaxSimLength_Input,Sim_Parms.Maximum_Simulation_Duration))
                        if(Collect_Input(ui->MinSimLength_Input,Sim_Parms.Minimum_Simulation_Duration))
                            if(Collect_Input(ui->TestInt_Input,Sim_Parms.Equilibrium_Test_Interval_Duration))
                                if(Collect_Input(ui->LenBetweenTest_Input,Sim_Parms.Time_between_Equilibrium_Tests))
                                    if(Collect_Input(ui->CommonPatternFreq_Input,Sim_Parms.Cycle_Params.Min_Prop_for_Common_Pattern))
                                        if(Collect_Input(ui->MinPatternObs_Input,Sim_Parms.Cycle_Params.Min_Occs_for_Common_Pattern))
                                            if(Collect_Input(ui->CycleSwitch_Input,Sim_Parms.Cycle_Params.Max_Cycle_Switch_Threshold))
                                                if(Collect_Input(ui->DiffErrTol_Input,Sim_Parms.Cycle_Params.Difference_Error_Tolerance))
                                                    if(Collect_Input(ui->MinExtremaPropDiff_Input,Sim_Parms.Cycle_Params.Minimum_Extrema_Proportion_Difference))
                                                        if(Collect_Input(ui->SeqCenterCV_Input,Sim_Parms.Cycle_Params.Sequence_Center_Target_CV))
                                                            if(Collect_Input(ui->SeqThreshIncCenters_Input,Sim_Parms.Cycle_Params.Sequence_Threshold_to_Increase_Centers))
                                                                if(Collect_Input(ui->SeqConvergence_Input,Sim_Parms.Cycle_Params.Sequence_Convergence_Criteria))
                                                                    if(Collect_Input(ui->SeqRecenter_Input,Sim_Parms.Cycle_Params.Sequence_Recenter_Iterations))
                                                                        if(Collect_Input(ui->SeqNewCenter_Input,Sim_Parms.Cycle_Params.Sequence_New_Center_Iterations))
                                                                            if(Collect_Input(ui->SeqMinCenters_Input,Sim_Parms.Cycle_Params.Sequence_Min_Centers))
                                                                                if(Collect_Input(ui->SeqMaxCenters_Input,Sim_Parms.Cycle_Params.Sequence_Max_Centers))
                                                                                    if(Collect_Input(ui->StationaryCycleFreq_Input,Sim_Parms.Cycle_Params.Stationary_Cycle_Frequency_Threshold))
                                                                                        if(Collect_Input(ui->StationaryCycleProp_Input,Sim_Parms.Cycle_Params.Proportion_Series_for_Stationary))
                                                                                            if(Collect_Input(ui->RoundTimesTo_Input,Sim_Parms.Cycle_Params.Round_Times_To))
                                                                                                if(Collect_Input(ui->SimRes_Input,Sim_Parms.Simulation_Resolution))
                                                                                                    if(Collect_Input(ui->Extinction_Input,Sim_Parms.Extinction_Threshold))
                                                                                                        if(Collect_Input(ui->PhaseLockFreq_Input,Sim_Parms.Phase_Lock_Threshold))
                                                                                                        {

                                                                                                        Sim_Parms.Time_Series_Write = ui->Write_Combo->currentIndex();
                                                                                                        Sim_Parms.Limit_Cycle_Write = ui->Limit_Write_Combo->currentIndex();
                                                                                                        Sim_Parms.Species_Data_Write = ui->Write_Species_Data_Input->currentIndex();
                                                                                                        Sim_Parms.Write_Each_Set = ui->WritePerSet_Input->isChecked();
                                                                                                        Sim_Parms.Test_against_no_Structure = ui->TestNullStructure_Input->isChecked();
                                                                                                        Sim_Parms.Write_Measure_Diagnostics = ui->WriteDiagnostic_Input->isChecked();
                                                                                                        Sim_Parms.Write_Path = ui->Write_Path_Input->text().toStdString()+"/";

                                                                                                        Sim_Methods_Check = true;

                                                                                                        //Parameter_Generators Input;

                                                                                                        Food_Web_Method();
                                                                                                        Spatial_Structure_Method();
                                                                                                        Body_Mass_Method();
                                                                                                        Metabolism_Method();
                                                                                                        Functional_Response_Method();
                                                                                                        Dispersal_Method();
                                                                                                        Carrying_Capacity_Method();
                                                                                                        Initial_Abundance_Method();
                                                                                                        Extinction_Threshold_Method();

                                                                                                        QString Input_Error_Check = QString::fromStdString(Input_Cache.Check());

                                                                                                        if(Input_Error_Check.isEmpty())
                                                                                                        {
                                                                                                            ui->Simulation_Progress->setEnabled(true);

                                                                                                            Meta_Allo_Simulation Meta_Allo(Parameter_Generators(Input_Cache),Sim_Parms);

                                                                                                            ui->Simulation_Progress->setMaximum(Meta_Allo.nSimulations);

                                                                                                            boost::thread Sim_Thread(boost::bind(&Meta_Allo_Simulation::Run_Careful_Simulation,&Meta_Allo));

                                                                                                            while(Meta_Allo.counter<Meta_Allo.nSimulations)
                                                                                                            {
                                                                                                                ui->Simulation_Progress->setValue(Meta_Allo.counter);
                                                                                                                ui->Simulation_Progress->repaint();
                                                                                                                QCoreApplication::processEvents();
                                                                                                            }

                                                                                                            Sim_Thread.join();

                                                                                                            ui->Simulation_Progress->setValue(Meta_Allo.counter);
                                                                                                            ui->Simulation_Progress->repaint();

                                                                                                        }
                                                                                                        else
                                                                                                        {
                                                                                                            //Error Popup
                                                                                                            QMessageBox::critical(this,"Simulation Initilization Problem",
                                                                                                                                  Input_Error_Check,
                                                                                                                                  QMessageBox::Ok);
                                                                                                        }
                                                                                                    }
    if(!Sim_Methods_Check)
    {
        //Error Popup
        QMessageBox::critical(this,"Simulation Initilization Problem",
                              "Simulation Setup Parameters Missing",
                              QMessageBox::Ok);
    }

    ui->Simulate->setEnabled(true);
}

void Adjust_Table_Size(QTableWidget *Table_Widget,int new_rows,int new_cols)
{
    while(new_rows>Table_Widget->rowCount())
    {
       Table_Widget->insertRow(Table_Widget->rowCount());
    }
    while(new_rows<Table_Widget->rowCount())
    {
       Table_Widget->removeRow(Table_Widget->rowCount()-1);
    }
    while(new_cols>Table_Widget->columnCount())
    {
       Table_Widget->insertColumn(Table_Widget->columnCount());
    }
    while(new_cols<Table_Widget->columnCount())
    {
       Table_Widget->removeColumn(Table_Widget->columnCount()-1);
    }
}

void MainWindow::on_Write_Input_Prompt_clicked()
{
    QString Write_Directory = QFileDialog::getExistingDirectory(this,
                                                                tr("Write to Directory"),".");

    ui->Write_Path_Input->setText(Write_Directory);
}

//
//
//

void MainWindow::on_Food_Web_Read_Input_Prompt_clicked()
{
    QString Write_Directory = QFileDialog::getOpenFileName(this,
                                                           tr("Read from File"),".");

    ui->Food_Web_Read_Input->setText(Write_Directory);

    std::string Food_Web_Read_File;

    if(Collect_Input(ui->Food_Web_Read_Input,Food_Web_Read_File))
    {
        if(Read_CSV_Matrices(Food_Web_Read_File,Multi_Food_Web,ui->Species_Input->value(),ui->Species_Input->value()))
        {
            //ui->Species_Input->setValue(Multi_Food_Web[0].shape()[0]); //Sets row & colum count of Food_Web_Multi_Input

            for(int i=0;i<Multi_Food_Web[0].shape()[0];i++)
                for(int j=0;j<Multi_Food_Web[0].shape()[1];j++)
                    ui->Food_Web_Multi_Input->setItem(i,j,new QTableWidgetItem(QString::number(Multi_Food_Web[0][i][j])));

            ui->Multi_Food_Web_Spin->setValue(1);
            ui->Multi_Food_Web_Spin->setEnabled(true);
            ui->Multi_Food_Web_Spin->setMaximum(Multi_Food_Web.size());
        }
        else
        {
            QMessageBox::critical(this,"Read File Problem",
                                  "Couldn't read file",
                                  QMessageBox::Ok);
        }
    }
    else
    {
        QMessageBox::critical(this,"Read File Problem",
                              "Couldn't read file",
                              QMessageBox::Ok);
    }


}

void MainWindow::on_Multi_Food_Web_Spin_valueChanged(int arg1)
{
    for(int i=0;i<ui->Food_Web_Multi_Input->rowCount();i++)
    {
        if(Multi_Food_Web[arg1-1].shape()[0]>i)
        {
            for(int j=0;j<ui->Food_Web_Multi_Input->columnCount();j++)
            {
                if(Multi_Food_Web[arg1-1].shape()[1]>j)
                {
                    ui->Food_Web_Multi_Input->setItem(i,j,new QTableWidgetItem(QString::number(Multi_Food_Web[arg1-1][i][j])));
                }
            }
        }
    }
}

void MainWindow::on_Food_Web_Multi_Input_itemChanged(QTableWidgetItem *item)
{
    Multi_Food_Web[ui->Multi_Food_Web_Spin->value()-1][item->row()][item->column()] = item->text().toDouble();
}

void MainWindow::on_Food_Web_Add_clicked()
{
    Multi_Food_Web.push_back(boost_matrix(boost::extents[ui->Species_Input->value()][ui->Species_Input->value()]));

    for(int i=0;i<Multi_Food_Web.back().shape()[0];i++)
        for(int j=0;j<Multi_Food_Web.back().shape()[1];j++)
            Multi_Food_Web.back()[i][j] = 0;

    ui->Multi_Food_Web_Spin->setMaximum(Multi_Food_Web.size());
    ui->Multi_Food_Web_Spin->setValue(Multi_Food_Web.size());
}

void MainWindow::on_Food_Web_Remove_clicked()
{
    if(Multi_Food_Web.size()>1)
    {
        Multi_Food_Web.erase(Multi_Food_Web.begin()+ui->Multi_Food_Web_Spin->value()-1);

        if(ui->Multi_Food_Web_Spin->value()>Multi_Food_Web.size())
            ui->Multi_Food_Web_Spin->setValue(Multi_Food_Web.size());
        else
            on_Multi_Food_Web_Spin_valueChanged(ui->Multi_Food_Web_Spin->value());

        ui->Multi_Food_Web_Spin->setMaximum(Multi_Food_Web.size());
    }
}

//
//
//

void MainWindow::on_Spatial_Structure_Read_Input_Prompt_clicked()
{
    QString Write_Directory = QFileDialog::getOpenFileName(this,
                                                           tr("Read from File"),".");

    ui->Spatial_Structure_Read_Input->setText(Write_Directory);

    std::string Spatial_Structure_Read_File;

    if(Collect_Input(ui->Spatial_Structure_Read_Input,Spatial_Structure_Read_File))
    {
        if(Read_CSV_Matrices(Spatial_Structure_Read_File,Multi_Spatial_Structure,ui->Patches_Input->value(),ui->Patches_Input->value()))
        {
            //ui->Patches_Input->setValue(Multi_Spatial_Structure[0].shape()[0]); //Sets row & colum count of Spatial_Structure_Multi_Input

            for(int i=0;i<Multi_Spatial_Structure[0].shape()[0];i++)
                for(int j=0;j<Multi_Spatial_Structure[0].shape()[1];j++)
                    ui->Spatial_Structure_Multi_Input->setItem(i,j,new QTableWidgetItem(QString::number(Multi_Spatial_Structure[0][i][j])));

            ui->Multi_Spatial_Structure_Spin->setValue(1);
            ui->Multi_Spatial_Structure_Spin->setEnabled(true);
            ui->Multi_Spatial_Structure_Spin->setMaximum(Multi_Spatial_Structure.size());
        }
        else
        {
            QMessageBox::critical(this,"Read File Problem",
                                  "Couldn't read file",
                                  QMessageBox::Ok);
        }
    }
    else
    {
        QMessageBox::critical(this,"Read File Problem",
                              "Couldn't read file",
                              QMessageBox::Ok);
    }


}

void MainWindow::on_Multi_Spatial_Structure_Spin_valueChanged(int arg1)
{
    for(int i=0;i<ui->Spatial_Structure_Multi_Input->rowCount();i++)
    {
        if(Multi_Spatial_Structure[arg1-1].shape()[0]>i)
        {
            for(int j=0;j<ui->Spatial_Structure_Multi_Input->columnCount();j++)
            {
                if(Multi_Spatial_Structure[arg1-1].shape()[1]>j)
                {
                    ui->Spatial_Structure_Multi_Input->setItem(i,j,new QTableWidgetItem(QString::number(Multi_Spatial_Structure[arg1-1][i][j])));
                }
            }
        }
    }
}

void MainWindow::on_Spatial_Structure_Multi_Input_itemChanged(QTableWidgetItem *item)
{
    Multi_Spatial_Structure[ui->Multi_Spatial_Structure_Spin->value()-1][item->row()][item->column()] = item->text().toDouble();
}

void MainWindow::on_Spatial_Structure_Add_clicked()
{
    Multi_Spatial_Structure.push_back(boost_matrix(boost::extents[ui->Patches_Input->value()][ui->Patches_Input->value()]));

    for(int i=0;i<Multi_Spatial_Structure.back().shape()[0];i++)
        for(int j=0;j<Multi_Spatial_Structure.back().shape()[1];j++)
            Multi_Spatial_Structure.back()[i][j] = 0;

    ui->Multi_Spatial_Structure_Spin->setMaximum(Multi_Spatial_Structure.size());
    ui->Multi_Spatial_Structure_Spin->setValue(Multi_Spatial_Structure.size());
}

void MainWindow::on_Spatial_Structure_Remove_clicked()
{
    if(Multi_Spatial_Structure.size()>1)
    {
        Multi_Spatial_Structure.erase(Multi_Spatial_Structure.begin()+ui->Multi_Spatial_Structure_Spin->value()-1);

        if(ui->Multi_Spatial_Structure_Spin->value()>Multi_Spatial_Structure.size())
            ui->Multi_Spatial_Structure_Spin->setValue(Multi_Spatial_Structure.size());
        else
            on_Multi_Spatial_Structure_Spin_valueChanged(ui->Multi_Spatial_Structure_Spin->value());

        ui->Multi_Spatial_Structure_Spin->setMaximum(Multi_Spatial_Structure.size());
    }
}

//
//
//

void MainWindow::on_Species_Input_valueChanged(int arg1)
{
    for(int i=0;i<Multi_Food_Web.size();i++)
        Multi_Food_Web[i].resize(boost::extents[arg1][arg1]);

    Adjust_Table_Size(ui->Food_Web_Multi_Input,arg1,arg1);
    Adjust_Table_Size(ui->Body_Mass_Input,arg1,ui->Body_Mass_Input->columnCount());
    Adjust_Table_Size(ui->Metabolic_Input,arg1,ui->Metabolic_Input->columnCount());
    Adjust_Table_Size(ui->Functional_Response_Input,arg1,ui->Functional_Response_Input->columnCount());
    Adjust_Table_Size(ui->Dispersal_Input,arg1,ui->Functional_Response_Input->columnCount());
    Adjust_Table_Size(ui->Carrying_Capacity_Input,arg1,ui->Carrying_Capacity_Input->columnCount());
    Adjust_Table_Size(ui->Initial_Abundance_Input,arg1,ui->Initial_Abundance_Input->columnCount());
    Adjust_Table_Size(ui->Initial_Abundance_Random_Mean_Input,arg1,ui->Initial_Abundance_Random_Mean_Input->columnCount());
    Adjust_Table_Size(ui->Initial_Abundance_Random_Dev_Input,arg1,ui->Initial_Abundance_Random_Dev_Input->columnCount());
}

void MainWindow::on_Patches_Input_valueChanged(int arg1)
{
    for(int i=0;i<Multi_Spatial_Structure.size();i++)
        Multi_Spatial_Structure[i].resize(boost::extents[arg1][arg1]);

    //Adjust_Table_Size(ui->Spatial_Structure_Multi_Input,arg1,arg1);
    Adjust_Table_Size(ui->Spatial_Structure_Multi_Input,arg1,arg1);
    Adjust_Table_Size(ui->Carrying_Capacity_Input,ui->Carrying_Capacity_Input->rowCount(),arg1);
    Adjust_Table_Size(ui->Initial_Abundance_Input,ui->Initial_Abundance_Input->rowCount(),arg1);
    Adjust_Table_Size(ui->Patch_Size_Input,ui->Patch_Size_Input->rowCount(),arg1);
    Adjust_Table_Size(ui->Patch_Size_Input_2,ui->Patch_Size_Input_2->rowCount(),arg1);
}

void MainWindow::on_Scale_by_Size_Iterates_Input_valueChanged(int arg1)
{
    Adjust_Table_Size(ui->Patch_Size_Input_2,arg1,ui->Patch_Size_Input_2->columnCount());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write XML Values

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::on_actionSave_triggered()
{
    //QString Save_File = ui->Default_Input->text();

    QString Save_File = QFileDialog::getSaveFileName(this,
                                           tr("Save Xml"), ".",
                                           tr("Xml files (*.xml)"));

    QFile Write_File(Save_File);

    //if (!Defaults_File.open(QIODevice::ReadWrite | QIODevice::Text))
    if (!Write_File.open(QIODevice::WriteOnly))
    {
        QMessageBox::critical(this,"Write XML File Problem",
                              "Couldn't write file",
                              QMessageBox::Ok);
        return;
    }

    QXmlStreamWriter xmlWriter(&Write_File);
    xmlWriter.setAutoFormatting(true);
    xmlWriter.writeStartDocument();

    xmlWriter.writeStartElement("Allometric_Metacommunity");

    xmlWriter.writeTextElement("nSpecies",ui->Species_Input->text());
    xmlWriter.writeTextElement("nPatches",ui->Patches_Input->text());

    xmlWriter.writeStartElement("Simulation_Properties");
    {

    Simulation_Properties Dummy_Sim_Parms;
;

    if(Collect_Input(ui->MinSimLength_Input,Dummy_Sim_Parms.Minimum_Simulation_Duration))
        xmlWriter.writeTextElement("Min_Simulation_Length",QString::number(Dummy_Sim_Parms.Minimum_Simulation_Duration));

    if(Collect_Input(ui->MaxSimLength_Input,Dummy_Sim_Parms.Maximum_Simulation_Duration))
        xmlWriter.writeTextElement("Max_Simulation_Length",QString::number(Dummy_Sim_Parms.Maximum_Simulation_Duration));

    if(Collect_Input(ui->TestInt_Input,Dummy_Sim_Parms.Equilibrium_Test_Interval_Duration))
        xmlWriter.writeTextElement("Test_Interval_Length",QString::number(Dummy_Sim_Parms.Equilibrium_Test_Interval_Duration));

    if(Collect_Input(ui->LenBetweenTest_Input,Dummy_Sim_Parms.Time_between_Equilibrium_Tests))
        xmlWriter.writeTextElement("Length_Between_Test_Intervals",QString::number(Dummy_Sim_Parms.Time_between_Equilibrium_Tests));

    if(Collect_Input(ui->nIterations_Input,Dummy_Sim_Parms.nIterations))
        xmlWriter.writeTextElement("Simulation_Iterations",QString::number(Dummy_Sim_Parms.nIterations));

    if(Collect_Input(ui->SimRes_Input,Dummy_Sim_Parms.Simulation_Resolution))
        xmlWriter.writeTextElement("Simulation_Resolution",QString::number(Dummy_Sim_Parms.Simulation_Resolution));

    if(Collect_Input(ui->nThreads_Input,Dummy_Sim_Parms.nThreads))
        xmlWriter.writeTextElement("Simulation_Threads",QString::number(Dummy_Sim_Parms.nThreads));

    //

    if(Collect_Input(ui->PhaseLockFreq_Input,Dummy_Sim_Parms.Phase_Lock_Threshold))
        xmlWriter.writeTextElement("Phase_Lock_Threshold",QString::number(Dummy_Sim_Parms.Phase_Lock_Threshold));

    //

    if(Collect_Input(ui->StationaryCycleFreq_Input,Dummy_Sim_Parms.Cycle_Params.Stationary_Cycle_Frequency_Threshold))
        xmlWriter.writeTextElement("Stationary_Cycle_Frequency",QString::number(Dummy_Sim_Parms.Cycle_Params.Stationary_Cycle_Frequency_Threshold));

    if(Collect_Input(ui->StationaryCycleProp_Input,Dummy_Sim_Parms.Cycle_Params.Proportion_Series_for_Stationary))
        xmlWriter.writeTextElement("Proportion_Stationary_Cycles",QString::number(Dummy_Sim_Parms.Cycle_Params.Proportion_Series_for_Stationary));

    if(Collect_Input(ui->DiffErrTol_Input,Dummy_Sim_Parms.Cycle_Params.Difference_Error_Tolerance))
        xmlWriter.writeTextElement("Difference_Error_Tolerance",QString::number(Dummy_Sim_Parms.Cycle_Params.Difference_Error_Tolerance));

    if(Collect_Input(ui->MinExtremaPropDiff_Input,Dummy_Sim_Parms.Cycle_Params.Minimum_Extrema_Proportion_Difference))
        xmlWriter.writeTextElement("Minimum_Extrema_Proportion_Difference",QString::number(Dummy_Sim_Parms.Cycle_Params.Minimum_Extrema_Proportion_Difference));

    if(Collect_Input(ui->SeqCenterCV_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_Center_Target_CV))
        xmlWriter.writeTextElement("Sequence_Center_Target_CV",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_Center_Target_CV));

    if(Collect_Input(ui->SeqThreshIncCenters_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_Threshold_to_Increase_Centers))
        xmlWriter.writeTextElement("Sequence_Threshold_to_Increase_Centers",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_Threshold_to_Increase_Centers));

    if(Collect_Input(ui->SeqConvergence_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_Convergence_Criteria))
        xmlWriter.writeTextElement("Sequence_Convergence_Criteria",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_Convergence_Criteria));

    if(Collect_Input(ui->SeqRecenter_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_Recenter_Iterations))
        xmlWriter.writeTextElement("Sequence_Recenter_Iterations",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_Recenter_Iterations));

    if(Collect_Input(ui->SeqNewCenter_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_New_Center_Iterations))
        xmlWriter.writeTextElement("Sequence_New_Center_Iterations",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_New_Center_Iterations));

    if(Collect_Input(ui->SeqMinCenters_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_Min_Centers))
        xmlWriter.writeTextElement("Sequence_Min_Centers",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_Min_Centers));

    if(Collect_Input(ui->SeqMaxCenters_Input,Dummy_Sim_Parms.Cycle_Params.Sequence_Max_Centers))
        xmlWriter.writeTextElement("Sequence_Max_Centers",QString::number(Dummy_Sim_Parms.Cycle_Params.Sequence_Max_Centers));

    if(Collect_Input(ui->CycleSwitch_Input,Dummy_Sim_Parms.Cycle_Params.Max_Cycle_Switch_Threshold))
        xmlWriter.writeTextElement("Max_Cycle_Switch_Threshold",QString::number(Dummy_Sim_Parms.Cycle_Params.Max_Cycle_Switch_Threshold));

    if(Collect_Input(ui->CommonPatternFreq_Input,Dummy_Sim_Parms.Cycle_Params.Min_Prop_for_Common_Pattern))
        xmlWriter.writeTextElement("Min_Prop_for_Common_Pattern",QString::number(Dummy_Sim_Parms.Cycle_Params.Min_Prop_for_Common_Pattern));

    if(Collect_Input(ui->MinPatternObs_Input,Dummy_Sim_Parms.Cycle_Params.Min_Occs_for_Common_Pattern))
        xmlWriter.writeTextElement("Min_Occs_for_Common_Pattern",QString::number(Dummy_Sim_Parms.Cycle_Params.Min_Occs_for_Common_Pattern));

    if(Collect_Input(ui->RoundTimesTo_Input,Dummy_Sim_Parms.Cycle_Params.Round_Times_To))
        xmlWriter.writeTextElement("Round_Times_To",QString::number(Dummy_Sim_Parms.Cycle_Params.Round_Times_To));

    //

    QString Dummy_Ext_Thresh;

    if(Collect_Input(ui->Extinction_Input,Dummy_Ext_Thresh))
        xmlWriter.writeTextElement("Simulation_Extinction_Threshold",Dummy_Ext_Thresh);

    if(Collect_Input(ui->nFood_Webs_Input,Dummy_Sim_Parms.nFood_Webs))
        xmlWriter.writeTextElement("N_Unique_Food_Webs",QString::number(Dummy_Sim_Parms.nFood_Webs));

    if(Collect_Input(ui->nSpatial_Structures_Input,Dummy_Sim_Parms.nSpatial_Structures))
        xmlWriter.writeTextElement("N_Unique_Spatial_Structures",QString::number(Dummy_Sim_Parms.nSpatial_Structures));

    xmlWriter.writeTextElement("Simulation_Species_Data_Write",QString::number(ui->Write_Species_Data_Input->currentIndex()));

    xmlWriter.writeTextElement("Simulation_Time_Series_Write",QString::number(ui->Write_Combo->currentIndex()));

    xmlWriter.writeTextElement("Simulation_Limit_Cycle_Write",QString::number(ui->Limit_Write_Combo->currentIndex()));

    xmlWriter.writeTextElement("Simulation_Write_Each_Set", QString::number(ui->WritePerSet_Input->isChecked()));

    xmlWriter.writeTextElement("Test_Against_Null_Structure", QString::number(ui->TestNullStructure_Input->isChecked()));

    xmlWriter.writeTextElement("Write_Measure_Diagnostics", QString::number(ui->WriteDiagnostic_Input->isChecked()));

    }
    xmlWriter.writeEndElement(); //Simulation_Properties

    xmlWriter.writeStartElement("Food_Web");
    xmlWriter.writeAttribute("Method",ui->Food_Web_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Food_Web_Method->currentIndex()));

    switch(ui->Food_Web_Method->currentIndex())
    {
    case 0:
    {
        xmlWriter.writeStartElement("Fixed_Structures");

        for(int i=0;i<Multi_Food_Web.size();i++)
        {
            xmlWriter.writeStartElement("structure");
            for(int j=0;j<Multi_Food_Web[i].shape()[0];j++)
            {
                xmlWriter.writeStartElement("row");
                for(int k=0;k<Multi_Food_Web[i].shape()[1];k++)
                    xmlWriter.writeTextElement("col",QString::number(Multi_Food_Web[i][j][k]));

                 xmlWriter.writeEndElement(); //row
            }
            xmlWriter.writeEndElement(); //structure
        }

        xmlWriter.writeEndElement(); //Fixed_Structures
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->Connectance_Input,Parameters))
            xmlWriter.writeTextElement("Connectance",Parameters);
    }
        break;
    }
    xmlWriter.writeEndElement(); //Food_Web

    xmlWriter.writeStartElement("Spatial_Structure");
    xmlWriter.writeAttribute("Method",ui->Spatial_Structure_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Spatial_Structure_Method->currentIndex()));

    switch(ui->Spatial_Structure_Method->currentIndex())
    {
    case 0:
    {
        xmlWriter.writeStartElement("Fixed_Structures");

        for(int i=0;i<Multi_Spatial_Structure.size();i++)
        {
            xmlWriter.writeStartElement("structure");
            for(int j=0;j<Multi_Spatial_Structure[i].shape()[0];j++)
            {
                xmlWriter.writeStartElement("row");
                for(int k=0;k<Multi_Spatial_Structure[i].shape()[1];k++)
                    xmlWriter.writeTextElement("col",QString::number(Multi_Spatial_Structure[i][j][k]));

                 xmlWriter.writeEndElement(); //row
            }
            xmlWriter.writeEndElement(); //structure
        }

        xmlWriter.writeEndElement(); //Fixed_Structures
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->Connectance_Spatial_Input,Parameters))
            xmlWriter.writeTextElement("Connectance",Parameters);
    }
        break;
    case 2:
    {
        QString Parameters;

        if(Collect_Input(ui->Branching_Input,Parameters))
            xmlWriter.writeTextElement("Branching_Probability",Parameters);
    }
        break;
    case 3:
    {
        QString Parameters;

        if(Collect_Input(ui->nBranches_Input,Parameters))
            xmlWriter.writeTextElement("Radial_Tree_nBranches",Parameters);

        if(Collect_Input(ui->Neighbor_Distance_Input,Parameters))
            xmlWriter.writeTextElement("Ring_Lattice_Neighbor_Distance",Parameters);
    }
        break;
    case 4:
    {
        QString Parameters;

        if(Collect_Input(ui->SS_N_Rewire_Input,Parameters))
            xmlWriter.writeTextElement("nRewire",Parameters);
    }
        break;
    }
    xmlWriter.writeEndElement(); //Spatial_Structure

    xmlWriter.writeStartElement("Body_Mass");
    xmlWriter.writeAttribute("Method",ui->Body_Mass_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Body_Mass_Method->currentIndex()));

    switch(ui->Body_Mass_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> Parameters;

        if(Collect_Input_Table_Column(0,ui->Body_Mass_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_Body_Mass");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_Body_Mass
        }
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->Base_mass_Input,Parameters))
            xmlWriter.writeTextElement("Body_Mass_Base",Parameters);

        if(Collect_Input(ui->Mass_trophic_Input,Parameters))
            xmlWriter.writeTextElement("Body_Mass_Trophic_Scaling",Parameters);

        if(Collect_Input(ui->Mass_variation_Input,Parameters))
            xmlWriter.writeTextElement("Body_Mass_Variation",Parameters);

        xmlWriter.writeTextElement("Body_Mass_Max_Trophic",QString::number(ui->Body_Mass_Scaling_Input->currentIndex()));

    }
        break;
    }

    xmlWriter.writeEndElement(); //Body_Mass

    xmlWriter.writeStartElement("Metabolism");
    xmlWriter.writeAttribute("Method",ui->Metabolism_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Metabolism_Method->currentIndex()));

    switch(ui->Metabolism_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> Parameters;

        if(Collect_Input_Table_Column(0,ui->Metabolic_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_ar");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_ar
        }

        if(Collect_Input_Table_Column(1,ui->Metabolic_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_ax");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_ax
        }

        if(Collect_Input_Table_Column(2,ui->Metabolic_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_y");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_y
        }

        if(Collect_Input_Table_Column(3,ui->Metabolic_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_e");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_e
        }
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->ar_Input,Parameters))
            xmlWriter.writeTextElement("ar",Parameters);

        if(Collect_Input(ui->ax_Input,Parameters))
            xmlWriter.writeTextElement("ax",Parameters);

        if(Collect_Input(ui->y_Input,Parameters))
            xmlWriter.writeTextElement("y",Parameters);

        if(Collect_Input(ui->e_Input,Parameters))
            xmlWriter.writeTextElement("e",Parameters);
    }
        break;
    }

    xmlWriter.writeEndElement(); //Metabolism

    xmlWriter.writeStartElement("Functional_Response");
    xmlWriter.writeAttribute("Method",ui->Functional_Response_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Functional_Response_Method->currentIndex()));

    switch(ui->Functional_Response_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> Parameters;

        if(Collect_Input_Table_Column(0,ui->Functional_Response_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_B0");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_B0
        }

        if(Collect_Input_Table_Column(1,ui->Functional_Response_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_h");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_h
        }
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->B0_Input,Parameters))
            xmlWriter.writeTextElement("B0",Parameters);

        if(Collect_Input(ui->h_Input,Parameters))
            xmlWriter.writeTextElement("h",Parameters);
    }
        break;
    }

    xmlWriter.writeEndElement(); //Functional_Response

    xmlWriter.writeStartElement("Dispersal");
    xmlWriter.writeAttribute("Method",ui->Dispersal_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Dispersal_Method->currentIndex()));

    switch(ui->Dispersal_Method->currentIndex())
    {
    case 0:
    {
        std::vector<double> Parameters;

        if(Collect_Input_Table_Column(0,ui->Dispersal_Input,Parameters))
        {
            xmlWriter.writeStartElement("Fixed_Dispersal");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Fixed_Dispersal
        }
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->m_Input,Parameters))
            xmlWriter.writeTextElement("Dispersal_Base",Parameters);

        if(Collect_Input(ui->m_Trophic_Input,Parameters))
            xmlWriter.writeTextElement("Dispersal_Trophic_Scaling",Parameters);

        if(Collect_Input(ui->m_Competitive_Input,Parameters))
            xmlWriter.writeTextElement("Dispersal_Competitive_Scaling",Parameters);

        if(Collect_Input(ui->m_Random_Input,Parameters))
            xmlWriter.writeTextElement("Dispersal_Random_Scaling",Parameters);
    }
        break;
    case 2:
    {
        QString Parameters;

        if(Collect_Input(ui->m_Input_2,Parameters))
            xmlWriter.writeTextElement("Dispersal_Base",Parameters);

        if(Collect_Input(ui->m_Trophic_Input_2,Parameters))
            xmlWriter.writeTextElement("Dispersal_Trophic_Scaling",Parameters);

        if(Collect_Input(ui->m_Competitive_Input_2,Parameters))
            xmlWriter.writeTextElement("Dispersal_Competitive_Scaling",Parameters);

        if(Collect_Input(ui->m_Random_Input_2,Parameters))
            xmlWriter.writeTextElement("Dispersal_Random_Scaling",Parameters);
    }
        break;
    }
    xmlWriter.writeEndElement(); //Dispersal

    xmlWriter.writeStartElement("Carrying_Capacity");
    xmlWriter.writeAttribute("Method",ui->Carrying_Capacity_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Carrying_Capacity_Method->currentIndex()));

    switch(ui->Carrying_Capacity_Method->currentIndex())
    {
    case 0:
    {
        boost_matrix Structure;

        if(Collect_Input_Table(ui->Carrying_Capacity_Input,Structure))
        {

            xmlWriter.writeStartElement("Fixed_Carrying_Capacity");
            for(int i=0;i<ui->Species_Input->value();i++)
            {
                xmlWriter.writeStartElement("row");
                for(int i2=0;i2<ui->Patches_Input->value();i2++)
                    xmlWriter.writeTextElement("col",QString::number(Structure[i][i2]));

                xmlWriter.writeEndElement(); //row
            }

            xmlWriter.writeEndElement(); //Fixed_Carrying_Capacity
        }

        boost_matrix Patch_Size;

        if(Collect_Input_Table(ui->Patch_Size_Input,Patch_Size))
        {
            xmlWriter.writeStartElement("Fixed_Patch_Size");
            for(int i=0;i<ui->Patches_Input->value();i++)
            {
                xmlWriter.writeTextElement("col",QString::number(Structure[0][i]));
            }

            xmlWriter.writeEndElement(); //Fixed_Patch_Size
        }
    }
        break;
    case 1:
    {
        QString Parameters;

        if(Collect_Input(ui->Carrying_Capacity_Fixed_Community_Input,Parameters))
            xmlWriter.writeTextElement("Carrying_Capacity",Parameters);
    }
        break;
    case 2:
    {
        boost_matrix Structure;

        if(Collect_Input_Table(ui->Patch_Size_Input_2,Structure))
        {
            xmlWriter.writeStartElement("Patch_Size");
            for(int i=0;i<Structure.shape()[0];i++)
            {
                xmlWriter.writeStartElement("row");
                for(int i2=0;i2<Structure.shape()[1];i2++)
                    xmlWriter.writeTextElement("col",QString::number(Structure[i][i2]));

                xmlWriter.writeEndElement(); //row
            }

            xmlWriter.writeEndElement(); //Patch_Size
        }
    }
        break;
    }

    xmlWriter.writeEndElement(); //Carrying_Capacity

    xmlWriter.writeStartElement("Initial_Abundance");
    xmlWriter.writeAttribute("Method",ui->Initial_Abundance_Method->currentText());
    xmlWriter.writeAttribute("Method_Index",QString::number(ui->Initial_Abundance_Method->currentIndex()));

    switch(ui->Initial_Abundance_Method->currentIndex())
    {
    case 0:
    {

        boost_matrix Structure;

        if(Collect_Input_Table(ui->Initial_Abundance_Input,Structure))
        {
            xmlWriter.writeStartElement("Fixed_Initial_Abundance");
            for(int i=0;i<ui->Species_Input->value();i++)
            {
                xmlWriter.writeStartElement("row");
                for(int i2=0;i2<ui->Patches_Input->value();i2++)
                    xmlWriter.writeTextElement("col",QString::number(Structure[i][i2]));

                xmlWriter.writeEndElement(); //row
            }

            xmlWriter.writeEndElement(); //Fixed_Initial_Abundance
        }
    }
        break;
    case 1:
    {
        std::vector<double> Parameters;

        if(Collect_Input_Table_Column(0,ui->Initial_Abundance_Random_Mean_Input,Parameters))
        {
            xmlWriter.writeStartElement("Random_Initial_Abundance_Mean");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Random_Initial_Abundance_Mean
        }

        if(Collect_Input_Table_Column(0,ui->Initial_Abundance_Random_Dev_Input,Parameters))
        {
            xmlWriter.writeStartElement("Random_Initial_Abundance_Dev");
            for(int i=0;i<ui->Species_Input->value();i++)
                xmlWriter.writeTextElement("row",QString::number(Parameters[i]));

            xmlWriter.writeEndElement(); //Random_Initial_Abundance_Dev
        }
    }
        break;
    }
    xmlWriter.writeEndElement(); //Initial_Abundance

    xmlWriter.writeEndElement(); //Allometric_Metacommunity

    xmlWriter.writeEndDocument();

    Write_File.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Read XML Values

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::on_actionLoad_triggered()
{
    QString Read_File = QFileDialog::getOpenFileName(this,
                                           tr("Import Xml"), ".",
                                           tr("Xml files (*.xml)"));

    QFile Defaults_File(Read_File);

    if (!Defaults_File.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QMessageBox::critical(this,"Load XML File Problem",
                              "Couldn't load from file",
                              QMessageBox::Ok);
        return;
    }

    QXmlStreamReader xmlReader(Defaults_File.readAll());

    Defaults_File.close();

    while(!xmlReader.atEnd() && !xmlReader.hasError())
    {
        xmlReader.readNext();

        if(xmlReader.isStartDocument())
            continue;

        if(xmlReader.isStartElement())
        {
            if(xmlReader.name()=="Allometric_Metacommunity")
                continue;

            if(xmlReader.name()=="nSpecies")
            {
                ui->Species_Input->setValue(xmlReader.readElementText().toInt());
                continue;
            }
            if(xmlReader.name()=="nPatches")
            {
                ui->Patches_Input->setValue(xmlReader.readElementText().toInt());
                continue;
            }
            if(xmlReader.name()=="Simulation_Properties")
            {
                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Simulation_Properties"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Max_Simulation_Length")
                    {
                        ui->MaxSimLength_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Min_Simulation_Length")
                    {
                        ui->MinSimLength_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Test_Interval_Length")
                    {
                        ui->TestInt_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Length_Between_Test_Intervals")
                    {
                        ui->LenBetweenTest_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Iterations")
                    {
                        ui->nIterations_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Resolution")
                    {
                        ui->SimRes_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Threads")
                    {
                        ui->nThreads_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Extinction_Threshold")
                    {
                        ui->Extinction_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="N_Unique_Food_Webs")
                    {
                        ui->nFood_Webs_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="N_Unique_Spatial_Structures")
                    {
                        ui->nSpatial_Structures_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Time_Series_Write")
                    {
                        ui->Write_Combo->setCurrentIndex(xmlReader.readElementText().toInt());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Limit_Cycle_Write")
                    {
                        ui->Limit_Write_Combo->setCurrentIndex(xmlReader.readElementText().toInt());
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Species_Data_Write")
                    {
                        ui->Write_Species_Data_Input->setCurrentIndex(xmlReader.readElementText().toInt());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Phase_Lock_Threshold")
                    {
                        ui->PhaseLockFreq_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    //

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Stationary_Cycle_Frequency")
                    {
                        ui->StationaryCycleFreq_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Proportion_Stationary_Cycles")
                    {
                        ui->StationaryCycleProp_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Difference_Error_Tolerance")
                    {
                        ui->DiffErrTol_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Minimum_Extrema_Proportion_Difference")
                    {
                        ui->MinExtremaPropDiff_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_Center_Target_CV")
                    {
                        ui->SeqCenterCV_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_Threshold_to_Increase_Centers")
                    {
                        ui->SeqThreshIncCenters_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_Convergence_Criteria")
                    {
                        ui->SeqConvergence_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_Recenter_Iterations")
                    {
                        ui->SeqRecenter_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_New_Center_Iterations")
                    {
                        ui->SeqNewCenter_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_Min_Centers")
                    {
                        ui->SeqMinCenters_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Sequence_Max_Centers")
                    {
                        ui->SeqMaxCenters_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Max_Cycle_Switch_Threshold")
                    {
                        ui->CycleSwitch_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Min_Occs_for_Common_Pattern")
                    {
                        ui->MinPatternObs_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Min_Prop_for_Common_Pattern")
                    {
                        ui->CommonPatternFreq_Input->setText(xmlReader.readElementText());
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Round_Times_To")
                    {
                        ui->RoundTimesTo_Input->setText(xmlReader.readElementText());
                        continue;
                    }


                    //

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Simulation_Write_Each_Set")
                    {
                        bool dummy_check;

                        int input = xmlReader.readElementText().toInt();

                        if(input==1)
                            dummy_check = true;
                        else
                            if(input==0)
                                dummy_check=false;
                            else
                                continue;

                        ui->WritePerSet_Input->setChecked(dummy_check);
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Test_Against_Null_Structure")
                    {
                        bool dummy_check;

                        int input = xmlReader.readElementText().toInt();

                        if(input==1)
                            dummy_check = true;
                        else
                            if(input==0)
                                dummy_check=false;
                            else
                                continue;

                        ui->TestNullStructure_Input->setChecked(dummy_check);
                        continue;
                    }

                    if(xmlReader.isStartElement()&&xmlReader.name()=="Write_Measure_Diagnostics")
                    {
                        bool dummy_check;

                        int input = xmlReader.readElementText().toInt();

                        if(input==1)
                            dummy_check = true;
                        else
                            if(input==0)
                                dummy_check=false;
                            else
                                continue;

                        ui->WriteDiagnostic_Input->setChecked(dummy_check);
                        continue;
                    }
                }
                continue;
            }
            if(xmlReader.name()=="Food_Web")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Food_Web_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Food_Web"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_Structures")
                            {
                                int i=0;
                                while(!xmlReader.hasError())
                                {
                                    xmlReader.readNext();
                                    if(xmlReader.isStartElement()&&xmlReader.name()=="structure")
                                    {
                                        if(Multi_Food_Web.size()>(i+1))
                                            Multi_Food_Web.resize(i+1);
                                        else if(Multi_Food_Web.size()==i)
                                            Multi_Food_Web.push_back(boost_matrix());
                                        else if(Multi_Food_Web.size()<i)
                                        {
                                            while(Multi_Food_Web.size()<i)
                                            {
                                                Multi_Food_Web.push_back(boost_matrix(boost::extents[0][0]));
                                                Multi_Food_Web.back()[0][0] = 0;
                                            }
                                        }

                                        Read_XML_Matrix(&xmlReader,Multi_Food_Web[i]);

                                        i++;
                                    }
                                    else if(xmlReader.isEndElement()&&!(xmlReader.name()=="structure"))
                                            break;
                                }

                                ui->Multi_Food_Web_Spin->setValue(1);
                                on_Multi_Food_Web_Spin_valueChanged(1);
                                ui->Multi_Food_Web_Spin->setMaximum(Multi_Food_Web.size());

                                continue;
                            }
                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="Connectance")
                            {
                                ui->Connectance_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        }
                    }

                }
                continue;
            }
            if(xmlReader.name()=="Spatial_Structure")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Spatial_Structure_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Spatial_Structure"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_Structures")
                            {   
                                int i=0;
                                while(!xmlReader.hasError())
                                {
                                    xmlReader.readNext();
                                    if(xmlReader.isStartElement()&&xmlReader.name()=="structure")
                                    {
                                        if(Multi_Spatial_Structure.size()>(i+1))
                                            Multi_Spatial_Structure.resize(i+1);
                                        else if(Multi_Spatial_Structure.size()==i)
                                            Multi_Spatial_Structure.push_back(boost_matrix());
                                        else if(Multi_Spatial_Structure.size()<i)
                                        {
                                            while(Multi_Spatial_Structure.size()<i)
                                            {
                                                Multi_Spatial_Structure.push_back(boost_matrix(boost::extents[0][0]));
                                                Multi_Spatial_Structure.back()[0][0] = 0;
                                            }
                                        }

                                        Read_XML_Matrix(&xmlReader,Multi_Spatial_Structure[i]);

                                        i++;
                                    }
                                    else if(xmlReader.isEndElement()&&!(xmlReader.name()=="structure"))
                                            break;
                                }

                                ui->Multi_Spatial_Structure_Spin->setValue(1);
                                on_Multi_Spatial_Structure_Spin_valueChanged(1);
                                ui->Multi_Spatial_Structure_Spin->setMaximum(Multi_Spatial_Structure.size());

                                continue;
                            }

                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="Connectance")
                            {
                                ui->Connectance_Spatial_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        case 2:
                        {
                            if(xmlReader.name()=="Branching_Probability")
                            {
                                ui->Branching_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        case 3:
                        {
                            if(xmlReader.name()=="Radial_Tree_nBranches")
                            {
                                ui->nBranches_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                            if(xmlReader.name()=="Ring_Lattice_Neighbor_Distance")
                            {
                                ui->Neighbor_Distance_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        case 4:
                        {
                            if(xmlReader.name()=="nRewire")
                            {
                                ui->SS_N_Rewire_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                        }
                            break;
                        }
                    }

                }
                continue;
            }

            if(xmlReader.name()=="Body_Mass")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Body_Mass_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Body_Mass"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_Body_Mass")
                            {
                                Read_XML_Column(&xmlReader,ui->Body_Mass_Input,0);

                                continue;
                            }

                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="Body_Mass_Base")
                            {
                                ui->Base_mass_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                            if(xmlReader.name()=="Body_Mass_Trophic_Scaling")
                            {
                                ui->Mass_trophic_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                            if(xmlReader.name()=="Body_Mass_Variation")
                            {
                                ui->Mass_variation_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Body_Mass_Max_Trophic")
                            {
                                ui->Body_Mass_Scaling_Input->setCurrentIndex(xmlReader.readElementText().toInt());

                                continue;
                            }


                        }
                            break;
                        }
                    }

                }
                continue;
            }
            if(xmlReader.name()=="Metabolism")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Metabolism_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Metabolism"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_ar")
                            {
                                Read_XML_Column(&xmlReader,ui->Metabolic_Input,0);

                                continue;
                            }

                            if(xmlReader.name()=="Fixed_ax")
                            {
                                Read_XML_Column(&xmlReader,ui->Metabolic_Input,1);

                                continue;
                            }

                            if(xmlReader.name()=="Fixed_y")
                            {
                                Read_XML_Column(&xmlReader,ui->Metabolic_Input,2);

                                continue;
                            }

                            if(xmlReader.name()=="Fixed_e")
                            {
                                Read_XML_Column(&xmlReader,ui->Metabolic_Input,3);

                                continue;
                            }

                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="ar")
                            {
                                ui->ar_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="ax")
                            {
                                ui->ax_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="y")
                            {
                                ui->y_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="e")
                            {
                                ui->e_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        }
                    }

                }
                continue;
            }
            if(xmlReader.name()=="Functional_Response")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Functional_Response_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Functional_Response"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_B0")
                            {
                                Read_XML_Column(&xmlReader,ui->Functional_Response_Input,0);

                                continue;
                            }

                            if(xmlReader.name()=="Fixed_h")
                            {
                                Read_XML_Column(&xmlReader,ui->Functional_Response_Input,1);

                                continue;
                            }
                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="B0")
                            {
                                ui->B0_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="h")
                            {
                                ui->h_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        }
                    }

                }
                continue;
            }

            if(xmlReader.name()=="Dispersal")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Dispersal_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Dispersal"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_Dispersal")
                            {
                                Read_XML_Column(&xmlReader,ui->Dispersal_Input,0);

                                continue;
                            }
                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="Dispersal_Base")
                            {
                                ui->m_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Dispersal_Trophic_Scaling")
                            {
                                ui->m_Trophic_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Dispersal_Competitive_Scaling")
                            {
                                ui->m_Competitive_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Dispersal_Random_Scaling")
                            {
                                ui->m_Random_Input->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        case 2:
                        {
                            if(xmlReader.name()=="Dispersal_Base")
                            {
                                ui->m_Input_2->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Dispersal_Trophic_Scaling")
                            {
                                ui->m_Trophic_Input_2->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Dispersal_Competitive_Scaling")
                            {
                                ui->m_Competitive_Input_2->setText(xmlReader.readElementText());

                                continue;
                            }
                            if(xmlReader.name()=="Dispersal_Random_Scaling")
                            {
                                ui->m_Random_Input_2->setText(xmlReader.readElementText());

                                continue;
                            }

                        }
                            break;
                        }
                    }

                }
                continue;
            }

            if(xmlReader.name()=="Carrying_Capacity")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Carrying_Capacity_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Carrying_Capacity"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_Carrying_Capacity")
                            {
                                Read_XML_Matrix(&xmlReader,ui->Carrying_Capacity_Input);

                                continue;
                            }

                            if(xmlReader.name()=="Fixed_Patch_Size")
                            {
                                Read_XML_Matrix(&xmlReader,ui->Patch_Size_Input);
                            }

                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="Carrying_Capacity")
                            {
                                ui->Carrying_Capacity_Fixed_Community_Input->setText(xmlReader.readElementText());

                                continue;
                            }
                        }
                            break;
                        case 2:
                        {
                            if(xmlReader.name()=="Patch_Size")
                            {
                                Read_XML_Matrix(&xmlReader,ui->Patch_Size_Input_2);
                            }
                        }
                            break;
                        }
                    }

                }
                continue;
            }

            if(xmlReader.name()=="Initial_Abundance")
            {
                int Method_Index = 0;

                foreach(const QXmlStreamAttribute &attr, xmlReader.attributes())
                {
                    if(attr.name().toString() == QLatin1String("Method_Index"))
                        Method_Index = attr.value().toInt();
                }

                ui->Initial_Abundance_Method->setCurrentIndex(Method_Index);

                while(!xmlReader.hasError()&&!(xmlReader.isEndElement()&&xmlReader.name()=="Initial_Abundance"))
                {
                    xmlReader.readNext();

                    if(xmlReader.isStartElement())
                    {
                        switch(Method_Index)
                        {
                        case 0:
                        {
                            if(xmlReader.name()=="Fixed_Initial_Abundance")
                            {
                                Read_XML_Matrix(&xmlReader,ui->Initial_Abundance_Input);

                                continue;
                            }
                        }
                            break;
                        case 1:
                        {
                            if(xmlReader.name()=="Random_Initial_Abundance_Mean")
                            {
                                Read_XML_Column(&xmlReader,ui->Initial_Abundance_Random_Mean_Input,0);

                                continue;
                            }

                            if(xmlReader.name()=="Random_Initial_Abundance_Dev")
                            {
                                Read_XML_Column(&xmlReader,ui->Initial_Abundance_Random_Dev_Input,0);

                                continue;
                            }
                        }
                            break;
                        }
                    }

                }
                continue;
            }


        }

    }

    if(xmlReader.hasError())
        QMessageBox::critical(this,"Load XML File Problem",
                          xmlReader.errorString(),
                          QMessageBox::Ok);
}












