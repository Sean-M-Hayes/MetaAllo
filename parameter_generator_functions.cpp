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
#include "parameter_generator_functions.h"

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

//General

std::vector<std::string> Parameter_Generators::Get_Generator_Set_Method_Labels() const
{
    std::vector<std::string> to_Return;

    std::string Labels[] = {"H","AMI","Hc"};

    to_Return.assign(Labels,Labels+3);

    return to_Return;
}

std::vector<double> Parameter_Generators::Get_Generator_Set_Method_Data() const
{
    std::vector<double> to_Return;

    if(this->Check().empty())
    {
        boost_matrix Complete_Food_Web(this->Food_Web->Get_Value().Food_Web);

        for(int i1 =0;i1<Complete_Food_Web.shape()[0];i1++){
            for(int i2=0;i2<Complete_Food_Web.shape()[1];i2++){

            Complete_Food_Web[i2][i1] = this->Metabolism->Get_ax()[i1]*pow(this->Body_Mass->Get_Value()[i1],-.25)*this->Metabolism->Get_y()[i1]/this->Metabolism->Get_e()[i1]*Complete_Food_Web[i2][i1];

            }
        }

        Food_Web_Metadata Complete_Food_Web_Metadata(Complete_Food_Web);

        double Values[] = {Complete_Food_Web_Metadata.H,Complete_Food_Web_Metadata.AMI,Complete_Food_Web_Metadata.Hc};

        to_Return.assign(Values,Values+3);
    }

return to_Return;
}


std::string Parameter_Generators::Check() const
{
    std::string to_Return;

    if(Food_Web==NULL)
        to_Return+="Food Web, ";
    if(Spatial_Structure==NULL)
        to_Return+="Spatial Structure, ";
    if(Body_Mass==NULL)
        to_Return+="Body Mass, ";
    if(Metabolism==NULL)
        to_Return+="Metabolism, ";
    if(Functional_Response==NULL)
        to_Return+="Functional Response, ";
    if(Dispersal==NULL)
        to_Return+="Dispersal, ";
    if(Carrying_Capacity==NULL)
        to_Return+="Carrying Capacity, ";
    if(Initial_Abundance==NULL)
        to_Return+="Initial Abundance, ";
    if(Extinction_Threshold==NULL)
        to_Return+="Extinction_Threshold, ";

    if(!to_Return.empty())
    {
        to_Return.erase(to_Return.end()-2,to_Return.end());
        to_Return+=" parameters are missing.";
    }

    return to_Return;
}

//order matters in current iteration - values generated in order and some have dependencies (Body_Mass/Dispersal depend on Food_Web sometimes)
Parameter_Generator* Parameter_Generators::operator[](int i) const
{
    switch(i)
    {
    case 0:
    {
        return Food_Web;
    }
        break;
    case 1:
    {
        return Spatial_Structure;
    }
        break;
    case 2:
    {
        return Body_Mass;
    }
        break;
    case 3:
    {
        return Metabolism;
    }
        break;
    case 4:
    {
        return Functional_Response;
    }
        break;
    case 5:
    {
        return Dispersal;
    }
        break;
    case 6:
    {
        return Carrying_Capacity;
    }
        break;
    case 7:
    {
        return Initial_Abundance;
    }
        break;
    case 8:
    {
        return Extinction_Threshold;
    }
    }
}

int Parameter_Generators::nGenerators() const
{
    return 9;
}

Parameter_Generators::Parameter_Generators()
{
    Food_Web = NULL;

    Spatial_Structure = NULL;

    Body_Mass = NULL;

    Metabolism = NULL;

    Functional_Response = NULL;

    Dispersal = NULL;

    Carrying_Capacity = NULL;

    Initial_Abundance = NULL;

    Extinction_Threshold = NULL;

}

Parameter_Generators::Parameter_Generators(const Parameter_Generators &obj)
{
    Food_Web = obj.Food_Web->Clone();
    Food_Web->Set_Container(this);

    Spatial_Structure = obj.Spatial_Structure->Clone();
    Spatial_Structure->Set_Container(this);

    Body_Mass = obj.Body_Mass->Clone();
    Body_Mass->Set_Container(this);

    Metabolism = obj.Metabolism->Clone();
    Metabolism->Set_Container(this);

    Functional_Response = obj.Functional_Response->Clone();
    Functional_Response->Set_Container(this);

    Dispersal = obj.Dispersal->Clone();
    Dispersal->Set_Container(this);

    Carrying_Capacity = obj.Carrying_Capacity->Clone();
    Carrying_Capacity->Set_Container(this);

    Initial_Abundance = obj.Initial_Abundance->Clone();
    Initial_Abundance->Set_Container(this);

    Extinction_Threshold = obj.Extinction_Threshold->Clone();
    Extinction_Threshold->Set_Container(this);
}

Parameter_Generators::~Parameter_Generators()
{
    delete Food_Web;
    delete Spatial_Structure;
    delete Body_Mass;
    delete Metabolism;
    delete Functional_Response;
    delete Dispersal;
    delete Carrying_Capacity;
    delete Initial_Abundance;
    delete Extinction_Threshold;
}


Parameter_Generator::Parameter_Generator(Parameter_Generators* Container,int In_nSets):
    Container_Pointer(Container),nSets(In_nSets),iSet(0){}

int Parameter_Generator::Current_Set() const
{
    return iSet;
}

int Parameter_Generator::Get_nSets() const
{
    return nSets;
}

std::vector<std::string> Parameter_Generator::Get_Generator_Method_Labels() const
{
    return Generator_Data_Labels;
}

std::vector<double> Parameter_Generator::Get_Generator_Method_Data() const
{
    return Generator_Data;
}

bool Parameter_Generator::Increment_Set()
{
    if(iSet<(nSets-1))
    {
        iSet++;
        Update_Generator_Methods();
        return true;
    }
    else
    {
        return false;
    }
}

void Parameter_Generator::Reset()
{
    iSet=0;
    Update_Generator_Methods();
}

void Parameter_Generator::Set_Container(Parameter_Generators* New_Container)
{
    Container_Pointer = New_Container;
}

//Food_Webs

Food_Web_Generator::Food_Web_Generator(Parameter_Generators* Container,boost_matrix Input):
    Parameter_Generator(Container,1),Food_Web_Data(Input)
{
    Update_Generator_Methods();
}

Food_Web_Metadata Food_Web_Generator::Get_Value()
{
    return Food_Web_Data;
}

void Food_Web_Generator::Set_Value(const Food_Web_Metadata &Input)
{
    Food_Web_Data(Input);
}

void Food_Web_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Food_Web");
    Generator_Data_Labels.push_back("nSpecies");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Food_Web_Data.Food_Web.shape()[0]);
}

Food_Web_Generator* Food_Web_Generator::Clone()
{
    return new Food_Web_Generator(*this);
}

//Food Web Random Gen

Food_Web_Niche_Model_Generator::Food_Web_Niche_Model_Generator(Parameter_Generators* Container,int Input_nSpecies, std::vector<double> Input_Connectance):
    Food_Web_Generator(Container,Input_Connectance.size()),nSpecies(Input_nSpecies), Connectance(Input_Connectance)
{    
    Generate_New();

    Update_Generator_Methods();
}


void Food_Web_Niche_Model_Generator::Generate_New()
{
    Food_Web_Data(Niche_Model(nSpecies, Connectance[iSet]));
}

void Food_Web_Niche_Model_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Niche_Model_Food_Web");
    Generator_Data_Labels.push_back("nSpecies");
    Generator_Data_Labels.push_back("Connectance");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(nSpecies);
    Generator_Data.push_back(Connectance[iSet]);
}

Food_Web_Niche_Model_Generator* Food_Web_Niche_Model_Generator::Clone()
{
    return new Food_Web_Niche_Model_Generator(*this);
}

//Food Web Derived Gen - Read

Food_Web_Read_Generator::Food_Web_Read_Generator(Parameter_Generators* Container,std::vector<boost_matrix> Input_Food_Webs):
    Food_Web_Generator(Container,1),Read_Food_Webs(Input_Food_Webs),Interset_Iter(0)
{
    Generate_New();

    Interset_Iter = 0; //So the first value it generates will be the first value read in, and is also primed with first value read in

    Update_Generator_Methods();
}

void Food_Web_Read_Generator::Generate_New()
{
    if(Interset_Iter>=Read_Food_Webs.size())
    {
        Interset_Iter=0;
    }

    Food_Web_Data(Read_Food_Webs[Interset_Iter]);

    Interset_Iter++;
}

void Food_Web_Read_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Read_Food_Web");

    Generator_Data.push_back(Interset_Iter);
}

Food_Web_Read_Generator* Food_Web_Read_Generator::Clone()
{
    return new Food_Web_Read_Generator(*this);
}

//Spatial_Structure

Spatial_Structure_Generator::Spatial_Structure_Generator(Parameter_Generators* Container,boost_matrix Input):
    Parameter_Generator(Container,1), Spatial_Structure_Data(Calculate_Unweighted_Symmetrical_Spatial_Structure_Metadata(Input))
{    
    Update_Generator_Methods();
}

Spatial_Structure_Metadata Spatial_Structure_Generator::Get_Value()
{
    return Spatial_Structure_Data;
}

void Spatial_Structure_Generator::Set_Value(const Spatial_Structure_Metadata &Input)
{
    Spatial_Structure_Data(Input);
}

void Spatial_Structure_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Spatial_Structure");
    Generator_Data_Labels.push_back("nPatches");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Spatial_Structure_Data.Spatial_Structure.shape()[0]);
}

Spatial_Structure_Generator* Spatial_Structure_Generator::Clone()
{
    return new Spatial_Structure_Generator(*this);
}

//Spatial Structure Derived

Spatial_Structure_Undirected_Erdos_Renyi_Generator::Spatial_Structure_Undirected_Erdos_Renyi_Generator(Parameter_Generators* Container,int Input_nPatches,std::vector<double> Input_Connectance):
    Spatial_Structure_Generator(Container,Input_Connectance.size()), nPatches(Input_nPatches),Connectance(Input_Connectance)
{
    Generate_New();

    Update_Generator_Methods();
}
void Spatial_Structure_Undirected_Erdos_Renyi_Generator::Generate_New()
{
    Spatial_Structure_Data(Erdos_Renyi_Undirected(nPatches, Connectance[iSet]));
}

void Spatial_Structure_Undirected_Erdos_Renyi_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Undirected_Erdos_Renyi_Spatial_Structure");
    Generator_Data_Labels.push_back("nPatches");
    Generator_Data_Labels.push_back("Connectance");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(nPatches);
    Generator_Data.push_back(Connectance[iSet]);
}

Spatial_Structure_Undirected_Erdos_Renyi_Generator* Spatial_Structure_Undirected_Erdos_Renyi_Generator::Clone()
{
    return new Spatial_Structure_Undirected_Erdos_Renyi_Generator(*this);
}

//Spatial Structure Derived II - Dendritic

Spatial_Structure_Undirected_Dendritic_Generator::Spatial_Structure_Undirected_Dendritic_Generator(Parameter_Generators* Container,int Input_nPatches,std::vector<double> Input_Branching_Probability):
    Spatial_Structure_Generator(Container,Input_Branching_Probability.size()),nPatches(Input_nPatches),Branching_Probability(Input_Branching_Probability)
{
    Spatial_Structure_Data(Dendritic_Undirected(nPatches, Branching_Probability[0]));

    Update_Generator_Methods();
}
void Spatial_Structure_Undirected_Dendritic_Generator::Generate_New()
{
    Spatial_Structure_Data(Dendritic_Undirected(nPatches, Branching_Probability[iSet]));
}

void Spatial_Structure_Undirected_Dendritic_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Undirected_Dendritic_Spatial_Structure");
    Generator_Data_Labels.push_back("nPatches");
    Generator_Data_Labels.push_back("Connectance");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(nPatches);
    Generator_Data.push_back(Branching_Probability[iSet]);
}

Spatial_Structure_Undirected_Dendritic_Generator* Spatial_Structure_Undirected_Dendritic_Generator::Clone()
{
    return new Spatial_Structure_Undirected_Dendritic_Generator(*this);
}

//Spatial Structure Derived III - Kurt Comparison

Spatial_Structure_Kurt_Comparison_Generator::Spatial_Structure_Kurt_Comparison_Generator(Parameter_Generators* Container,int Input_nPatches,std::vector<int> Input_nBranches,std::vector<int> Input_Neighbor_Distance):
    Spatial_Structure_Generator(Container,Input_nBranches.size()+Input_Neighbor_Distance.size()),nPatches(Input_nPatches),nBranches(Input_nBranches), Neighbor_Distance(Input_Neighbor_Distance)
{
    Generate_New();

    Update_Generator_Methods();
}

void Spatial_Structure_Kurt_Comparison_Generator::Generate_New()
{
    if(iSet<nBranches.size())
        Spatial_Structure_Data(Radial_Tree_Undirected(nPatches,nBranches[iSet]));
    else
        Spatial_Structure_Data(Ring_Lattice_Undirected(nPatches,Neighbor_Distance[iSet-nBranches.size()]));
}

void Spatial_Structure_Kurt_Comparison_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Kurt_Comparison_Spatial_Structure");
    Generator_Data_Labels.push_back("nPatches");
    Generator_Data_Labels.push_back("nBranches or Neighbor_Distance");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(nPatches);
    if(iSet<nBranches.size())
        Generator_Data.push_back(nBranches[iSet]);
    else
        Generator_Data.push_back(Neighbor_Distance[iSet-nBranches.size()]);
}

Spatial_Structure_Kurt_Comparison_Generator* Spatial_Structure_Kurt_Comparison_Generator::Clone()
{
    return new Spatial_Structure_Kurt_Comparison_Generator(*this);
}

//Spatial Structure Derived IV - Shipley-Skinner

void Combinations(int offset, int k, const std::vector<int> source_vector, std::vector<int> combination, std::vector<std::vector<int>> &Output_Vector)
{
    if(k==0)
    {
        Output_Vector.push_back(combination);
        return;
    }
    for(int i=offset;i<=source_vector.size()-k;++i)
    {
        combination.push_back(source_vector[i]);
        Combinations(i+1,k-1,source_vector,combination,Output_Vector);
        combination.pop_back();
    }
}

Spatial_Structure_Shipley_Skinner_Generator::Spatial_Structure_Shipley_Skinner_Generator(Parameter_Generators* Container,std::vector<int> Input_nRewire):
    Spatial_Structure_Generator(Container,Input_nRewire.size()),nRewire(Input_nRewire),Interset_Iter(0),Previous_iSet(0)
{
    //15 Removeable edges

    //std::vector<std::vector<std::vector<double>>> Remove_sets;

    for(int i=0; i<nRewire.size();i++)
    {
        std::vector<std::vector<int>> to_Add;
        std::vector<int> blank_workspace;
        std::vector<int> source;
        for(int i=0;i<15;i++){source.push_back(i);}//i<15 for n removeable edges

        Combinations(0,nRewire[i],source,blank_workspace,to_Add);

        Remove_Sets.push_back(to_Add);
    }

    Generate_New();

    Update_Generator_Methods();
}

void Spatial_Structure_Shipley_Skinner_Generator::Generate_New()
{
    if(Previous_iSet!=iSet||Interset_Iter>=Remove_Sets[iSet].size())
    {
        Interset_Iter=0;
        Previous_iSet=iSet;
    }

    Spatial_Structure_Data(Permuted_Shipley_Skinner(Remove_Sets[iSet][Interset_Iter]));

    Interset_Iter++;
}

void Spatial_Structure_Shipley_Skinner_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Shipley_Skinner_Spatial_Structure");
    Generator_Data_Labels.push_back("nRewire");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(nRewire[iSet]);
}

Spatial_Structure_Shipley_Skinner_Generator* Spatial_Structure_Shipley_Skinner_Generator::Clone()
{
    return new Spatial_Structure_Shipley_Skinner_Generator(*this);
}

//Spatial Structure Derived V - Read

Spatial_Structure_Read_Generator::Spatial_Structure_Read_Generator(Parameter_Generators* Container,std::vector<boost_matrix> Input_Spatial_Structures):
    Spatial_Structure_Generator(Container,1),Read_Spatial_Structures(Input_Spatial_Structures),Interset_Iter(0)
{
    Generate_New();

    Interset_Iter = 0; //So the first value it generates will be the first value read in, and is also primed with first value read in

    Update_Generator_Methods();
}

void Spatial_Structure_Read_Generator::Generate_New()
{
    if(Interset_Iter>=Read_Spatial_Structures.size())
    {
        Interset_Iter=0;
    }

    Spatial_Structure_Data(Calculate_Unweighted_Symmetrical_Spatial_Structure_Metadata(Read_Spatial_Structures[Interset_Iter]));

    Interset_Iter++;
}

void Spatial_Structure_Read_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Read_Spatial_Structure");

    Generator_Data.push_back(Interset_Iter);
}

Spatial_Structure_Read_Generator* Spatial_Structure_Read_Generator::Clone()
{
    return new Spatial_Structure_Read_Generator(*this);
}


//Body Mass

Body_Mass_Generator::Body_Mass_Generator(Parameter_Generators* Container,std::vector<double> Input):
    Parameter_Generator(Container,1),Body_Mass(Input)
{
    Update_Generator_Methods();
}

std::vector<double> Body_Mass_Generator::Get_Value()
{
    return Body_Mass;
}

void Body_Mass_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Body_Mass");

    Generator_Data.push_back(iSet);
}

Body_Mass_Generator* Body_Mass_Generator::Clone()
{
    return new Body_Mass_Generator(*this);
}

//Changing Input_Max_Trophic_Position from bool to int
Body_Mass_Trophic_Scaling_Generator::Body_Mass_Trophic_Scaling_Generator(Parameter_Generators* Container, int Input_Body_Mass_Scaling_Method, std::vector<double> Input_Base_Mass, std::vector<double> Input_Trophic_Scaling, std::vector<double> Input_Mass_Variation):
    Body_Mass_Generator(Container,Input_Base_Mass.size()*Input_Trophic_Scaling.size()*Input_Mass_Variation.size()),Body_Mass_Scaling_Method(Input_Body_Mass_Scaling_Method)
{    
    Body_Mass_Set.resize(boost::extents[nSets][3]);

    int i_row = 0;

    for (int i=0; i<Input_Base_Mass.size();i++)
    {
        for(int i2=0; i2<Input_Trophic_Scaling.size();i2++)
        {
            for(int i3=0; i3<Input_Mass_Variation.size();i3++)
            {
                Body_Mass_Set[i_row][0] = Input_Base_Mass[i];
                Body_Mass_Set[i_row][1] = Input_Trophic_Scaling[i2];
                Body_Mass_Set[i_row][2] = Input_Mass_Variation[i3];
                i_row++;
            }
        }
    }

    Generate_New();

    Update_Generator_Methods();
}

void Body_Mass_Trophic_Scaling_Generator::Generate_New()
{ 
    if(!(Container_Pointer->Food_Web==NULL))
    {
        std::vector<double> Trophic_Position;

        switch(Body_Mass_Scaling_Method)
        {
        case(0):
        {
            Trophic_Position = Container_Pointer->Food_Web->Get_Value().Max_Trophic_Position;
        }
            break;
        case(1):
        {
            Trophic_Position = Container_Pointer->Food_Web->Get_Value().Min_Trophic_Position;
        }
            break;
        case(2):
        {
            Trophic_Position = Container_Pointer->Food_Web->Get_Value().Avg_Trophic_Position;
        }
            break;
        }

        Body_Mass.resize(Trophic_Position.size());

        if(Body_Mass_Set[iSet][2]>0)
        {
            boost::random::random_device true_rand;
            boost::random::mt19937 mt_generator(true_rand());
            boost::random::uniform_real_distribution< double > rand_mass(-Body_Mass_Set[iSet][2],Body_Mass_Set[iSet][2]);

            for(int i=0;i<Trophic_Position.size();i++)
            {
                Body_Mass[i]=Body_Mass_Set[iSet][0]*pow(pow(10,Body_Mass_Set[iSet][1]),Trophic_Position[i])*(1+rand_mass(mt_generator));
            }
        }
        else
        {
            for(int i=0;i<Trophic_Position.size();i++)
            {
                Body_Mass[i]=Body_Mass_Set[iSet][0]*pow(pow(10,Body_Mass_Set[iSet][1]),Trophic_Position[i]);
            }
        }
    }
    else
    {
        Body_Mass.resize(1);
        Body_Mass[0]=0;
    }
}

void Body_Mass_Trophic_Scaling_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Trophic_Scaling_Body_Mass");
    Generator_Data_Labels.push_back("Base_Body_Mass");
    Generator_Data_Labels.push_back("Trophic_Scaling");
    Generator_Data_Labels.push_back("Body_Mass_Variation");
    Generator_Data_Labels.push_back("Body_Mass_Scaling_Method");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Body_Mass_Set[iSet][0]);
    Generator_Data.push_back(Body_Mass_Set[iSet][1]);
    Generator_Data.push_back(Body_Mass_Set[iSet][2]);
    Generator_Data.push_back(Body_Mass_Scaling_Method);
}

Body_Mass_Trophic_Scaling_Generator* Body_Mass_Trophic_Scaling_Generator::Clone()
{
    return new Body_Mass_Trophic_Scaling_Generator(*this);
}

//Metabolism

Metabolism_Generator::Metabolism_Generator(Parameter_Generators* Container,std::vector<double> Input_ar, std::vector<double> Input_ax, std::vector<double> Input_y, std::vector<double> Input_e):
    Parameter_Generator(Container,1),ar(Input_ar),ax(Input_ax),y(Input_y),e(Input_e)
{
    Update_Generator_Methods();
}

std::vector<double> Metabolism_Generator::Get_ar()
{
    return ar;
}

std::vector<double> Metabolism_Generator::Get_ax()
{
    return ax;
}

std::vector<double> Metabolism_Generator::Get_y()
{
    return y;
}

std::vector<double> Metabolism_Generator::Get_e()
{
    return e;
}

void Metabolism_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Metabolism");

    Generator_Data.push_back(iSet);
}

Metabolism_Generator* Metabolism_Generator::Clone()
{
    return new Metabolism_Generator(*this);
}

Metabolism_Fixed_Species_Generator::Metabolism_Fixed_Species_Generator(Parameter_Generators* Container,int Input_nSpecies, std::vector<double> Input_ar, std::vector<double> Input_ax, std::vector<double> Input_y, std::vector<double> Input_e):
    Metabolism_Generator(Container,Input_ar.size()*Input_ax.size()*Input_y.size()*Input_e.size()),nSpecies(Input_nSpecies)
{
    Metabolism_Set.resize(boost::extents[nSets][4]);

    int i_row = 0;

    for(int i=0;i<Input_ar.size();i++)
    {
        for(int i2=0; i2<Input_ax.size();i2++)
        {
            for(int i3=0; i3<Input_y.size();i3++)
            {
                for(int i4=0; i4<Input_e.size();i4++)
                {
                    Metabolism_Set[i_row][0] = Input_ar[i];
                    Metabolism_Set[i_row][1] = Input_ax[i2];
                    Metabolism_Set[i_row][2] = Input_y[i3];
                    Metabolism_Set[i_row][3] = Input_e[i4];
                    i_row++;
                }
            }
        }
    }

    Generate_New();

    Update_Generator_Methods();
}

void Metabolism_Fixed_Species_Generator::Generate_New()
{
    ar.resize(nSpecies);
    ax.resize(nSpecies);
    y.resize(nSpecies);
    e.resize(nSpecies);

    for(int i=0;i<nSpecies;i++)
    {
        ar[i] = Metabolism_Set[iSet][0];
        ax[i] = Metabolism_Set[iSet][1];
        y[i] = Metabolism_Set[iSet][2];
        e[i] = Metabolism_Set[iSet][3];
    }
}

void Metabolism_Fixed_Species_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Species_Fixed_Metabolism");
    Generator_Data_Labels.push_back("ar");
    Generator_Data_Labels.push_back("ax");
    Generator_Data_Labels.push_back("y");
    Generator_Data_Labels.push_back("e");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Metabolism_Set[iSet][0]);
    Generator_Data.push_back(Metabolism_Set[iSet][1]);
    Generator_Data.push_back(Metabolism_Set[iSet][2]);
    Generator_Data.push_back(Metabolism_Set[iSet][3]);
}

Metabolism_Fixed_Species_Generator* Metabolism_Fixed_Species_Generator::Clone()
{
    return new Metabolism_Fixed_Species_Generator(*this);
}

//Functional Response

Functional_Response_Generator::Functional_Response_Generator(Parameter_Generators* Container,std::vector<double> Input_B0, std::vector<double> Input_h):
    Parameter_Generator(Container,1),B0(Input_B0), h(Input_h)
{
    Update_Generator_Methods();
}

std::vector<double> Functional_Response_Generator::Get_B0()
{
    return B0;
}

std::vector<double> Functional_Response_Generator::Get_h()
{
    return h;
}

void Functional_Response_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Functional_Response");

    Generator_Data.push_back(iSet);
}

Functional_Response_Generator* Functional_Response_Generator::Clone()
{
    return new Functional_Response_Generator(*this);
}

//Functional Response Derived

Functional_Response_Fixed_Species_Generator::Functional_Response_Fixed_Species_Generator(Parameter_Generators* Container,int Input_nSpecies, std::vector<double> Input_B0, std::vector<double> Input_h):
    Functional_Response_Generator(Container,Input_B0.size()*Input_h.size()),nSpecies(Input_nSpecies)
{
    Functional_Response_Set.resize(boost::extents[nSets][2]);

    int i_row = 0;

    for(int i=0;i<Input_B0.size();i++)
    {
        for(int i2=0;i2<Input_h.size();i2++)
        {
            Functional_Response_Set[i_row][0] = Input_B0[i];
            Functional_Response_Set[i_row][1] = Input_h[i2];
            i_row++;
        }
    }

    Generate_New();

    Update_Generator_Methods();
}

void Functional_Response_Fixed_Species_Generator::Generate_New()
{
    B0.resize(nSpecies);
    h.resize(nSpecies);

    for(int i=0;i<nSpecies;i++)
    {
        B0[i] = Functional_Response_Set[iSet][0];
        h[i] = Functional_Response_Set[iSet][1];
    }
}

void Functional_Response_Fixed_Species_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Species_Fixed_Functional_Response");
    Generator_Data_Labels.push_back("B0");
    Generator_Data_Labels.push_back("h");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Functional_Response_Set[iSet][0]);
    Generator_Data.push_back(Functional_Response_Set[iSet][1]);
}

Functional_Response_Fixed_Species_Generator* Functional_Response_Fixed_Species_Generator::Clone()
{
    return new Functional_Response_Fixed_Species_Generator(*this);
}

//Dispersal

Dispersal_Generator::Dispersal_Generator(Parameter_Generators* Container,std::vector<double> Input_Dispersal):
    Parameter_Generator(Container,1),Dispersal(Input_Dispersal)
{
    Update_Generator_Methods();
}

std::vector<double> Dispersal_Generator::Get_Value()
{
    return Dispersal;
}

void Dispersal_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Metabolism");

    Generator_Data.push_back(iSet);
}

Dispersal_Generator* Dispersal_Generator::Clone()
{
    return new Dispersal_Generator(*this);
}

//Dispersal Derived

Dispersal_Allometric_Linear_Scaling_Generator::Dispersal_Allometric_Linear_Scaling_Generator(Parameter_Generators* Container, std::vector<double> Input_Dispersal_Base, std::vector<double> Input_Trophic_Scaling, std::vector<double> Input_Competitive_Scaling, std::vector<double> Input_Random_Scaling):
    Dispersal_Generator(Container,Input_Dispersal_Base.size()*Input_Trophic_Scaling.size()*Input_Competitive_Scaling.size()*Input_Random_Scaling.size())
{
    Dispersal_Set.resize(boost::extents[nSets][4]);

    int i_row=0;

    for(int i=0;i<Input_Dispersal_Base.size();i++)
    {
        for(int i2=0;i2<Input_Trophic_Scaling.size();i2++)
        {
            for(int i3=0;i3<Input_Competitive_Scaling.size();i3++)
            {
                for(int i4=0;i4<Input_Random_Scaling.size();i4++)
                {
                    Dispersal_Set[i_row][0] = Input_Dispersal_Base[i];
                    Dispersal_Set[i_row][1] = Input_Trophic_Scaling[i2];
                    Dispersal_Set[i_row][2] = Input_Competitive_Scaling[i3];
                    Dispersal_Set[i_row][3] = Input_Random_Scaling[i4];
                    i_row++;
                }
            }
        }
    }

    Generate_New();

    Update_Generator_Methods();
}

void Dispersal_Allometric_Linear_Scaling_Generator::Generate_New()
{
    if(!(Container_Pointer->Food_Web==NULL))
    {
        std::vector<double> Trophic_Position = Container_Pointer->Food_Web->Get_Value().Max_Trophic_Position;
        std::vector<double> Generality = Container_Pointer->Food_Web->Get_Value().Generality;

        boost::random::random_device true_rand;
        boost::random::mt19937 mt_generator(true_rand());
        boost::random::uniform_real_distribution< double > rand_dispersal(0,1);


        std::vector<double> z_Random;

        Dispersal.resize(Trophic_Position.size());

        for(int i=0;i<Trophic_Position.size();i++)
        {

            z_Random.push_back((rand_dispersal(mt_generator)-.5)/sqrt((double)1/12));

            Dispersal[i] = Dispersal_Set[iSet][0]*(1+Dispersal_Set[iSet][1]*Trophic_Position[i]+Dispersal_Set[iSet][2]*Generality[i]+Dispersal_Set[iSet][3]*z_Random[i]+1);

            if(Dispersal[i]<0)
                Dispersal[i] = 0;
        }
    }
    else
    {
        Dispersal.resize(1);
        Dispersal[0]=0;
    }
}

void Dispersal_Allometric_Linear_Scaling_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Linear_Scaling_Dispersal");
    Generator_Data_Labels.push_back("Dispersal_Base");
    Generator_Data_Labels.push_back("Trophic_Scaling");
    Generator_Data_Labels.push_back("Generality_Scaling");
    Generator_Data_Labels.push_back("Random_Scaling");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Dispersal_Set[iSet][0]);
    Generator_Data.push_back(Dispersal_Set[iSet][1]);
    Generator_Data.push_back(Dispersal_Set[iSet][2]);
    Generator_Data.push_back(Dispersal_Set[iSet][3]);
}

Dispersal_Allometric_Linear_Scaling_Generator* Dispersal_Allometric_Linear_Scaling_Generator::Clone()
{
    return new Dispersal_Allometric_Linear_Scaling_Generator(*this);
}

// Dispersal Derived II

Dispersal_Allometric_Exponential_Scaling_Generator::Dispersal_Allometric_Exponential_Scaling_Generator(Parameter_Generators* Container, std::vector<double> Input_Dispersal_Base, std::vector<double> Input_Trophic_Scaling, std::vector<double> Input_Competitive_Scaling, std::vector<double> Input_Random_Scaling):
    Dispersal_Generator(Container,Input_Dispersal_Base.size()*Input_Trophic_Scaling.size()*Input_Competitive_Scaling.size()*Input_Random_Scaling.size())
{
    Dispersal_Set.resize(boost::extents[nSets][4]);

    int i_row=0;

    for(int i=0;i<Input_Dispersal_Base.size();i++)
    {
        for(int i2=0;i2<Input_Trophic_Scaling.size();i2++)
        {
            for(int i3=0;i3<Input_Competitive_Scaling.size();i3++)
            {
                for(int i4=0;i4<Input_Random_Scaling.size();i4++)
                {
                    Dispersal_Set[i_row][0] = Input_Dispersal_Base[i];
                    Dispersal_Set[i_row][1] = Input_Trophic_Scaling[i2];
                    Dispersal_Set[i_row][2] = Input_Competitive_Scaling[i3];
                    Dispersal_Set[i_row][3] = Input_Random_Scaling[i4];
                    i_row++;
                }
            }
        }
    }

    Generate_New();

    Update_Generator_Methods();
}

void Dispersal_Allometric_Exponential_Scaling_Generator::Generate_New()
{
    if(!(Container_Pointer->Food_Web==NULL))
    {
        std::vector<double> Trophic_Position = Container_Pointer->Food_Web->Get_Value().Max_Trophic_Position;
        std::vector<double> Generality = Container_Pointer->Food_Web->Get_Value().Generality;

        boost::random::random_device true_rand;
        boost::random::mt19937 mt_generator(true_rand());
        boost::random::uniform_real_distribution< double > rand_dispersal(0,1);

        std::vector<double> z_Trophic;
        std::vector<double> z_Generality;
        std::vector<double> z_Random;

        Dispersal.resize(Trophic_Position.size());

        for(int i=0;i<Trophic_Position.size();i++)
        {
            z_Trophic.push_back((Trophic_Position[i]-.5)/sqrt((double)1/12));
            z_Generality.push_back((Generality[i]-.5)/sqrt((double)1/12));
            z_Random.push_back((rand_dispersal(mt_generator)-.5)/sqrt((double)1/12));

            Dispersal[i] = Dispersal_Set[iSet][0]*(1+pow(1+Dispersal_Set[iSet][1],Trophic_Position[i])-1+pow(1+Dispersal_Set[iSet][2],Generality[i])-1+pow(1+Dispersal_Set[iSet][3],z_Random[i]+1)-1);

            if(Dispersal[i]<0)
                Dispersal[i] = 0;
        }
    }
    else
    {
        Dispersal.resize(1);
        Dispersal[0]=0;
    }
}

void Dispersal_Allometric_Exponential_Scaling_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Exponential_Scaling_Dispersal");
    Generator_Data_Labels.push_back("Dispersal_Base");
    Generator_Data_Labels.push_back("Trophic_Scaling");
    Generator_Data_Labels.push_back("Generality_Scaling");
    Generator_Data_Labels.push_back("Random_Scaling");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Dispersal_Set[iSet][0]);
    Generator_Data.push_back(Dispersal_Set[iSet][1]);
    Generator_Data.push_back(Dispersal_Set[iSet][2]);
    Generator_Data.push_back(Dispersal_Set[iSet][3]);
}


Dispersal_Allometric_Exponential_Scaling_Generator* Dispersal_Allometric_Exponential_Scaling_Generator::Clone()
{
    return new Dispersal_Allometric_Exponential_Scaling_Generator(*this);
}

//Carrying Capactiy

Carrying_Capacity_Generator::Carrying_Capacity_Generator(Parameter_Generators* Container,boost_matrix Input_Carrying_Capacity,std::vector<double> Input_Patch_Size):
    Parameter_Generator(Container,1),Carrying_Capacity(Input_Carrying_Capacity),Patch_Size(Input_Patch_Size)
{
    Update_Generator_Methods();
}

boost_matrix Carrying_Capacity_Generator::Get_K()
{
    if(Container_Pointer->Food_Web->Get_Value().Food_Web.shape()[0]!=Carrying_Capacity.shape()[0]||Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]!=Carrying_Capacity.shape()[1])
    {
        boost_matrix k_temp(boost::extents[Container_Pointer->Food_Web->Get_Value().Food_Web.shape()[0]][Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]]);
        for(int i=0;i<k_temp.shape()[0];i++)
            for(int i2=0;i2<k_temp.shape()[1];i2++)
            {
                if(i>Carrying_Capacity.shape()[0]||i2>Carrying_Capacity.shape()[1])
                    k_temp[i][i2] = 1;
                else
                    k_temp[i][i2] = Carrying_Capacity[i][i2];
            }

        return k_temp;
    }

    return Carrying_Capacity;
}

std::vector<double> Carrying_Capacity_Generator::Get_Patch_Size()
{
    if(Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]!=Patch_Size.size())
    {
        if(Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]<Patch_Size.size())
        {
            std::vector<double> Patch_Size_temp;
            for(int i=0;i<Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0];i++)
                Patch_Size_temp.push_back(Patch_Size[i]);

            return Patch_Size_temp;
        }
        else
        {
            std::vector<double> Patch_Size_temp = Patch_Size;
            for(int i=0;i<Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]-Patch_Size.size();i++)
                Patch_Size_temp.push_back(1);

            return Patch_Size_temp;
        }
    }

    return Patch_Size;
}

void Carrying_Capacity_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Carrying_Capacity");

    Generator_Data.push_back(iSet);
}

Carrying_Capacity_Generator* Carrying_Capacity_Generator::Clone()
{
    return new Carrying_Capacity_Generator(*this);
}

//Carrying Capacity Derived

Carrying_Capacity_Fixed_Community_Generator::Carrying_Capacity_Fixed_Community_Generator(Parameter_Generators* Container,int Input_nSpecies, int Input_nPatches, std::vector<double> Input_Carrying_Capacity):
   Carrying_Capacity_Generator(Container,Input_Carrying_Capacity.size()),nSpecies(Input_nSpecies), nPatches(Input_nPatches), Fixed_Carrying_Capacity(Input_Carrying_Capacity)
{    
    Generate_New();

    Update_Generator_Methods();
}

void Carrying_Capacity_Fixed_Community_Generator::Generate_New()
{
    Carrying_Capacity.resize(boost::extents[nSpecies][nPatches]);
    Patch_Size.resize(nPatches);

    for(int i=0;i<nSpecies;i++)
    {
        for(int i2=0;i2<nPatches;i2++)
        {
            Carrying_Capacity[i][i2] = Fixed_Carrying_Capacity[iSet];
        }
    }

    for(int i=0;i<nPatches;i++)
        Patch_Size[i] = Fixed_Carrying_Capacity[iSet];
}

void Carrying_Capacity_Fixed_Community_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Community_Fixed_Carrying_Capacity");
    Generator_Data_Labels.push_back("Carrying_Capacity");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Fixed_Carrying_Capacity[iSet]);
}

Carrying_Capacity_Fixed_Community_Generator* Carrying_Capacity_Fixed_Community_Generator::Clone()
{
    return new Carrying_Capacity_Fixed_Community_Generator(*this);
}

//Carrying Capacity Derived II

Carrying_Capacity_Scale_By_Size_Generator::Carrying_Capacity_Scale_By_Size_Generator(Parameter_Generators* Container,int Input_nSpecies, boost_matrix Input_Patch_Size_Set):
   Carrying_Capacity_Generator(Container,Input_Patch_Size_Set.shape()[0]),nSpecies(Input_nSpecies), Patch_Size_Set(Input_Patch_Size_Set)
{
    Generate_New();

    Update_Generator_Methods();
}

void Carrying_Capacity_Scale_By_Size_Generator::Generate_New()
{
    Carrying_Capacity.resize(boost::extents[nSpecies][Patch_Size_Set.shape()[1]]);
    Patch_Size.resize(Patch_Size_Set.shape()[1]);

    for(int i=0;i<nSpecies;i++)
    {
        for(int i2=0;i2<Patch_Size_Set.shape()[1];i2++)
        {
            Carrying_Capacity[i][i2] = 1;
        }
    }

    for(int i=0;i<Patch_Size_Set.shape()[1];i++)
        Patch_Size[i] = Patch_Size_Set[iSet][i];
}

void Carrying_Capacity_Scale_By_Size_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Scale_By_Size_Carrying_Capacity");

    Generator_Data.push_back(iSet);
}

Carrying_Capacity_Scale_By_Size_Generator* Carrying_Capacity_Scale_By_Size_Generator::Clone()
{
    return new Carrying_Capacity_Scale_By_Size_Generator(*this);
}

//Initial Abundances

Initial_Abundance_Generator::Initial_Abundance_Generator(Parameter_Generators* Container,boost_matrix Input_Initial_Abundance):
    Parameter_Generator(Container,1),Initial_Abundance(Input_Initial_Abundance)
{
    Update_Generator_Methods();
}

boost_matrix Initial_Abundance_Generator::Get_Value()
{
    if(Container_Pointer->Food_Web->Get_Value().Food_Web.shape()[0]!=Initial_Abundance.shape()[0]||Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]!=Initial_Abundance.shape()[1])
    {
        boost_matrix ini_temp(boost::extents[Container_Pointer->Food_Web->Get_Value().Food_Web.shape()[0]][Container_Pointer->Spatial_Structure->Get_Value().Spatial_Structure.shape()[0]]);
        for(int i=0;i<ini_temp.shape()[0];i++)
            for(int i2=0;i2<ini_temp.shape()[1];i2++)
            {
                if(i>Initial_Abundance.shape()[0]||i2>Initial_Abundance.shape()[1])
                    ini_temp[i][i2] = 1;
                else
                    ini_temp[i][i2] = Initial_Abundance[i][i2];
            }

        return ini_temp;
    }

    return Initial_Abundance;
}

void Initial_Abundance_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Abundance");

    Generator_Data.push_back(iSet);
}

Initial_Abundance_Generator* Initial_Abundance_Generator::Clone()
{
    return new Initial_Abundance_Generator(*this);
}

//Initial Abundance Derived

Initial_Abundance_Random_Generator::Initial_Abundance_Random_Generator(Parameter_Generators* Container, int Input_nPatches, std::vector<double> Input_Initial_Abundance_Mean, std::vector<double> Input_Initial_Abundance_Dev):
    Initial_Abundance_Generator(Container,1),nSpecies(Input_Initial_Abundance_Mean.size()),nPatches(Input_nPatches),Initial_Abundance_Mean(Input_Initial_Abundance_Mean),Initial_Abundance_Dev(Input_Initial_Abundance_Dev)
{    
    Generate_New();
    Update_Generator_Methods();
}

void Initial_Abundance_Random_Generator::Generate_New()
{
    Initial_Abundance.resize(boost::extents[nSpecies][nPatches]);

    for(int i=0;i<nSpecies;i++)
    {
        double Minimum = Initial_Abundance_Mean[i]-Initial_Abundance_Dev[i];

        if(Minimum<0)
            Minimum=0;

        double Maximum = Initial_Abundance_Mean[i]+Initial_Abundance_Dev[i];

        boost::random::random_device true_rand;
        boost::random::mt19937 mt_generator(true_rand());
        boost::random::uniform_real_distribution< double > rand_abundance(Minimum,Maximum);

        for(int i2=0;i2<nPatches;i2++)
        {
            Initial_Abundance[i][i2] = rand_abundance(mt_generator);
        }
    }
}

void Initial_Abundance_Random_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Random_Abundance");

    Generator_Data.push_back(iSet);
}

Initial_Abundance_Random_Generator* Initial_Abundance_Random_Generator::Clone()
{
    return new Initial_Abundance_Random_Generator(*this);
}

//Extinction Threshold

Extinction_Threshold_Generator::Extinction_Threshold_Generator(Parameter_Generators* Container,double Input_Extinction_Threshold):
    Parameter_Generator(Container,1), Extinction_Threshold(Input_Extinction_Threshold)
{
    Update_Generator_Methods();
}

double Extinction_Threshold_Generator::Get_Value()
{
    return Extinction_Threshold;
}

void Extinction_Threshold_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Fixed_Extinction_Threshold");

    Generator_Data.push_back(iSet);
}

Extinction_Threshold_Generator* Extinction_Threshold_Generator::Clone()
{
    return new Extinction_Threshold_Generator(*this);
}

//Extinction Threshold Derived

Extinction_Threshold_Set_Generator::Extinction_Threshold_Set_Generator(Parameter_Generators* Container, std::vector<double> Extinction_Threshold_Set_Input):
    Extinction_Threshold_Generator(Container,(int)Extinction_Threshold_Set_Input.size()),Extinction_Threshold_Set(Extinction_Threshold_Set_Input)
{
    Generate_New();

    Update_Generator_Methods();
}

void Extinction_Threshold_Set_Generator::Generate_New()
{
    Extinction_Threshold = Extinction_Threshold_Set[iSet];
}

void Extinction_Threshold_Set_Generator::Update_Generator_Methods()
{
    Generator_Data_Labels.clear();
    Generator_Data.clear();

    Generator_Data_Labels.push_back("Set_Extinction");
    Generator_Data_Labels.push_back("Extinction_Threshold");

    Generator_Data.push_back(iSet);
    Generator_Data.push_back(Extinction_Threshold_Set[iSet]);
}

Extinction_Threshold_Set_Generator* Extinction_Threshold_Set_Generator::Clone()
{
    return new Extinction_Threshold_Set_Generator(*this);
}
