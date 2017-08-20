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
#ifndef PARAMETER_GENERATOR_FUNCTIONS
#define PARAMETER_GENERATOR_FUNCTIONS

#include "food_web_methods.h"
#include "spatial_structure_methods.h"
// <boost/multi_array.hpp> & <vector> from "food_web_methods.h"
// boost_matrix definition from "food_web_methods.h";
#include <string>

//Container

//class Parameter_Generator_Container
//{
//protected:
//    int nGenerators;
//    std::vector<Parameter_Generator*> Parameter_Generators;

//public:

//    Parameter_Generator_Container()
//    {
//        nGenerators = 0;
//    }

//    Parameter_Generator_Container(const Parameter_Generator_Container &obj)
//    {
//        for(int i=0;i<nGenerators;i++)
//        {
//            Parameter_Generators.push_back(obj.Parameter_Generators[i]->Clone());
//        }
//    }

//    void Add_Generator(Parameter_Generator* to_Add)
//    {
//        nGenerators++;
//        Parameter_Generators.push_back(to_Add);
//    }

//    Parameter_Generator* operator[](int i)
//    {
//        if(i<nGenerators)
//            return Parameter_Generators[i];
//        else
//        {
//            Parameter_Generator* Empty_Return = NULL;
//            return Empty_Return;
//        }
//    }
//};

//Make this into a class, make it generalizable, and make something more generalizable for Eqn_Parameters to bridge it

struct Parameter_Generators;

//Base Class

class Parameter_Generator
{
public:
    //Parameter_Generator();
    Parameter_Generator(Parameter_Generators* Container,int In_nSets);
    int Current_Set() const;
    int Get_nSets() const;
    std::vector<std::string> Get_Generator_Method_Labels() const;
    std::vector<double> Get_Generator_Method_Data() const;
    virtual bool Increment_Set();
    virtual void Reset();
    virtual void Update_Generator_Methods()=0;
    void Set_Container(Parameter_Generators* New_Container);
    virtual void Generate_New(){}
    virtual Parameter_Generator* Clone()=0;
    virtual ~Parameter_Generator(){}


protected:

    int nSets;
    int iSet;

    Parameter_Generators* Container_Pointer;

    std::vector<std::string> Generator_Data_Labels;
    std::vector<double> Generator_Data;

};

//Food Web Default

class Food_Web_Generator: public Parameter_Generator
{

public:
    Food_Web_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container,In_nSets){}
    Food_Web_Generator(Parameter_Generators* Container,boost_matrix Input);
    Food_Web_Metadata Get_Value();
    void Set_Value(const Food_Web_Metadata &Input);
    void Update_Generator_Methods();
    Food_Web_Generator* Clone();
    ~Food_Web_Generator(){}

protected:

    Food_Web_Metadata Food_Web_Data;

};

//Food Web Derived

class Food_Web_Niche_Model_Generator: public Food_Web_Generator
{

public:

    Food_Web_Niche_Model_Generator(Parameter_Generators* Container,int Input_nSpecies, std::vector<double> Input_Connectance);
    void Generate_New();
    void Update_Generator_Methods();
    Food_Web_Niche_Model_Generator* Clone();
    ~Food_Web_Niche_Model_Generator(){}

protected:

    int nSpecies;

    std::vector<double> Connectance;
};

//Food Web Derived II

class Food_Web_Read_Generator: public Food_Web_Generator
{

public:

    Food_Web_Read_Generator(Parameter_Generators* Container,std::vector<boost_matrix> Input_Food_Webs);
    void Generate_New();
    void Update_Generator_Methods();
    Food_Web_Read_Generator* Clone();
    ~Food_Web_Read_Generator(){}

protected:

    int Interset_Iter;
    std::vector<boost_matrix> Read_Food_Webs;
};

// Spatial Structure Default

class Spatial_Structure_Generator: public Parameter_Generator
{
public:
    Spatial_Structure_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container, In_nSets){}
    Spatial_Structure_Generator(Parameter_Generators* Container,boost_matrix Input);
    Spatial_Structure_Metadata Get_Value();
    void Set_Value(const Spatial_Structure_Metadata &Input);
    void Update_Generator_Methods();
    Spatial_Structure_Generator* Clone();
    ~Spatial_Structure_Generator(){}

protected:

    Spatial_Structure_Metadata Spatial_Structure_Data;

};

// Spatial Structure Derived

class Spatial_Structure_Undirected_Erdos_Renyi_Generator: public Spatial_Structure_Generator
{

public:

    Spatial_Structure_Undirected_Erdos_Renyi_Generator(Parameter_Generators* Container,int Input_nPatches,std::vector<double> Input_Connectance);
    void Generate_New();
    void Update_Generator_Methods();
    Spatial_Structure_Undirected_Erdos_Renyi_Generator* Clone();
    ~Spatial_Structure_Undirected_Erdos_Renyi_Generator(){}

protected:

    int nPatches;

    std::vector<double> Connectance;
};

//Spatial Structure Derived II

class Spatial_Structure_Undirected_Dendritic_Generator: public Spatial_Structure_Generator
{

public:

    Spatial_Structure_Undirected_Dendritic_Generator(Parameter_Generators* Container,int Input_nPatches,std::vector<double> Input_Branching_Probability);
    void Generate_New();
    void Update_Generator_Methods();
    Spatial_Structure_Undirected_Dendritic_Generator* Clone();
    ~Spatial_Structure_Undirected_Dendritic_Generator(){}

protected:

    int nPatches;

    std::vector<double> Branching_Probability;
};

//Spatial Structure Derived III

class Spatial_Structure_Kurt_Comparison_Generator: public Spatial_Structure_Generator
{

public:

    Spatial_Structure_Kurt_Comparison_Generator(Parameter_Generators* Container,int Input_nPatches,std::vector<int> Input_nBranches,std::vector<int> Input_Neighbor_Distance);
    void Generate_New();
    void Update_Generator_Methods();
    Spatial_Structure_Kurt_Comparison_Generator* Clone();
    ~Spatial_Structure_Kurt_Comparison_Generator(){}

protected:

    int nPatches;

    std::vector<int> nBranches;
    std::vector<int> Neighbor_Distance;
};

//Spatial Structure Derived IV

//Spatial_Structure_Parameter_Method = new Spatial_Structure_Shipley_Skinner_Generator(&Input_Cache,nRewire,nRescale,Rescale_Factor);

class Spatial_Structure_Shipley_Skinner_Generator: public Spatial_Structure_Generator
{

public:

    Spatial_Structure_Shipley_Skinner_Generator(Parameter_Generators* Container,std::vector<int> Input_nRewire);
    void Generate_New();
    void Update_Generator_Methods();
    Spatial_Structure_Shipley_Skinner_Generator* Clone();
    ~Spatial_Structure_Shipley_Skinner_Generator(){}

protected:

    std::vector<int> nRewire;

    int Interset_Iter;
    int Previous_iSet;
    std::vector<std::vector<std::vector<int>>> Remove_Sets;
};

//Spatial Structure Derived V

//Spatial_Structure_Parameter_Method = new Spatial_Structure_Shipley_Skinner_Generator(&Input_Cache,nRewire,nRescale,Rescale_Factor);

class Spatial_Structure_Read_Generator: public Spatial_Structure_Generator
{

public:

    Spatial_Structure_Read_Generator(Parameter_Generators* Container,std::vector<boost_matrix> Input_Spatial_Structures);
    void Generate_New();
    void Update_Generator_Methods();
    Spatial_Structure_Read_Generator* Clone();
    ~Spatial_Structure_Read_Generator(){}

protected:

    int Interset_Iter;
    std::vector<boost_matrix> Read_Spatial_Structures;
};

//Body Mass Default

class Body_Mass_Generator: public Parameter_Generator
{
public:
    Body_Mass_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container,In_nSets){}
    Body_Mass_Generator(Parameter_Generators* Container,std::vector<double> Input);
    std::vector<double> Get_Value();
    void Update_Generator_Methods();
    Body_Mass_Generator* Clone();
    ~Body_Mass_Generator(){}

protected:

    std::vector<double> Body_Mass;

};

//Body Mass Derived

class Body_Mass_Trophic_Scaling_Generator: public Body_Mass_Generator
{
public:
    Body_Mass_Trophic_Scaling_Generator(Parameter_Generators* Container, int Input_Body_Mass_Scaling_Method, std::vector<double> Input_Base_Mass, std::vector<double> Input_Trophic_Scaling, std::vector<double> Input_Mass_Variation);
    void Generate_New();
    void Update_Generator_Methods();
    Body_Mass_Trophic_Scaling_Generator* Clone();
    ~Body_Mass_Trophic_Scaling_Generator(){}

protected:

    boost_matrix Body_Mass_Set;
    int Body_Mass_Scaling_Method;
};

//Metabolism Default

class Metabolism_Generator: public Parameter_Generator
{
public:
    Metabolism_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container, In_nSets){}
    Metabolism_Generator(Parameter_Generators* Container,std::vector<double> Input_ar, std::vector<double> Input_ax, std::vector<double> Input_y, std::vector<double> Input_e);
    std::vector<double> Get_ar();
    std::vector<double> Get_ax();
    std::vector<double> Get_y();
    std::vector<double> Get_e();
    void Update_Generator_Methods();
    Metabolism_Generator* Clone();
    ~Metabolism_Generator(){}

protected:
    std::vector<double> ar;
    std::vector<double> ax;
    std::vector<double> y;
    std::vector<double> e;

};

// Metabolism Derived

class Metabolism_Fixed_Species_Generator: public Metabolism_Generator
{
public:
    Metabolism_Fixed_Species_Generator(Parameter_Generators* Container,int Input_nSpecies,std::vector<double> Input_ar, std::vector<double> Input_ax, std::vector<double> Input_y, std::vector<double> Input_e);
    void Generate_New();
    void Update_Generator_Methods();
    Metabolism_Fixed_Species_Generator* Clone();
    ~Metabolism_Fixed_Species_Generator(){}

protected:

    int nSpecies;

    boost_matrix Metabolism_Set;

};

//Functional Response Default

class Functional_Response_Generator: public Parameter_Generator
{
public:
    Functional_Response_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container, In_nSets){}
    Functional_Response_Generator(Parameter_Generators* Container,std::vector<double> Input_B0, std::vector<double> Input_h);
    std::vector<double> Get_B0();
    std::vector<double> Get_h();
    void Update_Generator_Methods();
    Functional_Response_Generator* Clone();
    ~Functional_Response_Generator(){}

protected:
    std::vector<double> B0;
    std::vector<double> h;
};

//Functional Response Derived

class Functional_Response_Fixed_Species_Generator: public Functional_Response_Generator
{
public:
    Functional_Response_Fixed_Species_Generator(Parameter_Generators* Container,int Input_nSpecies, std::vector<double> Input_B0, std::vector<double> Input_h);
    void Generate_New();
    void Update_Generator_Methods();
    Functional_Response_Fixed_Species_Generator* Clone();
    ~Functional_Response_Fixed_Species_Generator(){}
protected:
    int nSpecies;

    boost_matrix Functional_Response_Set;
};

//Dispersal Default

class Dispersal_Generator: public Parameter_Generator
{
public:
    Dispersal_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container, In_nSets){}
    Dispersal_Generator(Parameter_Generators* Container,std::vector<double> Input_Dispersal);
    std::vector<double> Get_Value();
    void Update_Generator_Methods();
    Dispersal_Generator* Clone();
    ~Dispersal_Generator(){}

protected:
    std::vector<double> Dispersal;
};

//Dispersal Derived

class Dispersal_Allometric_Linear_Scaling_Generator: public Dispersal_Generator
{
public:
    Dispersal_Allometric_Linear_Scaling_Generator(Parameter_Generators* Container, std::vector<double> Input_Dispersal_Base, std::vector<double> Input_Trophic_Scaling, std::vector<double> Input_Competitive_Scaling, std::vector<double> Input_Random_Scaling);
    void Generate_New();
    void Update_Generator_Methods();
    Dispersal_Allometric_Linear_Scaling_Generator* Clone();
    ~Dispersal_Allometric_Linear_Scaling_Generator(){}

protected:

    boost_matrix Dispersal_Set;

    //    std::vector<double> Dispersal_Base;
    //    std::vector<double> Trophic_Scaling;
    //    std::vector<double> Competitive_Scaling;
    //    std::vector<double> Random_Scaling;

};


class Dispersal_Allometric_Exponential_Scaling_Generator: public Dispersal_Generator
{
public:
    Dispersal_Allometric_Exponential_Scaling_Generator(Parameter_Generators* Container, std::vector<double> Input_Dispersal_Base, std::vector<double> Input_Trophic_Scaling, std::vector<double> Input_Competitive_Scaling, std::vector<double> Input_Random_Scaling);
    void Generate_New();
    void Update_Generator_Methods();
    Dispersal_Allometric_Exponential_Scaling_Generator* Clone();
    ~Dispersal_Allometric_Exponential_Scaling_Generator(){}

protected:

    boost_matrix Dispersal_Set;

    //    std::vector<double> Dispersal_Base;
    //    std::vector<double> Trophic_Scaling;
    //    std::vector<double> Competitive_Scaling;
    //    std::vector<double> Random_Scaling;

};

//Carrying Capacity Default

class Carrying_Capacity_Generator: public Parameter_Generator
{
public:
    Carrying_Capacity_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container,In_nSets){}
    Carrying_Capacity_Generator(Parameter_Generators* Container,boost_matrix Input_Carrying_Capacity,std::vector<double> Input_Patch_Size);
    boost_matrix Get_K();
    std::vector<double> Get_Patch_Size();
    void Update_Generator_Methods();
    Carrying_Capacity_Generator* Clone();
    ~Carrying_Capacity_Generator(){}

protected:
    boost_matrix Carrying_Capacity;
    std::vector<double> Patch_Size;
};

//Carrying Capacity Derived

class Carrying_Capacity_Fixed_Community_Generator: public Carrying_Capacity_Generator
{
public:
    Carrying_Capacity_Fixed_Community_Generator(Parameter_Generators* Container,int Input_nSpecies, int Input_nPatches, std::vector<double> Input_Carrying_Capacity);
    void Generate_New();
    void Update_Generator_Methods();
    Carrying_Capacity_Fixed_Community_Generator* Clone();
    ~Carrying_Capacity_Fixed_Community_Generator(){}

protected:

    int nSpecies;
    int nPatches;
    std::vector<double> Fixed_Carrying_Capacity;

};

class Carrying_Capacity_Scale_By_Size_Generator: public Carrying_Capacity_Generator
{
public:
    Carrying_Capacity_Scale_By_Size_Generator(Parameter_Generators* Container,int Input_nSpecies, boost_matrix Input_Patch_Size_Set);
    void Generate_New();
    void Update_Generator_Methods();
    Carrying_Capacity_Scale_By_Size_Generator* Clone();
    ~Carrying_Capacity_Scale_By_Size_Generator(){}

protected:

    int nSpecies;
    boost_matrix Patch_Size_Set;

};

//Initial Abundances

class Initial_Abundance_Generator: public Parameter_Generator
{
public:
    Initial_Abundance_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container,In_nSets){}
    Initial_Abundance_Generator(Parameter_Generators* Container,boost_matrix Input_Initial_Abundance);
    boost_matrix Get_Value();
    void Update_Generator_Methods();
    Initial_Abundance_Generator* Clone();
    ~Initial_Abundance_Generator(){}

protected:
    boost_matrix Initial_Abundance;

};

//Initial Abundance Derived

class Initial_Abundance_Random_Generator: public Initial_Abundance_Generator
{
public:
    Initial_Abundance_Random_Generator(Parameter_Generators* Container, int Input_nPatches, std::vector<double> Input_Initial_Abundance_Mean, std::vector<double> Input_Initial_Abundance_Dev);
    void Generate_New();
    void Update_Generator_Methods();
    Initial_Abundance_Random_Generator* Clone();
    ~Initial_Abundance_Random_Generator(){}

protected:

    int nSpecies;
    int nPatches;

    std::vector<double> Initial_Abundance_Mean;
    std::vector<double> Initial_Abundance_Dev;

};

//Extinction Threshold Default

class Extinction_Threshold_Generator: public Parameter_Generator
{
public:
    Extinction_Threshold_Generator(Parameter_Generators* Container, int In_nSets):
        Parameter_Generator(Container,In_nSets){}
    Extinction_Threshold_Generator(Parameter_Generators* Container,double Input_Extinction_Threshold);
    double Get_Value();
    void Update_Generator_Methods();
    Extinction_Threshold_Generator* Clone();
    ~Extinction_Threshold_Generator(){}

protected:
    double Extinction_Threshold;

};

//Extinction Threshold Derived

class Extinction_Threshold_Set_Generator: public Extinction_Threshold_Generator
{
public:
    Extinction_Threshold_Set_Generator(Parameter_Generators* Container, std::vector<double> Input_Extinction_Threshold_Set);
    void Generate_New();
    void Update_Generator_Methods();
    Extinction_Threshold_Set_Generator* Clone();
    ~Extinction_Threshold_Set_Generator(){}

protected:

    std::vector<double> Extinction_Threshold_Set;

};

struct Parameter_Generators{

    Food_Web_Generator* Food_Web;
    Spatial_Structure_Generator* Spatial_Structure;
    Body_Mass_Generator* Body_Mass;
    Metabolism_Generator* Metabolism;
    Functional_Response_Generator* Functional_Response;
    Dispersal_Generator* Dispersal;
    Carrying_Capacity_Generator* Carrying_Capacity;
    Initial_Abundance_Generator* Initial_Abundance;
    Extinction_Threshold_Generator* Extinction_Threshold;

    std::vector<std::string> Get_Generator_Set_Method_Labels() const;

    std::vector<double> Get_Generator_Set_Method_Data() const;

    std::string Check() const;

    //order matters in current iteration - values generated in order and some have dependencies (Body_Mass/Dispersal depend on Food_Web sometimes)
    Parameter_Generator* operator[](int i) const;

    int nGenerators() const;

    Parameter_Generators();


    Parameter_Generators(const Parameter_Generators &obj);


    ~Parameter_Generators();


};



#endif // PARAMETER_GENERATOR_FUNCTIONS

