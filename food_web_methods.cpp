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
#include "food_web_methods.h"
#include "spatial_structure_methods.h"

#include <algorithm>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <utility>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/connected_components.hpp>

#include <boost/graph/erdos_renyi_generator.hpp>

#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>

#include <boost/graph/isomorphism.hpp>

#include <Eigen/Dense>

#include <boost/multi_array.hpp>

#include <iostream>


double Log2( double n )
{
    // log(n)/log(2) is log2.
    return log( n ) / log( (double)2 );
}

Food_Web_Metadata::Food_Web_Metadata()
{

}

void Food_Web_Metadata::operator()(const Food_Web_Metadata &New_Values)
{
    H = New_Values.H;
    AMI = New_Values.AMI;
    Hc = New_Values.Hc;

    //H_in = New_Values.H_in;
    //H_out = New_Values.H_out;
    N_in = New_Values.N_in;
    N_out = New_Values.N_out;

    is_Producer = New_Values.is_Producer;
    is_Consumer = New_Values.is_Consumer;

    Min_Trophic_Position = New_Values.Min_Trophic_Position;
    Avg_Trophic_Position = New_Values.Avg_Trophic_Position;
    Max_Trophic_Position = New_Values.Max_Trophic_Position;
    Trophic_Positions = New_Values.Trophic_Positions;

    Generality = New_Values.Generality;
    Positional_Index = New_Values.Positional_Index;

    Food_Web.resize(boost::extents[New_Values.Food_Web.shape()[0]][New_Values.Food_Web.shape()[1]]);

    for(int i=0;i<New_Values.Food_Web.shape()[0];i++)
        for(int i2=0;i2<New_Values.Food_Web.shape()[1];i2++)
            Food_Web[i][i2]=New_Values.Food_Web[i][i2];
}

bool Food_Web_Metadata::is_Equal(const Food_Web_Metadata &Compare_To)
{   
    if(abs(H-Compare_To.H)<.00001)
        if(abs(AMI-Compare_To.AMI)<.00001)
            if(abs(Hc-Compare_To.Hc)<.00001)
                if(is_Isomorphic(Food_Web,Compare_To.Food_Web))
                    return true;
    return false;
}

boost_matrix Niche_Model(int nSpecies, double Connectance)
{
    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, double> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    int Test_Connected_Components = 0;

    Graph Gen_Graph;

    boost_matrix Food_Web;

    Food_Web.resize(boost::extents[nSpecies][nSpecies]);

    std::vector< double > Trophic_Value(nSpecies,0);
    std::vector< double > Generality(nSpecies,0);
    std::vector< double > ci(nSpecies,0);

    double beta  = (1/(2*Connectance))-1;

    boost::random::random_device true_rand;
    boost::random::mt19937 mt_generator(true_rand());
    boost::random::uniform_real_distribution< double > uni_0_1(0,1);

    while(Test_Connected_Components != 1)
    {
        Test_Connected_Components = 0;

        Gen_Graph.clear();

        std::fill(Food_Web.data(),Food_Web.data()+Food_Web.num_elements(),0);

        for(int i=0;i<nSpecies;i++)
        {
            Trophic_Value[i] = uni_0_1(mt_generator);

            Generality[i] = (1-pow((1-uni_0_1(mt_generator)),(1/beta)))*Trophic_Value[i];

            boost::random::uniform_real_distribution< double > uni_ci((Generality[i]/2),Trophic_Value[i]);

            ci[i] = uni_ci(mt_generator);

            add_vertex(Gen_Graph);
        }

        //Base Web Generation
        for(int i1 = 0;i1<nSpecies;i1++)
        {
            for(int i2 = 0;i2<nSpecies;i2++)
            {
                if(Trophic_Value[i2]>=(ci[i1]-(Generality[i1]/2))&&Trophic_Value[i2]<=(ci[i1]+(Generality[i1]/2))&&i1!=i2)
                {
                    add_edge(i1,i2,Gen_Graph);
                    Food_Web[i1][i2]=1;
                }
            }
        }

        //Check for identical Species
        bool All_Species_Unique = true;

        for(int i=0;i<nSpecies;i++)
            for(int j=i+1;j<nSpecies;j++)
            {
                bool Species_Unique = false;

                for(int k=0;k<nSpecies;k++)
                    if(Food_Web[k][i]!=Food_Web[k][j])
                    {
                        Species_Unique = true;
                        break;
                    }

                if(!Species_Unique)
                {
                    All_Species_Unique = false;
                    break;
                }
            }

        if(All_Species_Unique)
        {

            //Check for disconnected Species
            std::vector<std::vector<int>> Food_Chains = Calculate_Food_Chains(Food_Web);

            bool Resource_Check = true;

            for(int i = 0;i<nSpecies;i++)
                if(Food_Chains[i].empty())
                    Resource_Check = false;

            //Check all species connected
            if(Resource_Check)
            {
                bool Connected_Check = true;

                std::vector<double> a_col_sum(nSpecies,0);
                std::vector<double> a_row_sum(nSpecies,0);

                for(int i1 = 0;i1<nSpecies;i1++)
                {
                    a_col_sum[i1] = 0;
                    a_row_sum[i1] = 0;

                    for(int i2 = 0;i2<nSpecies;i2++)
                    {
                        a_col_sum[i1] = a_col_sum[i1] + Food_Web[i2][i1];
                        a_row_sum[i1] = a_row_sum[i1] + Food_Web[i1][i2];
                    }

                    if(a_row_sum[i1]==0 && a_col_sum[i1]==0)
                    {
                        Connected_Check==false;
                        break;
                    }
                }

                if(Connected_Check)
                {
                    std::vector<int> connected_component_id(num_vertices(Gen_Graph));

                    Test_Connected_Components = boost::connected_components(Gen_Graph, &connected_component_id[0]);
                }
            }
        }
    }

    std::vector<double> a_col_sum(nSpecies,0);

    for(int i1 = 0;i1<nSpecies;i1++)
        for(int i2 = 0;i2<nSpecies;i2++)
                a_col_sum[i1] = a_col_sum[i1] + Food_Web[i2][i1];

    for(int i1 = 0;i1<nSpecies;i1++)
        for(int i2 = 0;i2<nSpecies;i2++)
            if(Food_Web[i1][i2]==1)
                Food_Web[i1][i2]=1/a_col_sum[i2];

    return Food_Web;
}

//boost_matrix Niche_Model(int nSpecies, double Connectance)
//{
//    //default_random_engine generator;
//    //uniform_real_distribution<double> distribution(0,1);
//    //double something  = distribution(generator);

//    boost_matrix Food_Web;

//    std::vector< double > Trophic_Value(nSpecies,0);
//    std::vector< double > Generality(nSpecies,0);
//    std::vector< double > ci(nSpecies,0);

//    double beta  = (1/(2*Connectance))-1;

//    //default_random_engine generator ((unsigned int) time(NULL));

//    boost::random::random_device true_rand;

//    boost::random::mt19937 mt_generator(true_rand());

//    boost::random::uniform_real_distribution< double > uni_0_1(0,1);

//    for(int i=0;i<nSpecies;i++)
//    {
//        Trophic_Value[i] = uni_0_1(mt_generator);

//        Generality[i] = (1-pow((1-uni_0_1(mt_generator)),(1/beta)))*Trophic_Value[i];

//        boost::random::uniform_real_distribution< double > uni_ci((Generality[i]/2),Trophic_Value[i]);

//        ci[i] = uni_ci(mt_generator);
//    }

//    Food_Web.resize(boost::extents[nSpecies][nSpecies]);

//    //Base Web Generation

//    for(int i1 = 0;i1<nSpecies;i1++)
//    {
//        for(int i2 = 0;i2<nSpecies;i2++)
//        {
//            if(Trophic_Value[i2]>=(ci[i1]-(Generality[i1]/2))&&Trophic_Value[i2]<=(ci[i1]+(Generality[i1]/2))&&i1!=i2)
//                Food_Web[i1][i2] = 1;
//            else
//                Food_Web[i1][i2] = 0;
//        }
//    }

//    std::vector<double> a_row_sum(nSpecies,0);
//    std::vector<double> a_col_sum(nSpecies,0);

//    while(1)
//    {

//        for(int i1 = 0;i1<nSpecies;i1++)
//        {
//            a_col_sum[i1] = 0;
//            a_row_sum[i1] = 0;

//            for(int i2 = 0;i2<nSpecies;i2++)
//            {
//                    a_col_sum[i1] = a_col_sum[i1] + Food_Web[i2][i1];
//                    a_row_sum[i1] = a_row_sum[i1] + Food_Web[i1][i2];
//            }
//        }

//        //Check for disconnected Species
//        std::vector<std::vector<int>> Food_Chains = Calculate_Food_Chains(Food_Web);

//        std::vector<int> which_empty_species;

//        for(int i = 0;i<nSpecies;i++)
//        {
//            if(a_row_sum[i]==0 && a_col_sum[i]==0)
//            {
//                //Species is not connected to others
//                which_empty_species.push_back(i);
//            }
//            else if(Food_Chains[i].empty())
//            {
//                //Species is not connected to a basal resource
//                which_empty_species.push_back(i);
//            }
//        }

//        //Escape if none
//        if(which_empty_species.size()==0)
//            break;

//        //Reassign values if disconnected
//        for(int i=0;i<which_empty_species.size();i++)
//        {
//                Trophic_Value[which_empty_species[i]] = uni_0_1(mt_generator);

//                Generality[which_empty_species[i]] = (1-pow((1-uni_0_1(mt_generator)),(1/beta)))*Trophic_Value[which_empty_species[i]];

//                boost::random::uniform_real_distribution< double > uni_ci((Generality[which_empty_species[i]]/2),Trophic_Value[which_empty_species[i]]);

//                ci[which_empty_species[i]] = uni_ci(mt_generator);
//        }

//        //Rebuild Web
//        for(int i1 = 0;i1<nSpecies;i1++)
//        {
//            for(int i2 = 0;i2<nSpecies;i2++)
//            {
//                if(Trophic_Value[i2]>=(ci[i1]-(Generality[i1]/2))&&Trophic_Value[i2]<=(ci[i1]+(Generality[i1]/2))&&i1!=i2)
//                    Food_Web[i1][i2] = 1;
//                else
//                    Food_Web[i1][i2] = 0;
//            }
//        }
//    }

//    for(int i1 = 0;i1<nSpecies;i1++)
//    {
//        for(int i2 = 0;i2<nSpecies;i2++)
//        {
//            if(Food_Web[i1][i2]==1)
//            {
//                Food_Web[i1][i2]=1/a_col_sum[i2];
//            }
//        }
//    }

//    return Food_Web;
//}

Food_Web_Metadata::Food_Web_Metadata(const boost_matrix &Input_Food_Web):
    Food_Web(Input_Food_Web)
{
    int nSpecies = Food_Web.shape()[0];

    std::vector<double> a_row_sum(nSpecies,0);
    std::vector<double> a_col_sum(nSpecies,0);
    double a_total_sum = 0;

    for(int i = 0;i<nSpecies;i++)
    {
        for(int i2 = 0;i2<nSpecies;i2++)
        {
            a_row_sum[i] += Food_Web[i][i2];

            a_col_sum[i] += Food_Web[i2][i];
        }

        a_total_sum = a_total_sum + a_row_sum[i];
    }

    std::vector<double> H_in(nSpecies,0);
    std::vector<double> H_out(nSpecies,0);
    N_in.assign(nSpecies,0);
    N_out.assign(nSpecies,0);

    H = 0;
    AMI = 0;
    Hc = 0;

    for(int i_Sp1 = 0;i_Sp1<nSpecies;i_Sp1++)
    {
        for(int i_Sp2 = 0;i_Sp2<nSpecies;i_Sp2++)
        {

            if(Food_Web[i_Sp2][i_Sp1]>0)
            {
                if(a_col_sum[i_Sp1]>0)
                    H_in[i_Sp1] -= Food_Web[i_Sp2][i_Sp1]/a_col_sum[i_Sp1]*Log2(Food_Web[i_Sp2][i_Sp1]/a_col_sum[i_Sp1]);
            }

            if(Food_Web[i_Sp1][i_Sp2]>0)
            {
                if(a_row_sum[i_Sp1]>0)
                    H_out[i_Sp1] -= Food_Web[i_Sp1][i_Sp2]/a_row_sum[i_Sp1]*Log2(Food_Web[i_Sp1][i_Sp2]/a_row_sum[i_Sp1]);

                if(a_row_sum[i_Sp1]>0 && a_col_sum[i_Sp2]>0)
                {
                    AMI += Food_Web[i_Sp1][i_Sp2]/a_total_sum*Log2((Food_Web[i_Sp1][i_Sp2]*a_total_sum)/(a_row_sum[i_Sp1]*a_col_sum[i_Sp2]));

                    Hc -= Food_Web[i_Sp1][i_Sp2]/a_total_sum*Log2((Food_Web[i_Sp1][i_Sp2]*Food_Web[i_Sp1][i_Sp2])/(a_row_sum[i_Sp1]*a_col_sum[i_Sp2]));

                }

                H -= Food_Web[i_Sp1][i_Sp2]/a_total_sum*Log2(Food_Web[i_Sp1][i_Sp2]/a_total_sum);
            }
        }

        if(a_col_sum[i_Sp1]!=0)
            N_in[i_Sp1] = pow(2,H_in[i_Sp1]);

        if(a_row_sum[i_Sp1]!=0)
            N_out[i_Sp1] = pow(2,H_out[i_Sp1]);

        if(N_in[i_Sp1]>=1)
        {
            Generality.push_back(N_in[i_Sp1]-1);
            is_Producer.push_back(0);
            is_Consumer.push_back(1);
        }
        else
        {
            Generality.push_back(0);
            is_Producer.push_back(1);
            is_Consumer.push_back(0);
        }

        Positional_Index.push_back(N_in[i_Sp1]/(N_in[i_Sp1]+N_out[i_Sp1]));
    }

    //Trophic Bullshit

    Trophic_Positions = Calculate_Food_Chains(Food_Web);

    for(int i=0;i<nSpecies;i++)
    {
        Min_Trophic_Position.push_back(*std::min_element(Trophic_Positions[i].begin(),Trophic_Positions[i].end()));
        Max_Trophic_Position.push_back(*std::max_element(Trophic_Positions[i].begin(),Trophic_Positions[i].end()));

        if(Max_Trophic_Position[i]>nSpecies)
            Max_Trophic_Position[i]=(double)nSpecies;
    }

    Avg_Trophic_Position = Min_Trophic_Position;

    double Max_Min_Trophic_Position = *std::max_element(Min_Trophic_Position.begin(),Min_Trophic_Position.end());

    for(double i_Troph=1;i_Troph<Max_Min_Trophic_Position;i_Troph++)
    {
        for(int i_Sp=0;i_Sp<nSpecies;i_Sp++)
        {
            if(Min_Trophic_Position[i_Sp]==i_Troph)
            {
                std::vector<double> Avg_Troph;

                for(int i_Prey=0;i_Prey<nSpecies;i_Prey++)
                {
                    if(Food_Web[i_Prey][i_Sp]>0)
                        if(Min_Trophic_Position[i_Prey]==i_Troph)
                            Avg_Troph.push_back(Min_Trophic_Position[i_Prey]+1);
                        else
                            Avg_Troph.push_back(Avg_Trophic_Position[i_Prey]+1);
                }

                Avg_Trophic_Position[i_Sp] = std::accumulate(Avg_Troph.begin(),Avg_Troph.end(),0.0)/Avg_Troph.size();
            }
        }
    }


}

Food_Web_Metadata Calculate_Food_Web_Metadata (const boost_matrix &Food_Web)
{
    Food_Web_Metadata to_Return(Food_Web);

//    int nSpecies = Food_Web.shape()[0];

//    std::vector<double> a_row_sum(nSpecies,0);
//    std::vector<double> a_col_sum(nSpecies,0);
//    double a_total_sum = 0;

//    for(int i = 0;i<nSpecies;i++)
//    {
//        for(int i2 = 0;i2<nSpecies;i2++)
//        {
//            a_row_sum[i] += Food_Web[i][i2];

//            a_col_sum[i] += Food_Web[i2][i];
//        }

//        a_total_sum = a_total_sum + a_row_sum[i];
//    }

//    to_Return.H_in.assign(nSpecies,0);
//    to_Return.H_out.assign(nSpecies,0);
//    to_Return.N_in.assign(nSpecies,0);
//    to_Return.N_out.assign(nSpecies,0);

//    to_Return.H = 0;
//    to_Return.AMI = 0;
//    to_Return.Hc = 0;

//    for(int i_Sp1 = 0;i_Sp1<nSpecies;i_Sp1++)
//    {
//        for(int i_Sp2 = 0;i_Sp2<nSpecies;i_Sp2++)
//        {

//            if(Food_Web[i_Sp2][i_Sp1]>0)
//            {
//                if(a_row_sum[i_Sp1]>0)
//                    to_Return.H_in[i_Sp1] -= Food_Web[i_Sp2][i_Sp1]/a_col_sum[i_Sp1]*Log2(Food_Web[i_Sp2][i_Sp1]/a_col_sum[i_Sp1]);
//            }

//            if(Food_Web[i_Sp1][i_Sp2]>0)
//            {
//                if(a_col_sum[i_Sp1]>0)
//                    to_Return.H_out[i_Sp1] -= Food_Web[i_Sp1][i_Sp2]/a_row_sum[i_Sp1]*Log2(Food_Web[i_Sp1][i_Sp2]/a_row_sum[i_Sp2]);

//                if(a_row_sum[i_Sp1]>0 && a_col_sum[i_Sp2]>0)
//                {
//                    to_Return.AMI += Food_Web[i_Sp1][i_Sp2]/a_total_sum*Log2((Food_Web[i_Sp1][i_Sp2]*a_total_sum)/(a_row_sum[i_Sp1]*a_col_sum[i_Sp2]));

//                    to_Return.Hc -= Food_Web[i_Sp1][i_Sp2]/a_total_sum*Log2((Food_Web[i_Sp1][i_Sp2]*Food_Web[i_Sp1][i_Sp2])/(a_row_sum[i_Sp1]*a_col_sum[i_Sp2]));

//                }

//                to_Return.H -= Food_Web[i_Sp1][i_Sp2]/a_total_sum*Log2(Food_Web[i_Sp1][i_Sp2]/a_total_sum);
//            }
//        }

//        if(a_col_sum[i_Sp1]!=0)
//            to_Return.N_in[i_Sp1] = pow(2,to_Return.H_in[i_Sp1]);

//        if(a_row_sum[i_Sp1]!=0)
//            to_Return.N_out[i_Sp1] = pow(2,to_Return.H_out[i_Sp1]);

//        //to_Return.Trophic_Position.push_back((a_col_sum[i_Sp1]*to_Return.N_in[i_Sp1])/(a_col_sum[i_Sp1]*to_Return.N_in[i_Sp1]+a_row_sum[i_Sp1]*to_Return.N_out[i_Sp1]));

//        if(to_Return.N_in[i_Sp1]>=1)
//        {
//            to_Return.Generality.push_back(to_Return.N_in[i_Sp1]-1);
//            to_Return.is_Producer.push_back(0);
//            to_Return.is_Consumer.push_back(1);
//        }
//        else
//        {
//            to_Return.Generality.push_back(0);
//            to_Return.is_Producer.push_back(1);
//            to_Return.is_Consumer.push_back(0);
//        }
//    }

//    //Trophic Bullshit

//    std::vector<std::vector<int>> Food_Chains = Calculate_Food_Chains(Food_Web);

//    for(int i=0;i<nSpecies;i++)
//    {

//        to_Return.Trophic_Position.push_back(*std::max_element(Food_Chains[i].begin(),Food_Chains[i].end()));
//    }

    return to_Return;
}

std::vector<std::vector<int>> Calculate_Food_Chains(const boost_matrix &Food_Web)
{
    std::vector<std::vector<int>> Trophic_Paths(Food_Web.shape()[0],std::vector<int>());
//    std::vector<std::vector<int>> Prey_Species(Food_Web.shape()[0],std::vector<int>());

//std::vector<double> a_col_sum(Food_Web.shape()[0],0);

        for(int i1 = 0;i1<Food_Web.shape()[0];i1++)
        {
            bool is_Producer = true;

            for(int i2 = 0;i2<Food_Web.shape()[0];i2++)
            {
                if(i1!=i2&&Food_Web[i2][i1]>0)
                {
                    is_Producer = false;
                    break;
                }
            }

            if(is_Producer)
            {
                Trophic_Paths[i1].push_back(0);
            }
        }

    int i_Troph=0;

    bool Search_Check = true;

    while(Search_Check&&i_Troph<Food_Web.shape()[0]*2)
    {
        Search_Check = false;

        for(int i_Sp=0;i_Sp<Food_Web.shape()[0];i_Sp++)
            if(!Trophic_Paths[i_Sp].empty())
                if(*std::min_element(Trophic_Paths[i_Sp].begin(),Trophic_Paths[i_Sp].end())==i_Troph) //changed to min element
                    for(int i_Pred=0;i_Pred<Food_Web.shape()[0];i_Pred++)
                        if(Food_Web[i_Sp][i_Pred]>0)
                        {
                            Search_Check=true;

                            //Only add trophic level if predator's not already given that trophic level
                            if(std::find(Trophic_Paths[i_Pred].begin(),Trophic_Paths[i_Pred].end(),i_Troph+1)==Trophic_Paths[i_Pred].end())
                                Trophic_Paths[i_Pred].push_back(i_Troph+1);

                            //Only add prey species if predator's not already given that prey
//                            if(Prey_Species[i_Pred].empty())
//                                Prey_Specis[i_Pred].push_back(i_Sp);
//                            else if(std::find(Prey_Species[i_Pred].begin(),Prey_Species[i_Pred].end(),i_Sp)==Prey_Species[i_Pred].end())
//                                Prey_Species[i_Pred].push_back(i_Sp);
                        }
        i_Troph++;
    }

    return Trophic_Paths;
}
