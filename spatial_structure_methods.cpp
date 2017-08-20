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
#include "spatial_structure_methods.h"

#include <utility>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/connected_components.hpp>

#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>

#include <boost/graph/random.hpp>

#include <boost/graph/isomorphism.hpp>

#include <Eigen/Dense>

#include <boost/multi_array.hpp>

Spatial_Structure_Metadata::Spatial_Structure_Metadata()
{
    Spatial_Structure.resize(boost::extents[1][1]);
    Spatial_Structure[0][0] = 0;

    Average_Degree = 0;
    Degree_Skewness = 0;
    Average_Path_Length = 0;
    Clustering_Coefficient = 0;
    Eigenratio = 0;

    Degree.assign(1,0);
}

Spatial_Structure_Metadata::Spatial_Structure_Metadata(int nPatches)
{
    Spatial_Structure.resize(boost::extents[nPatches][nPatches]);

    std::fill(Spatial_Structure.data(),Spatial_Structure.data()+Spatial_Structure.num_elements(),0);

    Average_Degree = 0;
    Degree_Skewness = 0;
    Average_Path_Length = 0;
    Clustering_Coefficient = 0;
    Eigenratio = 0;

    Degree.assign(1,0);
}

void Spatial_Structure_Metadata::operator()(Spatial_Structure_Metadata New_Values)
{
    Average_Degree = New_Values.Average_Degree;
    Degree_Skewness = New_Values.Degree_Skewness;
    Average_Path_Length = New_Values.Average_Path_Length;
    Clustering_Coefficient = New_Values.Clustering_Coefficient;
    Eigenratio = New_Values.Eigenratio;
    Degree = New_Values.Degree;

    Spatial_Structure.resize(boost::extents[New_Values.Spatial_Structure.shape()[0]][New_Values.Spatial_Structure.shape()[1]]);

    for(int i=0;i<New_Values.Spatial_Structure.shape()[0];i++)
        for(int i2=0;i2<New_Values.Spatial_Structure.shape()[1];i2++)
            Spatial_Structure[i][i2]=New_Values.Spatial_Structure[i][i2];
}

bool Spatial_Structure_Metadata::is_Equal(const Spatial_Structure_Metadata &Compare_To)
{
    if(Average_Degree==Compare_To.Average_Degree)
        if(Degree_Skewness==Compare_To.Degree_Skewness)
            if(Average_Path_Length==Compare_To.Average_Path_Length)
                if(Clustering_Coefficient==Compare_To.Clustering_Coefficient)
                    if(Eigenratio = Compare_To.Eigenratio)
                        if(is_Isomorphic(Spatial_Structure,Compare_To.Spatial_Structure))
                            return true;

    return false;
}

Quotient_Structure_Metadata::Quotient_Structure_Metadata()
{
    Quotient_Structure.resize(boost::extents[1][1]);
    Quotient_Structure[0][0] = 0;

    Average_In_Degree = 0;
    Average_Out_Degree = 0;
    In_Degree_Skewness = 0;
    Out_Degree_Skewness = 0;

    Average_Path_Length = 0;
    Clustering_Coefficient = 0;
    Eigenratio = 0;

    In_Degree.assign(1,0);
    Out_Degree.assign(1,0);
}

// type for weight/distance on each edge
// Changing this WILL break things until code is cleaned up
typedef double t_weight;

Spatial_Structure_Metadata Erdos_Renyi_Undirected(int nPatches, double Connectance)
{
    Spatial_Structure_Metadata to_Return;

    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    //Generate Random Graph
    typedef boost::sorted_erdos_renyi_iterator<boost::random::mt19937,Graph> ERGen;

    boost::random::random_device true_rand;

    boost::random::mt19937 mt_generator(true_rand());

    int Test_Connected_Components = 0;

    Graph Gen_Graph;

    while(Test_Connected_Components != 1)
    {

    Gen_Graph = Graph(ERGen(mt_generator,nPatches,Connectance),ERGen(),nPatches);

    std::vector<int> connected_component_id(num_vertices(Gen_Graph));

    Test_Connected_Components = boost::connected_components(Gen_Graph, &connected_component_id[0]);

    }


    //Setting up vertex & edge weight access

    typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    Index_Map index = get(boost::vertex_index, Gen_Graph);

    typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Translate into boost_matrix

    //base matrix filled with 0's
    to_Return.Spatial_Structure.resize(boost::extents[nPatches][nPatches]);
    std::fill(to_Return.Spatial_Structure.origin(),to_Return.Spatial_Structure.origin()+to_Return.Spatial_Structure.num_elements(),0);

    //iterate over all edges, add to matrix
    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            //Weight assignments, all weights = 1
            put(Gen_Graph_Weight,*ei,1);

            to_Return.Spatial_Structure[index[source(*ei,Gen_Graph)]][index[target(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
            to_Return.Spatial_Structure[index[target(*ei,Gen_Graph)]][index[source(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
    {
        double Emigration = 0;

        for(int j=0;j<nPatches;j++)
        {
            Emigration+=to_Return.Spatial_Structure[i][j];
        }

        to_Return.Spatial_Structure[i][i]=-Emigration;
    }

    //Structure Metrics

    //Degree, Average, and Skewness
    to_Return.Average_Degree = 0;
    to_Return.Degree.clear();

    {boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (boost::tie(vi,vi_end) = vertices(Gen_Graph); vi != vi_end; ++vi)
        {
            to_Return.Average_Degree += (double)out_degree(*vi,Gen_Graph);
            to_Return.Degree.push_back((double)out_degree(*vi,Gen_Graph));
        }
    }

    to_Return.Average_Degree = to_Return.Average_Degree/nPatches;

    double Degree_Second_Cumulant = 0;
    double Degree_Third_Cumulant = 0;
    for(int i=0;i<nPatches;i++)
    {
        Degree_Second_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,2);
        Degree_Third_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,3);
    }

    if(Degree_Second_Cumulant==0||Degree_Third_Cumulant==0)
        to_Return.Degree_Skewness = 0;
    else
        to_Return.Degree_Skewness = Degree_Third_Cumulant/pow(Degree_Second_Cumulant,3/2);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<nPatches-1;i++)
            for(int i2=i+1;i2<nPatches;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/((double)nPatches*((double)nPatches-1)))*((double)Paths_Sum);

    //Eigenratio

    //Setup/Transfer to Eigen Matrix

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

    eigen_matrix Eigen_Mat = eigen_matrix::Constant(nPatches,nPatches,0);

    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            Eigen_Mat(index[source(*ei,Gen_Graph)],index[target(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
            Eigen_Mat(index[target(*ei,Gen_Graph)],index[source(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
        Eigen_Mat(i,i) = to_Return.Spatial_Structure[i][i];

    Eigen::SelfAdjointEigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

    to_Return.Eigenratio = eigensolver.eigenvalues()(0,0)/eigensolver.eigenvalues()(eigensolver.eigenvalues().rows()-2,0);

    return to_Return;
}

Spatial_Structure_Metadata Dendritic_Undirected(int nPatches, double Branching_Probability)
{
    Spatial_Structure_Metadata to_Return;

    if(nPatches>0)
    {

    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    boost::random::random_device true_rand;

    boost::random::mt19937 mt_generator(true_rand());

    boost::random::uniform_real_distribution< double > uni_0_1(0,1);

    Graph Gen_Graph;

    int iPatches = 0;

    std::vector<int> Growing_Nodes(1,0);

    while(Growing_Nodes.size() > 0)
    {
        boost::random::uniform_int_distribution<> Rand_Grow(0,Growing_Nodes.size()-1);

        int Choose_Node = Rand_Grow(mt_generator);

        int Active_Node = Growing_Nodes[Choose_Node];
        Growing_Nodes.erase(Growing_Nodes.begin()+Choose_Node);

        if(uni_0_1(mt_generator)<Branching_Probability)
        {
            add_edge(Active_Node,++iPatches,Gen_Graph);
            Growing_Nodes.push_back(iPatches);

            if(iPatches==(nPatches-1))
                break;
        }

        add_edge(Active_Node,++iPatches,Gen_Graph);
        Growing_Nodes.push_back(iPatches);

        if(iPatches==(nPatches-1))
            break;
    }

    //Setting up vertex & edge weight access

    typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    Index_Map index = get(boost::vertex_index, Gen_Graph);

    typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Translate into boost_matrix

    //base matrix filled with 0's
    to_Return.Spatial_Structure.resize(boost::extents[nPatches][nPatches]);
    std::fill(to_Return.Spatial_Structure.origin(),to_Return.Spatial_Structure.origin()+to_Return.Spatial_Structure.num_elements(),0);

    //iterate over all edges, add to matrix
    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            //Weight assignments, all weights = 1
            put(Gen_Graph_Weight,*ei,1);

            to_Return.Spatial_Structure[index[source(*ei,Gen_Graph)]][index[target(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
            to_Return.Spatial_Structure[index[target(*ei,Gen_Graph)]][index[source(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
    {
        double Emigration = 0;

        for(int j=0;j<nPatches;j++)
        {
            Emigration+=to_Return.Spatial_Structure[i][j];
        }

        to_Return.Spatial_Structure[i][i]=-Emigration;
    }

    //Structure Metrics

    //Degree, Average, and Skewness
    to_Return.Average_Degree = 0;
    to_Return.Degree.clear();

    {boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (boost::tie(vi,vi_end) = vertices(Gen_Graph); vi != vi_end; ++vi)
        {
            to_Return.Average_Degree += (double)out_degree(*vi,Gen_Graph);
            to_Return.Degree.push_back((double)out_degree(*vi,Gen_Graph));
        }
    }

    to_Return.Average_Degree = to_Return.Average_Degree/nPatches;

    double Degree_Second_Cumulant = 0;
    double Degree_Third_Cumulant = 0;
    for(int i=0;i<nPatches;i++)
    {
        Degree_Second_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,2);
        Degree_Third_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,3);
    }

    if(Degree_Second_Cumulant==0||Degree_Third_Cumulant==0)
        to_Return.Degree_Skewness = 0;
    else
        to_Return.Degree_Skewness = Degree_Third_Cumulant/pow(Degree_Second_Cumulant,3/2);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<nPatches-1;i++)
            for(int i2=i+1;i2<nPatches;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/((double)nPatches*((double)nPatches-1)))*((double)Paths_Sum);

    //Eigenratio

    //Setup/Transfer to Eigen Matrix

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

    eigen_matrix Eigen_Mat = eigen_matrix::Constant(nPatches,nPatches,0);

    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            Eigen_Mat(index[source(*ei,Gen_Graph)],index[target(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
            Eigen_Mat(index[target(*ei,Gen_Graph)],index[source(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
        Eigen_Mat(i,i) = to_Return.Spatial_Structure[i][i];

    Eigen::SelfAdjointEigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

    to_Return.Eigenratio = eigensolver.eigenvalues()(0,0)/eigensolver.eigenvalues()(eigensolver.eigenvalues().rows()-2,0);

    }

    return to_Return;
}

Spatial_Structure_Metadata Radial_Tree_Undirected(int nPatches, int nBranches)
{
    Spatial_Structure_Metadata to_Return;

    if(nPatches>0)
    {

    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    Graph Gen_Graph;

    int iPatches = 0;

    std::vector<int> Growing_Nodes(1,0);

    while(Growing_Nodes.size() > 0&&iPatches<(nPatches-1))
    {
        int Active_Node = Growing_Nodes[0];
        Growing_Nodes.erase(Growing_Nodes.begin());

        for(int i=0;i<nBranches;i++)
        {
            add_edge(Active_Node,++iPatches,Gen_Graph);
            Growing_Nodes.push_back(iPatches);

            if(iPatches==(nPatches-1))
                break;
        }
    }

    //Setting up vertex & edge weight access

    typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    Index_Map index = get(boost::vertex_index, Gen_Graph);

    typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Translate into boost_matrix

    //base matrix filled with 0's
    to_Return.Spatial_Structure.resize(boost::extents[nPatches][nPatches]);
    std::fill(to_Return.Spatial_Structure.origin(),to_Return.Spatial_Structure.origin()+to_Return.Spatial_Structure.num_elements(),0);

    //iterate over all edges, add to matrix
    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            //Weight assignments, all weights = 1
            put(Gen_Graph_Weight,*ei,1);

            to_Return.Spatial_Structure[index[source(*ei,Gen_Graph)]][index[target(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
            to_Return.Spatial_Structure[index[target(*ei,Gen_Graph)]][index[source(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
    {
        double Emigration = 0;

        for(int j=0;j<nPatches;j++)
        {
            Emigration+=to_Return.Spatial_Structure[i][j];
        }

        to_Return.Spatial_Structure[i][i]=-Emigration;
    }

    //Structure Metrics

    //Degree, Average, and Skewness
    to_Return.Average_Degree = 0;
    to_Return.Degree.clear();

    {boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (boost::tie(vi,vi_end) = vertices(Gen_Graph); vi != vi_end; ++vi)
        {
            to_Return.Average_Degree += (double)out_degree(*vi,Gen_Graph);
            to_Return.Degree.push_back((double)out_degree(*vi,Gen_Graph));
        }
    }

    to_Return.Average_Degree = to_Return.Average_Degree/nPatches;

    double Degree_Second_Cumulant = 0;
    double Degree_Third_Cumulant = 0;
    for(int i=0;i<nPatches;i++)
    {
        Degree_Second_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,2);
        Degree_Third_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,3);
    }

    if(Degree_Second_Cumulant==0||Degree_Third_Cumulant==0)
        to_Return.Degree_Skewness = 0;
    else
        to_Return.Degree_Skewness = Degree_Third_Cumulant/pow(Degree_Second_Cumulant,3/2);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<nPatches-1;i++)
            for(int i2=i+1;i2<nPatches;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/((double)nPatches*((double)nPatches-1)))*((double)Paths_Sum);

    //Eigenratio

    //Setup/Transfer to Eigen Matrix

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

    eigen_matrix Eigen_Mat = eigen_matrix::Constant(nPatches,nPatches,0);

    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            Eigen_Mat(index[source(*ei,Gen_Graph)],index[target(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
            Eigen_Mat(index[target(*ei,Gen_Graph)],index[source(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
        Eigen_Mat(i,i) = to_Return.Spatial_Structure[i][i];

    Eigen::SelfAdjointEigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

    to_Return.Eigenratio = eigensolver.eigenvalues()(0,0)/eigensolver.eigenvalues()(eigensolver.eigenvalues().rows()-2,0);

    }

    return to_Return;
}

Spatial_Structure_Metadata Ring_Lattice_Undirected(int nPatches, int Neighbor_Distance)
{
    Spatial_Structure_Metadata to_Return;

    if(nPatches>0)
    {

    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    Graph Gen_Graph;

    for(int i=0;i<nPatches;i++)
    {
        for(int i2=(i+1);i2<=(i+Neighbor_Distance);i2++)
        {
            if(i2<nPatches)
                add_edge(i,i2,Gen_Graph);
            else
                add_edge(i,(i2-nPatches),Gen_Graph);
        }
    }

    //Setting up vertex & edge weight access

    typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    Index_Map index = get(boost::vertex_index, Gen_Graph);

    typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Translate into boost_matrix

    //base matrix filled with 0's
    to_Return.Spatial_Structure.resize(boost::extents[nPatches][nPatches]);
    std::fill(to_Return.Spatial_Structure.origin(),to_Return.Spatial_Structure.origin()+to_Return.Spatial_Structure.num_elements(),0);

    //iterate over all edges, add to matrix
    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            //Weight assignments, all weights = 1
            put(Gen_Graph_Weight,*ei,1);

            to_Return.Spatial_Structure[index[source(*ei,Gen_Graph)]][index[target(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
            to_Return.Spatial_Structure[index[target(*ei,Gen_Graph)]][index[source(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
    {
        double Emigration = 0;

        for(int j=0;j<nPatches;j++)
        {
            Emigration+=to_Return.Spatial_Structure[i][j];
        }

        to_Return.Spatial_Structure[i][i]=-Emigration;
    }

    //Structure Metrics

    //Degree, Average, and Skewness
    to_Return.Average_Degree = 0;
    to_Return.Degree.clear();

    {boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (boost::tie(vi,vi_end) = vertices(Gen_Graph); vi != vi_end; ++vi)
        {
            to_Return.Average_Degree += (double)out_degree(*vi,Gen_Graph);
            to_Return.Degree.push_back((double)out_degree(*vi,Gen_Graph));
        }
    }

    to_Return.Average_Degree = to_Return.Average_Degree/nPatches;

    double Degree_Second_Cumulant = 0;
    double Degree_Third_Cumulant = 0;
    for(int i=0;i<nPatches;i++)
    {
        Degree_Second_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,2);
        Degree_Third_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,3);
    }

    if(Degree_Second_Cumulant==0||Degree_Third_Cumulant==0)
        to_Return.Degree_Skewness = 0;
    else
        to_Return.Degree_Skewness = Degree_Third_Cumulant/pow(Degree_Second_Cumulant,3/2);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<nPatches-1;i++)
            for(int i2=i+1;i2<nPatches;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/((double)nPatches*((double)nPatches-1)))*((double)Paths_Sum);

    //Eigenratio

    //Setup/Transfer to Eigen Matrix

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

    eigen_matrix Eigen_Mat = eigen_matrix::Constant(nPatches,nPatches,0);

    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            Eigen_Mat(index[source(*ei,Gen_Graph)],index[target(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
            Eigen_Mat(index[target(*ei,Gen_Graph)],index[source(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
        Eigen_Mat(i,i) = to_Return.Spatial_Structure[i][i];

    Eigen::SelfAdjointEigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

    to_Return.Eigenratio = eigensolver.eigenvalues()(0,0)/eigensolver.eigenvalues()(eigensolver.eigenvalues().rows()-2,0);

    }

    return to_Return;
}

Spatial_Structure_Metadata Calculate_Unweighted_Symmetrical_Spatial_Structure_Metadata(boost_matrix Spatial_Structure)
{
    Spatial_Structure_Metadata to_Return;

    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    int nPatches = Spatial_Structure.shape()[0];

    Graph Gen_Graph(nPatches);

    for(int i=0;i<nPatches;i++)
        for(int i2=0;i2<nPatches;i2++)
            if(Spatial_Structure[i][i2]>0)
                add_edge(i,i2,Gen_Graph);

    //Setting up vertex & edge weight access

    typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    Index_Map index = get(boost::vertex_index, Gen_Graph);

    typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Translate into boost_matrix

    //base matrix filled with 0's
    to_Return.Spatial_Structure.resize(boost::extents[nPatches][nPatches]);
    std::fill(to_Return.Spatial_Structure.origin(),to_Return.Spatial_Structure.origin()+to_Return.Spatial_Structure.num_elements(),0);

    //iterate over all edges, add to matrix
    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            //Weight assignments, all weights = 1
            put(Gen_Graph_Weight,*ei,1);

            to_Return.Spatial_Structure[index[source(*ei,Gen_Graph)]][index[target(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
            to_Return.Spatial_Structure[index[target(*ei,Gen_Graph)]][index[source(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
    {
        double Emigration = 0;

        for(int j=0;j<nPatches;j++)
        {
            Emigration+=to_Return.Spatial_Structure[i][j];
        }

        to_Return.Spatial_Structure[i][i]=-Emigration;
    }

    //Structure Metrics

    //Degree, Average, and Skewness
    to_Return.Average_Degree = 0;
    to_Return.Degree.clear();

    {boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (boost::tie(vi,vi_end) = vertices(Gen_Graph); vi != vi_end; ++vi)
        {
            to_Return.Average_Degree += (double)out_degree(*vi,Gen_Graph);
            to_Return.Degree.push_back((double)out_degree(*vi,Gen_Graph));
        }
    }

    to_Return.Average_Degree = to_Return.Average_Degree/nPatches;

    double Degree_Second_Cumulant = 0;
    double Degree_Third_Cumulant = 0;
    for(int i=0;i<nPatches;i++)
    {
        Degree_Second_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,2);
        Degree_Third_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,3);
    }

    if(Degree_Second_Cumulant==0||Degree_Third_Cumulant==0)
        to_Return.Degree_Skewness = 0;
    else
        to_Return.Degree_Skewness = Degree_Third_Cumulant/pow(Degree_Second_Cumulant,3/2);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<nPatches-1;i++)
            for(int i2=i+1;i2<nPatches;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/(double)nPatches*((double)nPatches-1))*((double)Paths_Sum);

    //Eigenratio

    if(nPatches>1)
    {

    //Setup/Transfer to Eigen Matrix

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

    eigen_matrix Eigen_Mat = eigen_matrix::Constant(nPatches,nPatches,0);

    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            Eigen_Mat(index[source(*ei,Gen_Graph)],index[target(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
            Eigen_Mat(index[target(*ei,Gen_Graph)],index[source(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<nPatches;i++)
        Eigen_Mat(i,i) = to_Return.Spatial_Structure[i][i];

    Eigen::SelfAdjointEigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

    to_Return.Eigenratio = eigensolver.eigenvalues()(0,0)/eigensolver.eigenvalues()(eigensolver.eigenvalues().rows()-2,0);
    }
    else
    {
        to_Return.Eigenratio = 0;
    }

    return to_Return;
}

bool is_Isomorphic(boost_matrix Structure_A, boost_matrix Structure_B)
{
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::directedS,boost::no_property,boost::no_property> Graph;

    if(Structure_A.shape()[0]!=Structure_B.shape()[0])
        return false;

    int nVertices = Structure_A.shape()[0];

    Graph Graph_A(nVertices);

    for(int i=0;i<nVertices;i++)
        for(int i2=0;i2<nVertices;i2++)
            if(Structure_A[i][i2]>0)
                add_edge(i,i2,Graph_A);

    Graph Graph_B(nVertices);

    for(int i=0;i<nVertices;i++)
        for(int i2=0;i2<nVertices;i2++)
            if(Structure_B[i][i2]>0)
                add_edge(i,i2,Graph_B);

    //std::vector<boost::graph_traits<Graph>::vertex_descriptor> v1(Structure_A.shape()[0]);
    boost::property_map<Graph,boost::vertex_index_t>::type v1_index_map = get(boost::vertex_index,Graph_A);

    std::vector<boost::graph_traits<Graph>::vertex_descriptor> f(Structure_A.shape()[0]);

    return isomorphism(Graph_A,Graph_B, boost::isomorphism_map(boost::make_iterator_property_map(f.begin(),v1_index_map,f[0])));
}

Spatial_Structure_Metadata Permuted_Shipley_Skinner(std::vector<int> Remove_Which)
{
    Spatial_Structure_Metadata to_Return;

    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    Graph Gen_Graph;

    t_weight Distance_Decay_Coeff = 1;

    std::vector<t_weight> Distance_Weight;

    std::vector<std::pair<int,int>> Removeable_Edges;

    Removeable_Edges.push_back(std::pair<int,int>(0,1));
    Removeable_Edges.push_back(std::pair<int,int>(0,2));
    Removeable_Edges.push_back(std::pair<int,int>(1,2));
    Removeable_Edges.push_back(std::pair<int,int>(3,13));
    Removeable_Edges.push_back(std::pair<int,int>(4,7));
    Removeable_Edges.push_back(std::pair<int,int>(6,17));
    Removeable_Edges.push_back(std::pair<int,int>(7,11));
    Removeable_Edges.push_back(std::pair<int,int>(7,15));
    Removeable_Edges.push_back(std::pair<int,int>(8,12));
    Removeable_Edges.push_back(std::pair<int,int>(8,14));
    Removeable_Edges.push_back(std::pair<int,int>(10,17));
    Removeable_Edges.push_back(std::pair<int,int>(13,15));
    Removeable_Edges.push_back(std::pair<int,int>(14,15));
    Removeable_Edges.push_back(std::pair<int,int>(14,17));
    Removeable_Edges.push_back(std::pair<int,int>(16,17));

    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*1.48));//e_0
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*6.12));//e_1
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*1.12));//e_2
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*4.06));//e_3
    Distance_Weight.push_back(.152*exp(-Distance_Decay_Coeff*0));//e_4
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*9.27));//e_5
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*10.4));//e_6
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*4.11));//e_7
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*10.23));//e_8
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*5.59));//e_9
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*7.24));//e_10
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*7.09));//e_11
    Distance_Weight.push_back(.484*exp(-Distance_Decay_Coeff*0));//e_12
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*5.64));//e_13
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*.72));//e_14
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*3.57));//e_15
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*3.46));//e_16
    Distance_Weight.push_back(.392*exp(-Distance_Decay_Coeff*0));//e_17
    Distance_Weight.push_back(.555*exp(-Distance_Decay_Coeff*0));//e_18
    Distance_Weight.push_back(.830*exp(-Distance_Decay_Coeff*0));//e_19
    Distance_Weight.push_back(.3*exp(-Distance_Decay_Coeff*0));//e_20
    Distance_Weight.push_back(.445*exp(-Distance_Decay_Coeff*0));//e_21
    Distance_Weight.push_back(.222*exp(-Distance_Decay_Coeff*0));//e_22
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*1.53));//e_23
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*3.48));//e_24
    Distance_Weight.push_back(.175*exp(-Distance_Decay_Coeff*0));//e_25
    Distance_Weight.push_back(.1*exp(-Distance_Decay_Coeff*3.43));//e_26

    //Links that can be broken
    add_edge(0,1,Distance_Weight[0],Gen_Graph); //e_0
    add_edge(0,2,Distance_Weight[1],Gen_Graph); //e_1

    add_edge(1,2,Distance_Weight[2],Gen_Graph); //e_2

    add_edge(3,13,Distance_Weight[5],Gen_Graph); //e_5

    add_edge(4,7,Distance_Weight[6],Gen_Graph); //e_6

    add_edge(6,17,Distance_Weight[10],Gen_Graph); //e_10

    add_edge(7,11,Distance_Weight[11],Gen_Graph); //e_11
    add_edge(7,15,Distance_Weight[13],Gen_Graph); //e_13

    add_edge(8,12,Distance_Weight[14],Gen_Graph); //e_14
    add_edge(8,14,Distance_Weight[15],Gen_Graph); //e_15

    add_edge(10,17,Distance_Weight[22],Gen_Graph); //e_22

    add_edge(13,15,Distance_Weight[23],Gen_Graph); //e_23

    add_edge(14,15,Distance_Weight[24],Gen_Graph); //e_24
    add_edge(14,17,Distance_Weight[25],Gen_Graph); //e_25

    add_edge(16,17,Distance_Weight[26],Gen_Graph); //e_26

//    boost::random::random_device true_rand;
//    boost::random::mt19937 mt_generator(true_rand());

//    for(int i=0;i<nRewires;i++)
//        remove_edge(random_edge(Gen_Graph,mt_generator),Gen_Graph);

    for(int i=0;i<Remove_Which.size();i++)
        if(Remove_Which[i]<Removeable_Edges.size())
            remove_edge(Removeable_Edges[Remove_Which[i]].first,Removeable_Edges[Remove_Which[i]].second,Gen_Graph);

    //Guarunteed Links

    add_edge(1,5,Distance_Weight[3],Gen_Graph); //e_3

    add_edge(2,11,Distance_Weight[4],Gen_Graph); //e_4

    add_edge(4,11,Distance_Weight[7],Gen_Graph); //e_7
    add_edge(4,12,Distance_Weight[8],Gen_Graph); //e_8

    add_edge(5,6,Distance_Weight[9],Gen_Graph); //e_9

    add_edge(7,13,Distance_Weight[12],Gen_Graph); //e_12

    add_edge(8,17,Distance_Weight[16],Gen_Graph); //e_16

    add_edge(9,13,Distance_Weight[17],Gen_Graph); //e_17
    add_edge(9,14,Distance_Weight[18],Gen_Graph); //e_18
    add_edge(9,15,Distance_Weight[19],Gen_Graph); //e_19
    add_edge(9,17,Distance_Weight[20],Gen_Graph); //e_20

    add_edge(10,16,Distance_Weight[21],Gen_Graph); //e_21

    //Setting up vertex & edge weight access

    typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    Index_Map index = get(boost::vertex_index, Gen_Graph);

    typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Translate into boost_matrix

    //base matrix filled with 0's
    to_Return.Spatial_Structure.resize(boost::extents[18][18]);
    std::fill(to_Return.Spatial_Structure.origin(),to_Return.Spatial_Structure.origin()+to_Return.Spatial_Structure.num_elements(),0);

    //iterate over all edges, add to matrix
    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            //Weight assignments, all weights = 1
            //put(Gen_Graph_Weight,*ei,1);

            to_Return.Spatial_Structure[index[source(*ei,Gen_Graph)]][index[target(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
            to_Return.Spatial_Structure[index[target(*ei,Gen_Graph)]][index[source(*ei,Gen_Graph)]] = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<18;i++)
    {
        double Emigration = 0;

        for(int j=0;j<18;j++)
        {
            Emigration+=to_Return.Spatial_Structure[i][j];
        }

        to_Return.Spatial_Structure[i][i]=-Emigration;
    }

    //Structure Metrics

    //Degree, Average, and Skewness
    to_Return.Average_Degree = 0;
    to_Return.Degree.clear();

    {boost::graph_traits<Graph>::vertex_iterator vi,vi_end;
        for (boost::tie(vi,vi_end) = vertices(Gen_Graph); vi != vi_end; ++vi)
        {
            to_Return.Average_Degree += (double)out_degree(*vi,Gen_Graph);
            to_Return.Degree.push_back((double)out_degree(*vi,Gen_Graph));
        }
    }

    to_Return.Average_Degree = to_Return.Average_Degree/18;

    double Degree_Second_Cumulant = 0;
    double Degree_Third_Cumulant = 0;
    for(int i=0;i<18;i++)
    {
        Degree_Second_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,2);
        Degree_Third_Cumulant += pow(to_Return.Degree[i]-to_Return.Average_Degree,3);
    }

    to_Return.Degree_Skewness = Degree_Third_Cumulant/pow(Degree_Second_Cumulant,3/2);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<18-1;i++)
            for(int i2=i+1;i2<18;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/((double)18*((double)18-1)))*((double)Paths_Sum);

    //Eigenratio

    //Setup/Transfer to Eigen Matrix

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

    eigen_matrix Eigen_Mat = eigen_matrix::Constant(18,18,0);

    {boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(Gen_Graph); ei != ei_end; ++ei)
        {
            Eigen_Mat(index[source(*ei,Gen_Graph)],index[target(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
            Eigen_Mat(index[target(*ei,Gen_Graph)],index[source(*ei,Gen_Graph)]) = get(Gen_Graph_Weight,*ei);
        }
    }

    for(int i=0;i<18;i++)
        Eigen_Mat(i,i) = to_Return.Spatial_Structure[i][i];


    Eigen::SelfAdjointEigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

    to_Return.Eigenratio = eigensolver.eigenvalues()(0,0)/eigensolver.eigenvalues()(eigensolver.eigenvalues().rows()-2,0);

    return to_Return;
}

Quotient_Structure_Metadata Calculate_Quotient_Structure_Metadata(const boost_matrix &Spatial_Structure, std::vector<int> Synchronized_Clusters)
{
    Quotient_Structure_Metadata to_Return;

    //Calculate Quotient Structure

    std::vector<std::vector<double>> temp_Quotient_Structure;

    std::vector<int> Clusters_Found;

    std::vector<std::vector<int>> Cluster_Index;

    Clusters_Found.push_back(Synchronized_Clusters[0]);

    Cluster_Index.push_back(std::vector<int>(1,0));

    temp_Quotient_Structure.push_back(std::vector<double>());

    for(int i=0;i<Spatial_Structure.shape()[1];i++)
        temp_Quotient_Structure[0].push_back(Spatial_Structure[i][0]);

    for(int i=1;i<Synchronized_Clusters.size();i++)
    {
        bool Cluster_Found_Match = false;

        for(int j=0;j<Clusters_Found.size();j++)
        {
            if(Synchronized_Clusters[i]==Clusters_Found[j])
            {
                Cluster_Index[j].push_back(i);

                for(int k=0;k<Spatial_Structure.shape()[1];k++)
                    temp_Quotient_Structure[j][k]+=Spatial_Structure[k][i];

                Cluster_Found_Match = true;
                break;
            }
        }

        if(!Cluster_Found_Match)
        {
            Clusters_Found.push_back(Synchronized_Clusters[i]);
            Cluster_Index.push_back(std::vector<int>(1,i));
            temp_Quotient_Structure.push_back(std::vector<double>());

            for(int j=0;j<Spatial_Structure.shape()[1];j++)
                temp_Quotient_Structure.back().push_back(Spatial_Structure[j][i]);
        }
    }

   to_Return.Quotient_Structure.resize(boost::extents[temp_Quotient_Structure.size()][temp_Quotient_Structure.size()]);

   for(int i=0;i<temp_Quotient_Structure.size();i++)
   {
       for(int j=0;j<Cluster_Index.size();j++)
           to_Return.Quotient_Structure[j][i] = temp_Quotient_Structure[i][Cluster_Index[j][0]];
   }

   int nPatches = to_Return.Quotient_Structure.shape()[0];

   if(nPatches<2) //Escapes if quotient network is one patch
       return to_Return;

   for(int i=0;i<nPatches;i++)
   {
       double Emigration = 0;

       for(int j=0;j<nPatches;j++)
       {
           if(i!=j)
                Emigration+=to_Return.Quotient_Structure[i][j];
       }

       to_Return.Quotient_Structure[i][i]=-Emigration;
   }


   //Degree, Average, and Skewness
   to_Return.Average_In_Degree = 0;
   to_Return.Average_Out_Degree = 0;

   to_Return.In_Degree.clear();
   to_Return.Out_Degree.clear();

   for(int i=0;i<nPatches;i++)
   {
       to_Return.In_Degree.push_back(0);
       to_Return.Out_Degree.push_back(0);

       for(int j=0;j<nPatches;j++)
       {
           if(j!=i)
           {
                to_Return.In_Degree[i]+=to_Return.Quotient_Structure[i][j];
                to_Return.Out_Degree[i]+=to_Return.Quotient_Structure[j][i];
           }
       }

       to_Return.Average_In_Degree += to_Return.In_Degree[i];
       to_Return.Average_Out_Degree += to_Return.Out_Degree[i];
   }

   to_Return.Average_In_Degree = to_Return.Average_In_Degree/nPatches;
   to_Return.Average_Out_Degree = to_Return.Average_Out_Degree/nPatches;

   double In_Degree_Second_Cumulant = 0;
   double In_Degree_Third_Cumulant = 0;

   double Out_Degree_Second_Cumulant = 0;
   double Out_Degree_Third_Cumulant = 0;
   for(int i=0;i<nPatches;i++)
   {
       In_Degree_Second_Cumulant += pow(to_Return.In_Degree[i]-to_Return.Average_In_Degree,2);
       In_Degree_Third_Cumulant += pow(to_Return.In_Degree[i]-to_Return.Average_In_Degree,3);

       Out_Degree_Second_Cumulant += pow(to_Return.Out_Degree[i]-to_Return.Average_Out_Degree,2);
       Out_Degree_Third_Cumulant += pow(to_Return.Out_Degree[i]-to_Return.Average_Out_Degree,3);
   }

   if(In_Degree_Second_Cumulant==0||In_Degree_Third_Cumulant==0)
       to_Return.In_Degree_Skewness = 0;
   else
       to_Return.In_Degree_Skewness = In_Degree_Third_Cumulant/pow(In_Degree_Second_Cumulant,3/2);

   if(Out_Degree_Second_Cumulant==0||Out_Degree_Third_Cumulant==0)
       to_Return.Out_Degree_Skewness = 0;
   else
       to_Return.Out_Degree_Skewness = Out_Degree_Third_Cumulant/pow(Out_Degree_Second_Cumulant,3/2);

   //Eigenratio

   //Setup/Transfer to Eigen Matrix

   if(nPatches>2)
   {

       typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

       eigen_matrix Eigen_Mat = eigen_matrix::Constant(nPatches,nPatches,0);

       for(int i=0;i<nPatches;i++)
           for(int j=0;j<nPatches;j++)
           {
               Eigen_Mat(i,j) = to_Return.Quotient_Structure[i][j];
           }

       Eigen::EigenSolver<eigen_matrix> eigensolver(Eigen_Mat);

       double min_eigen = 0;
       double max_eigen = 0;

       for(int i=0;i<eigensolver.eigenvalues().rows();i++)
       {
           if(eigensolver.eigenvalues()(i,0).real()<min_eigen||min_eigen==0)
               min_eigen = eigensolver.eigenvalues()(i,0).real();
           if(eigensolver.eigenvalues()(i,0).real()>max_eigen&&eigensolver.eigenvalues()(i,0).real()<-1e-15||max_eigen==0) //1e-15 is lazy guess at machine epsilon
               max_eigen = eigensolver.eigenvalues()(i,0).real();
       }

       to_Return.Eigenratio = min_eigen/max_eigen;

   }
   else
       to_Return.Eigenratio = 0;

    //Graph Lib Measures
    //Undirected

    //graph types
    //Edge weight necessary for calculations later, default to 1
    typedef boost::property<boost::edge_weight_t, t_weight> Edge_Weight_Property;
    typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::no_property,Edge_Weight_Property> Graph;

    Graph Gen_Graph(nPatches);

    for(int i=0;i<nPatches;i++)
        for(int i2=0;i2<nPatches;i2++)
            if(to_Return.Quotient_Structure[i][i2]>0)
                add_edge(i,i2,Gen_Graph);

    //Setting up vertex & edge weight access

    //typedef boost::property_map<Graph, boost::vertex_index_t>::type Index_Map;
    //Index_Map index = get(boost::vertex_index, Gen_Graph);

    //typedef boost::property_map<Graph,boost::edge_weight_t>::type Weight_Map;
    //Weight_Map Gen_Graph_Weight = get(boost::edge_weight, Gen_Graph);

    //Clustering Coefficient

    typedef boost::exterior_vertex_property<Graph, double> Clustering_Property;
    typedef Clustering_Property::container_type Clustering_Container;
    typedef Clustering_Property::map_type Clustering_Map;

    Clustering_Container coefs(num_vertices(Gen_Graph));
    Clustering_Map Gen_Graph_Clustering(coefs,Gen_Graph);

    to_Return.Clustering_Coefficient = all_clustering_coefficients(Gen_Graph,Gen_Graph_Clustering);

    //Shortest Paths

    typedef boost::exterior_vertex_property<Graph, t_weight> Distance_Property;
    typedef Distance_Property::matrix_type Distance_Matrix;
    typedef Distance_Property::matrix_map_type Distance_Matrix_Map;

    Distance_Matrix Distances_Mat(num_vertices(Gen_Graph));
    Distance_Matrix_Map Gen_Graph_Distances(Distances_Mat,Gen_Graph);

    bool valid = floyd_warshall_all_pairs_shortest_paths(Gen_Graph,Gen_Graph_Distances);

    t_weight Paths_Sum = 0;

    if(valid)
    {
        for(int i=0;i<nPatches-1;i++)
            for(int i2=i+1;i2<nPatches;i2++)
                if(i!=i2)
                    Paths_Sum+=Distances_Mat[i][i2];
    }

    to_Return.Average_Path_Length = (double)(1/(double)nPatches*((double)nPatches-1))*((double)Paths_Sum);


    return to_Return;
}
