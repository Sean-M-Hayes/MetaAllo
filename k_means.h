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
#ifndef K_MEANS_H
#define K_MEANS_H
#include "forward_declarations.h"
#include <vector>

struct Normalized_Data
{
    std::vector<double> Means;
    std::vector<double> SDs;

    std::vector<double> Normalized_Means;
    std::vector<double> Normalized_SDs;

    boost_matrix Data;

    Normalized_Data(const std::vector<boost_matrix> &In_Data);

    Normalized_Data(const boost_matrix &In_Data);
};

struct Node
{
    //In violation of rule of three, be careful with those pointers
    std::vector<int> Point_Map;
    std::vector<double> Centroid;
    std::vector<double> Minima;
    std::vector<double> Maxima;

    Node* left_child;
    Node* right_child;

    Node(std::vector<int> Points, const boost_matrix &Point_List, int grow_limit, double target_cv); // target_cv is new
    void Update_Centers(std::vector<int> Candidates, boost_matrix &Centers, std::vector<std::vector<int> > &Centers_Map);

    ~Node();
};

struct k_means_tree_filtering
{
    bool Centers_Failed;
    bool Convergence_Failed;
    double S_Dbw_Center_Fit;
    boost_matrix Centers;
    std::vector<std::vector<int>> Centers_Map;

    k_means_tree_filtering(const boost_matrix &Test_Vectors,
                           const std::vector<double> &Test_SDs,
                           int k_Centers,
                           double Convergence_Criteria,
                           int Niter,
                           Node &kd_tree_root);

    void operator=(k_means_tree_filtering New_Values);
};

double S_Dbw_Cluster_Goodness_of_Fit(const boost_matrix &Centers, const std::vector<std::vector<int>> &Centers_Map, const boost_matrix &Points_List, const std::vector<double> &Point_SDs);

std::vector<int> Build_Element_Sequence(int Sequence_Length, const k_means_tree_filtering &Clustering_Fit);
std::vector<std::vector<int>> Build_Fuzzy_Element_Sequence(const boost_matrix &Data, const k_means_tree_filtering &Clustering_Fit);

void Write_Clustering_Data(std::string Write_Path, const boost_matrix &Extrema, const std::vector<int> &Cluster_Assignment);
void Write_Test_Clustering_Data(std::string Write_Path, const boost_matrix &Extrema, const std::vector<std::vector<int>> &Cluster_Assignment);
void Write_Cluster_Distance(std::string Write_Path, const boost_matrix &Extrema, const std::vector<int> &Cluster_Assignment, const k_means_tree_filtering &Cluster_Fit);

#endif // K_MEANS_H
