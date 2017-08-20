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
#ifndef SPATIAL_STRUCTURE_METHODS_H
#define SPATIAL_STRUCTURE_METHODS_H

#include "forward_declarations.h"

struct Spatial_Structure_Metadata{

    boost_matrix Spatial_Structure;

    double Average_Degree;
    double Degree_Skewness;

    double Average_Path_Length;
    double Clustering_Coefficient;
    double Eigenratio;

    double Connectance;

    std::vector<double> Degree;

    Spatial_Structure_Metadata();

    Spatial_Structure_Metadata(int nPatches);

    void operator()( Spatial_Structure_Metadata New_Values );

    bool is_Equal(const Spatial_Structure_Metadata &Compare_To);

};

struct Quotient_Structure_Metadata{

    boost_matrix Quotient_Structure;

    double Average_In_Degree;
    double Average_Out_Degree;

    double In_Degree_Skewness;
    double Out_Degree_Skewness;

    std::vector<double> In_Degree;
    std::vector<double> Out_Degree;

    double Average_Path_Length;
    double Clustering_Coefficient;
    double Eigenratio;

    Quotient_Structure_Metadata();
};

Spatial_Structure_Metadata Erdos_Renyi_Undirected(int nPatches, double Connectance);

Spatial_Structure_Metadata Dendritic_Undirected(int nPatches, double Branching_Probability);

Spatial_Structure_Metadata Radial_Tree_Undirected(int nPatches, int nBranches);

Spatial_Structure_Metadata Ring_Lattice_Undirected(int nPatches, int Neighbor_Distance);

Spatial_Structure_Metadata Calculate_Unweighted_Symmetrical_Spatial_Structure_Metadata(boost_matrix Spatial_Structure);

Spatial_Structure_Metadata Permuted_Shipley_Skinner(std::vector<int> Remove_Which);

Quotient_Structure_Metadata Calculate_Quotient_Structure_Metadata(const boost_matrix &Spatial_Structure, std::vector<int> Synchronized_Clusters);

bool is_Isomorphic(boost_matrix Structure_A, boost_matrix Structure_B);

#endif // SPATIAL_STRUCTURE_METHODS_H
