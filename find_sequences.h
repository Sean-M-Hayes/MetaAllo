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
#ifndef FIND_SEQUENCES_H
#define FIND_SEQUENCES_H

#include "forward_declarations.h"
#include <vector>

//class Element_Sequence_Set{

//private:
//    std::vector<std::vector<int>> Element_Sequence;
//    std::vector<std::vector<double>> Element_Key;

//public:
//    Element_Sequence_Set(){}

//    Element_Sequence_Set(const std::vector<std::vector<double>> &Element_Criteria, const double &Deviation_Threshold);

//    Element_Sequence_Set(const std::vector<boost_matrix> &Element_Criteria, const std::vector<double> &Deviation_Threshold);

//    void operator()(const Element_Sequence_Set &New_Values);

//    std::vector<double> Get_Element_Key(int index) const;

//    const std::vector<std::vector<double>>& Get_Element_Key() const;

//    const std::vector<std::vector<int>>& Get_Element_Sequence() const;
//};

//class Sequence_Set
//{

//private:
//    //std::vector<std::vector<int>> Sequence;

//    std::vector<std::vector<std::vector<int>>> Sequence;
//    std::vector<boost_matrix> Key;

//    std::vector<double> Base_CV;
//    std::vector<double> S_Dbw_Center_Fit;

//public:
//    Sequence_Set(){}
//    Sequence_Set(const std::vector<boost_matrix> &Element_Dimensions,
//                 const double Target_CV,
//                 const double Threshold_to_Increase_Centers,
//                 const double Convergence_Criteria,
//                 const int Recenter_Iterations,
//                 const int New_Center_Iterations,
//                 const int Min_Centers,
//                 const int Max_Centers);

//    void operator()(const Sequence_Set &New_Values);

//    std::vector<double> Get_Key (int index) const;

//    const boost_matrix& Get_Key() const;

//    //const std::vector<std::vector<int>>& Get_Sequence() const;

//    const std::vector<std::vector<std::vector<int>>>& Get_Sequence() const; //

//    const double& Get_Center_Fit() const;

//    const double& Get_Element_CV() const;

//};

class Sequence
{
private:
    std::vector<std::vector<int>> Sequence_List;
    boost_matrix Key;
    double Base_Variability; //
    double S_Dbw_Center_Fit;

public:
    Sequence(const boost_matrix &Element_Dimensions,
             const std::vector<bool> Scale_SD_by_Mean,
             const std::vector<double> Else_Scale_SD_by,
             const double Target_Var,
             const double Threshold_to_Increase_Centers,
             const double Convergence_Criteria,
             const int Recenter_Iterations,
             const int New_Center_Iterations,
             const int Min_Centers,
             const int Max_Centers);

    void operator()(const Sequence &New_Values);

    std::vector<double> Get_Key(int index) const;

    const boost_matrix& Get_Key() const;

    const std::vector<std::vector<int>>& Get_Sequence() const;

    const double& Get_Center_Fit() const;

    const double& Get_Element_Variability() const; //
};

#endif // FIND_SEQUENCES_H
