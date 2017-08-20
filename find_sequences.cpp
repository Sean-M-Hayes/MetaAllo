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
#include "find_sequences.h"
#include "k_means.h"

Sequence::Sequence(const boost_matrix &Element_Dimensions,
                           const std::vector<bool> Scale_SD_by_Mean,
                           const std::vector<double> Else_Scale_SD_by,
                           const double Target_Val,
                           const double Threshold_to_Increase_Centers,
                           const double Convergence_Criteria,
                           const int Recenter_Iterations,
                           const int New_Center_Iterations,
                           const int Min_Centers,
                           const int Max_Centers)
{    
    Normalized_Data Normalized_Dimensions(Element_Dimensions);

    Base_Variability = 0;

    for(int i=0;i<Normalized_Dimensions.Data.shape()[1];i++)
    {
        if(Scale_SD_by_Mean.size()>i)
            if(!Scale_SD_by_Mean[i]&&Else_Scale_SD_by.size()>i)
            {
                Base_Variability += pow(Normalized_Dimensions.SDs[i]/Else_Scale_SD_by[i],2);
                continue;
            }

        Base_Variability += pow(Normalized_Dimensions.SDs[i]/Normalized_Dimensions.Means[i],2);
    }

    Base_Variability = sqrt(Base_Variability);

    std::vector<int> Points;

    for(int i=0;i<Normalized_Dimensions.Data.shape()[0];i++)
        Points.push_back(i);

    //Fix this
    Node kd_tree_root(Points,Normalized_Dimensions.Data,1,Target_Val/Base_Variability); //Tree will terminate and force grouping of points when (CV nodes)*Base_CV = Target_Var

    k_means_tree_filtering Best_Fit(Normalized_Dimensions.Data,Normalized_Dimensions.Normalized_SDs,1,Convergence_Criteria,Recenter_Iterations,kd_tree_root);

    if(Base_Variability>Target_Val)
    {
        for(int i=Min_Centers;i<=Max_Centers;i++)
        {
            for(int j=0;j<New_Center_Iterations;j++)
            {
                k_means_tree_filtering Test(Normalized_Dimensions.Data,Normalized_Dimensions.Normalized_SDs,i,Convergence_Criteria,Recenter_Iterations,kd_tree_root);

                if(Test.Centers_Failed)
                    goto Escape;

                if((Best_Fit.S_Dbw_Center_Fit-Test.S_Dbw_Center_Fit)>Threshold_to_Increase_Centers)
                    Best_Fit = Test;
                else if((Best_Fit.S_Dbw_Center_Fit-Test.S_Dbw_Center_Fit)>0&&Test.Centers.shape()[0]<Best_Fit.Centers.shape()[0])
                    Best_Fit = Test;
            }

            if(Best_Fit.S_Dbw_Center_Fit*Base_Variability<=Target_Val)
                break;
        }
    }

    Escape:

    Sequence_List = Build_Fuzzy_Element_Sequence(Normalized_Dimensions.Data,Best_Fit);

    Key.resize(boost::extents[Best_Fit.Centers.shape()[0]][Best_Fit.Centers.shape()[1]]);

    for(int i=0;i<Best_Fit.Centers.shape()[0];i++)
        for(int j=0;j<Best_Fit.Centers.shape()[1];j++)
            Key[i][j] = Best_Fit.Centers[i][j]*Normalized_Dimensions.SDs[j]+Normalized_Dimensions.Means[j];

    S_Dbw_Center_Fit = Best_Fit.S_Dbw_Center_Fit;

}

void Sequence::operator ()(const Sequence &New_Values)
{
    Sequence_List = New_Values.Get_Sequence();

    boost_matrix New_Key = New_Values.Get_Key();

    Key.resize(boost::extents[New_Key.shape()[0]][New_Key.shape()[1]]);

    for(int i=0;i<Key.shape()[0];i++)
        for(int j=0;j<Key.shape()[1];j++)
            Key[i][j] = New_Key[i][j];

    Base_Variability = New_Values.Get_Element_Variability();
    S_Dbw_Center_Fit = New_Values.Get_Center_Fit();
}

std::vector<double> Sequence::Get_Key(int index) const
{
    std::vector<double> to_Return;

    for(int i=0;i<Key.shape()[0];i++)
        to_Return.push_back(Key[i][index]);

    return to_Return;
}

const boost_matrix& Sequence::Get_Key() const
{
    return Key;
}

const std::vector<std::vector<int>>& Sequence::Get_Sequence() const
{
    return Sequence_List;
}

const double& Sequence::Get_Center_Fit() const
{
    return S_Dbw_Center_Fit;
}

const double& Sequence::Get_Element_Variability() const
{
    return Base_Variability;
}



