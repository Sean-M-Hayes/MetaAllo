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
#ifndef FIND_EXTREMA_H
#define FIND_EXTREMA_H

#include "forward_declarations.h"
#include <vector>

std::pair<double,double> Quadratic_Interpolation(const double first, const double middle, const double last);

class Extrema_Stats{

private:

    struct Cycles
    {

        std::vector<std::pair<int,int>> Indicies;

        std::vector<double> Times;

        boost_matrix Stats; //Dist, First, Second, Middle (Opposite type of First/Second)

        //Temp for Debugging
        std::vector<double> Dist;
        std::vector<double> First;
        std::vector<double> Second;
        std::vector<double> Middle;

        Cycles(){}

        Cycles(std::vector<std::pair<int,int>> Input_Indicies, std::vector<double> Input_Times, boost_matrix Input_Stats):
            Indicies(Input_Indicies), Times(Input_Times), Stats(Input_Stats)
        {
            for(int i=0;i<Stats.shape()[0];i++)
            {
                Dist.push_back(Stats[i][0]);
                First.push_back(Stats[i][1]);
                Second.push_back(Stats[i][2]);
                Middle.push_back(Stats[i][3]);
            }
        }

    };

    std::vector<Cycles> Maxima_Cycles;

public:

    Extrema_Stats(){}

    Extrema_Stats(const std::pair<boost_matrix_1d_view,boost_matrix_2d_view> &Time_series,
                  double Simulation_Resolution,
                  double Minimum_Proportion_Difference,
                  double Minimum_Absolute_Difference);

    //Used

    std::vector<boost_matrix> Get_Cycle_Stats() const;

    std::vector<std::vector<double>> Get_Cycle_Times() const;

    std::vector<std::vector<std::pair<int,int>>> Get_Cycle_Indicies() const;
};

#endif // FIND_EXTREMA_H
