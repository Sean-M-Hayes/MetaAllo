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
#include "find_extrema.h"

std::pair<double,double> Quadratic_Interpolation(const double first, const double middle, const double last)
{
    double b = (last-first)/2;
    double c = (last-2*middle+first)/2;

    double time = -b/2/c;

    return std::pair<double,double>(time,middle+b*time+c*pow(time,2));
}

Extrema_Stats::Extrema_Stats(const std::pair<boost_matrix_1d_view, boost_matrix_2d_view> &Time_series,double Simulation_Resolution, double Minimum_Proportion_Difference, double Minimum_Absolute_Difference)
{


    for(int i=0;i<Time_series.second.shape()[1];i++)
    {
        std::vector<int> first_pre_max;

        std::vector<double> minima_time;
        std::vector<double> minima_val;

        std::vector<double> maxima_time;
        std::vector<double> maxima_val;

        for(int j=1;j<Time_series.second.shape()[0]-1;j++)
            if(Time_series.second[j-1][i]>=Time_series.second[j][i]&&Time_series.second[j+1][i]>Time_series.second[j][i]) //Is it a minima?
            {
                std::pair<double,double> test_x(Quadratic_Interpolation(Time_series.second[j-1][i],Time_series.second[j][i],Time_series.second[j+1][i]));

                test_x.first *= Simulation_Resolution;

                if(minima_time.size()==0) //no minima added yet
                {
                    if(maxima_time.size()>0) //no minima, but maxima added
                    {
                        if((maxima_val.back()-test_x.second)/((maxima_val.back()+test_x.second)/2)>Minimum_Proportion_Difference&&maxima_val.back()-test_x.second>Minimum_Absolute_Difference)
                        {
                            minima_time.push_back(Time_series.first[j]+test_x.first);
                            minima_val.push_back(test_x.second);
                        }
                    }
                    else //no minima or maxima added
                    {
                        minima_time.push_back(Time_series.first[j]+test_x.first);
                        minima_val.push_back(test_x.second);
                    }
                }
                else if(maxima_time.size()==0) //there is a preceding minima, but no maxima
                {
                    if(minima_val.back()>test_x.second)
                    {
                        minima_time.back() = Time_series.first[j]+test_x.first;
                        minima_val.back() = test_x.second;
                    }
                }
                else if(minima_time.back()< maxima_time.back()) //a maxima was found after the last minima
                {
                    double Test = (maxima_val.back()-test_x.second)/((maxima_val.back()+test_x.second)/2);

                    if((maxima_val.back()-test_x.second)/((maxima_val.back()+test_x.second)/2)>Minimum_Proportion_Difference&&maxima_val.back()-test_x.second>Minimum_Absolute_Difference)
                    {
                        minima_time.push_back(Time_series.first[j]+test_x.first);
                        minima_val.push_back(test_x.second);
                    }
                }
                else if( minima_val.back()>test_x.second) //the last extrema found was a minima
                {
                    minima_time.back() = Time_series.first[j]+test_x.first;
                    minima_val.back() = test_x.second;
                }
            }
            else if(Time_series.second[j-1][i]<=Time_series.second[j][i]&&Time_series.second[j+1][i]<Time_series.second[j][i]) //Is it a maxima?
            {
                std::pair<double,double> test_x(Quadratic_Interpolation(Time_series.second[j-1][i],Time_series.second[j][i],Time_series.second[j+1][i]));
                test_x.first *= Simulation_Resolution;

                if(test_x.second>Minimum_Absolute_Difference)
                {

                if(maxima_time.size()==0) //no maxima added yet
                {
                    if(minima_time.size()>0) //no maxima, but minima added
                    {
                        if((test_x.second-minima_val.back())/((minima_val.back()+test_x.second)/2)>Minimum_Proportion_Difference&&test_x.second-minima_val.back()>Minimum_Absolute_Difference)
                        {
                            maxima_time.push_back(Time_series.first[j]+test_x.first);
                            maxima_val.push_back(test_x.second);

                            if(test_x.first>0)
                                first_pre_max.push_back(j);
                            else
                                first_pre_max.push_back(j-1);
                        }
                    }
                    else //no maxima or minima added
                    {
                    maxima_time.push_back(Time_series.first[j]+test_x.first);
                    maxima_val.push_back(test_x.second);

                    if(test_x.first>0)
                        first_pre_max.push_back(j);
                    else
                        first_pre_max.push_back(j-1);
                    }
                }
                else if(minima_time.size()==0) //there is a preceding maxima, but no minima
                {
                    if(maxima_val.back()<test_x.second)
                    {
                        maxima_time.back()=Time_series.first[j]+test_x.first;
                        maxima_val.back()=test_x.second;

                        if(test_x.first>0)
                            first_pre_max.back()=j;
                        else
                            first_pre_max.back()=j-1;
                    }
                }
                else if(maxima_time.back()<minima_time.back()) //a minima was found after the last maxima
                {
                    if((test_x.second-minima_val.back())/((minima_val.back()+test_x.second)/2)>Minimum_Proportion_Difference&&test_x.second-minima_val.back()>Minimum_Absolute_Difference)
                    {

                        maxima_time.push_back(Time_series.first[j]+test_x.first);
                        maxima_val.push_back(test_x.second);

                        if(test_x.first>0)
                            first_pre_max.push_back(j);
                        else
                            first_pre_max.push_back(j-1);
                    }
                }
                else if(maxima_val.back()<test_x.second) //the last extrema found was a maxima
                {
                    maxima_time.back()=Time_series.first[j]+test_x.first;
                    maxima_val.back()=test_x.second;

                    if(test_x.first>0)
                        first_pre_max.back()=j;
                    else
                        first_pre_max.back()=j-1;
                }
                }
            }

        if(maxima_time.size()>3&&minima_time.size()>2)
        {
            //Set to skip first and last, which may be inaccurate in some cases
            boost_matrix Maxima_Cycle_Stats(boost::extents[maxima_time.size()-3][4]);

            std::vector<std::pair<int,int>> Maxima_Cycle_Indicies;

            std::vector<double> Maxima_Cycle_Times;

            for(int j=2;j<maxima_time.size()-1;j++)
            {
                Maxima_Cycle_Stats[j-2][0] = maxima_time[j]-maxima_time[j-1];
                Maxima_Cycle_Stats[j-2][1] = maxima_val[j-1];
                Maxima_Cycle_Stats[j-2][2] = maxima_val[j];

                if(minima_time[j-1]>maxima_time[j-1]&&minima_time[j-1]<maxima_time[j])
                    Maxima_Cycle_Stats[j-2][3] = minima_val[j-1];
                else if (minima_time.size()>j)
                {
                        if(minima_time[j]<maxima_time[j]&&minima_time[j]>maxima_time[j-1])
                            Maxima_Cycle_Stats[j-2][3] = minima_val[j];
                        else
                        {
                            //Shouldn't happen
                            Maxima_Cycle_Stats.resize(boost::extents[j-3][4]);
                            break;
                        }
                }
                else
                {
                    //Shouldn't happen
                    Maxima_Cycle_Stats.resize(boost::extents[j-3][4]);
                    break;
                }


                Maxima_Cycle_Indicies.push_back(std::pair<int,int>(first_pre_max[j-1],first_pre_max[j]));

                Maxima_Cycle_Times.push_back(maxima_time[j-1]);
            }

            Maxima_Cycles.push_back(Cycles(Maxima_Cycle_Indicies,Maxima_Cycle_Times,Maxima_Cycle_Stats));
        }
        else
        {
            Maxima_Cycles.push_back(Cycles());
        }
    }
}

//Used

std::vector<boost_matrix> Extrema_Stats::Get_Cycle_Stats() const
{
    std::vector<boost_matrix> to_Return;

    for(int i=0;i<Maxima_Cycles.size();i++)
        to_Return.push_back(Maxima_Cycles[i].Stats);

    return to_Return;
}

std::vector<std::vector<double>> Extrema_Stats::Get_Cycle_Times() const
{
    std::vector<std::vector<double>> to_Return;

    for(int i=0;i<Maxima_Cycles.size();i++)
        to_Return.push_back(Maxima_Cycles[i].Times);

    return to_Return;
}

std::vector<std::vector<std::pair<int,int>>> Extrema_Stats::Get_Cycle_Indicies() const
{
    std::vector<std::vector<std::pair<int,int>>> to_Return;

    for(int i=0;i<Maxima_Cycles.size();i++)
        to_Return.push_back(Maxima_Cycles[i].Indicies);

    return to_Return;
}
