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
#ifndef FIND_CYCLES_H
#define FIND_CYCLES_H

#include "forward_declarations.h"
#include "find_extrema.h"
#include "find_patterns.h"
#include "find_sequences.h"
#include <vector>

struct Cycle_Detection_Parameters
{
    //Threshold for detecting stationarity
    double Stationary_Cycle_Frequency_Threshold; //How common a cycle must be to be considered stationary
    double Proportion_Series_for_Stationary; //What proportion of series must be stationary

    double Difference_Error_Tolerance;

    double Minimum_Extrema_Proportion_Difference;

    double Sequence_Center_Target_CV;
    double Sequence_Threshold_to_Increase_Centers;
    double Sequence_Convergence_Criteria;
    int Sequence_Recenter_Iterations;
    int Sequence_New_Center_Iterations;
    int Sequence_Min_Centers;
    int Sequence_Max_Centers;

    //Biggest loss in frequency to tolerate when deciding to switch to a less frequent pattern for comparison
    double Max_Cycle_Switch_Threshold;

    int Min_Occs_for_Common_Pattern;
    double Min_Prop_for_Common_Pattern;

    double Round_Times_To;
};

struct Cycle_Details
{
    //Base Descriptors

    double Extrema_Variability;
    double Center_Fit;

    double Relative_Duration;

    std::vector<double> Maxima_Timing;
    std::vector<double> Maxima;
    std::vector<double> Minima;
    double Period;

    double Max_Amplitude;

    std::vector<std::pair<double,double>> Example_Time_series;

    double Average;
    double SD;
    double Skewness;

    //Frequency & Amplitude Comparisons

    double Max_Maxima_Diff;
    double Min_Minima_Diff;
    double Period_Diff;
    double Max_Amplitude_Diff;

    double Average_Diff;
    double SD_Diff;
    double Skewness_Diff;

    //Phase Comparisons

    std::vector<double> Common_Phase_Diff;
    double Common_Phase_Diff_Average;
    double Common_Phase_Diff_SD;

    double Common_Phase_Diff_Frequency;

    double Average_Phase_Diff;
    double SD_Phase_Diff;
};

class Find_Cycles{

public:

    int Base_Key;

//    double Extrema_CV;
//    double Center_Fit;

    bool Stationary; //Rest of cycle descriptors will not be populated if false

    std::vector<int> Fixed_Point; //for Fixed_Point[i]==1, Max_Cycle_Details[i] will be empty

    std::vector<Cycle_Details> Max_Cycle_Details;

    Extrema_Stats Extrema_to_Compare;
    std::vector<Sequence> Sequences_to_Compare;
    std::vector<Pattern_Sequence> Patterns_to_Compare;

    bool Key_is_Equal(std::vector<int> Compare_1, std::vector<int> Compare_2);

    void Calculate_Fixed_Point_Stats(const int column,
                                                  const std::pair<boost_matrix_1d_view, boost_matrix_2d_view> &Time_series,
                                                  Cycle_Details &Output);

    void Calculate_Cycle_Descriptive_Stats(const int column,
                                           const std::pair<boost_matrix_1d_view,boost_matrix_2d_view> &Time_series,
                                           const std::vector<std::pair<int,int>> &Raw_Cycle_Indicies,
                                           const std::vector<std::vector<int> > &Element_Sequence,
                                           const boost_matrix &Element_Key,
                                           const std::vector<int> &Cycle_Sequence,
                                           const double Cycle_Duration,
                                           Cycle_Details &Output);

    void Calculate_Cycle_Desc_Difference_Stats(const Cycle_Details &Base, Cycle_Details &Output);

    void Calculate_Cycle_Phase_Difference_Stats(boost_matrix &Phase_Diffs,
                                                const Cycle_Detection_Parameters &Algorithm_Parameters,
                                                const double Simulation_Resolution,
                                                Cycle_Details &Output);

    Find_Cycles(const std::pair<boost_matrix_1d_view,boost_matrix_2d_view> &Time_series,
                double Simulation_Resolution,
                const Cycle_Detection_Parameters &Algorithm_Parameters,
                bool Measure_if_Nonstationary);

};



std::pair<boost_matrix_1d_view,boost_matrix_2d_view> Subset_Time_Series(boost_matrix &Matrix, std::pair<int,int> Test_Interval, int Species, int nSpecies);


#endif // FIND_CYCLES_H
