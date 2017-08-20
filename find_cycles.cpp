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
#include "find_cycles.h"


//Constructor -> contains all other methods
Find_Cycles::Find_Cycles(const std::pair<boost_matrix_1d_view, boost_matrix_2d_view> &Time_series,
                         double Simulation_Resolution,
                         const Cycle_Detection_Parameters &Algorithm_Parameters,
                         bool Measure_if_Nonstationary):
    Max_Cycle_Details(Time_series.second.shape()[1]),
    Fixed_Point(Time_series.second.shape()[1],0),
    Extrema_to_Compare(Time_series,
                       Simulation_Resolution,
                       Algorithm_Parameters.Minimum_Extrema_Proportion_Difference,
                       Algorithm_Parameters.Difference_Error_Tolerance)
{   
    //Extrema_Stats Extrema_to_Compare(Time_series);
    std::vector<boost_matrix> Extrema_Stats_Matrix = Extrema_to_Compare.Get_Cycle_Stats();
    std::vector<std::vector<double>> Cycle_Times = Extrema_to_Compare.Get_Cycle_Times();
    std::vector<std::vector<std::pair<int,int>>> Raw_Cycle_Indicies = Extrema_to_Compare.Get_Cycle_Indicies();

    for(int i=0;i<Time_series.second.shape()[1];i++)
        if(Extrema_Stats_Matrix[i].shape()[0]==0)
            if(std::abs(Time_series.second[0][i]-Time_series.second[Time_series.second.shape()[0]-1][i])<Algorithm_Parameters.Difference_Error_Tolerance)
            {
                Fixed_Point[i] = true;

                Calculate_Fixed_Point_Stats(i,Time_series,Max_Cycle_Details[i]);
            }

    if(*std::min_element(Fixed_Point.begin(),Fixed_Point.end())>0) //Everything is Fixed & Stationary
    {
        Stationary = true;
    }
    else
    {
        for(int i=0;i<Extrema_Stats_Matrix.size();i++)
        {
            //Scales period variation by 10*Simulation resolution
            //Variation in period will not scale with period length, but simulation resolution
            //Estimates of everything else will scale with size

            std::vector<bool> Scale_SD_by_Mean(1,false);
            std::vector<double> Scale_SD(1,10*Simulation_Resolution);

            Sequences_to_Compare.push_back(Sequence(Extrema_Stats_Matrix[i],
                                                    Scale_SD_by_Mean,
                                                    Scale_SD,
                                                    Algorithm_Parameters.Sequence_Center_Target_CV,
                                                    Algorithm_Parameters.Sequence_Threshold_to_Increase_Centers,
                                                    Algorithm_Parameters.Sequence_Convergence_Criteria,
                                                    Algorithm_Parameters.Sequence_Recenter_Iterations,
                                                    Algorithm_Parameters.Sequence_New_Center_Iterations,
                                                    Algorithm_Parameters.Sequence_Min_Centers,
                                                    Algorithm_Parameters.Sequence_Max_Centers));

            Max_Cycle_Details[i].Center_Fit = Sequences_to_Compare[i].Get_Center_Fit();
            Max_Cycle_Details[i].Extrema_Variability = Sequences_to_Compare[i].Get_Element_Variability();
            //Max_Cycle_Details[i].Extrema_SD = Sequences_to_Compare[i].Get_Element_SD();
        }

//        Sequences_to_Compare(Sequence_Set(Extrema_Stats_Matrix,
//                                          Algorithm_Parameters.Sequence_Center_Target_CV,
//                                          Algorithm_Parameters.Sequence_Threshold_to_Increase_Centers,
//                                          Algorithm_Parameters.Sequence_Convergence_Criteria,
//                                          Algorithm_Parameters.Sequence_Recenter_Iterations,
//                                          Algorithm_Parameters.Sequence_New_Center_Iterations,
//                                          Algorithm_Parameters.Sequence_Min_Centers,
//                                          Algorithm_Parameters.Sequence_Max_Centers));

        //std::vector<double> Cycle_Periods = Sequences_to_Compare.Get_Key(0);

        double Proportion_Stationary_Series = 0;
        //std::vector<Pattern_Sequence> Patterns_to_Compare;

        for(int i=0;i<Sequences_to_Compare.size();i++)
        {
            Patterns_to_Compare.push_back(Pattern_Sequence(Sequences_to_Compare[i].Get_Sequence(),
                                                           Algorithm_Parameters.Min_Occs_for_Common_Pattern,
                                                           Algorithm_Parameters.Min_Prop_for_Common_Pattern,
                                                           Algorithm_Parameters.Max_Cycle_Switch_Threshold,
                                                           Sequences_to_Compare[i].Get_Sequence().size()));

            if(Algorithm_Parameters.Stationary_Cycle_Frequency_Threshold<=Patterns_to_Compare[i].Get_Max_Sequence_Relative_Duration())
                Proportion_Stationary_Series++;
        }

        Proportion_Stationary_Series /= (double)Sequences_to_Compare.size();

        if(Proportion_Stationary_Series>=Algorithm_Parameters.Proportion_Series_for_Stationary)
            Stationary = true;
        else
            Stationary = false;

        if(Stationary||Measure_if_Nonstationary)
        {
            Base_Key = -1;

            std::vector<double> Test_Durations;

            //Find best key of set for comparison
            for(int i=0;i<Patterns_to_Compare.size();i++)
            {
                Test_Durations.push_back(Patterns_to_Compare[i].Get_Max_Sequence_Relative_Duration());
            }

            Base_Key = std::distance(Test_Durations.begin(),std::max_element(Test_Durations.begin(),Test_Durations.end()));

            //End find best key

            //Base_Key will == -1 if there are no patterns

            if(Test_Durations[Base_Key]>0)
            {
                std::vector<int> Max_Key = Patterns_to_Compare[Base_Key].Get_Max_Sequence_Key();
                std::vector<int> Max_Map = Patterns_to_Compare[Base_Key].Get_Max_Sequence_Map();
                double Max_Duration = Patterns_to_Compare[Base_Key].Get_Max_Sequence_Relative_Duration();

                Calculate_Cycle_Descriptive_Stats(Base_Key,
                                                  Time_series,
                                                  Raw_Cycle_Indicies[Base_Key],
                                                  Sequences_to_Compare[Base_Key].Get_Sequence(),
                                                  Sequences_to_Compare[Base_Key].Get_Key(),
                                                  Max_Key,
                                                  Max_Duration,
                                                  Max_Cycle_Details[Base_Key]);

                Calculate_Cycle_Desc_Difference_Stats(Max_Cycle_Details[Base_Key],Max_Cycle_Details[Base_Key]);

                std::vector<double> Max_Cycle_Times;

                for(int j=0;j<Max_Map.size();j++)
                    for(int k=0;k<Max_Key.size();k++)
                        Max_Cycle_Times.push_back(Cycle_Times[Base_Key][Max_Map[j]+k]);

                //Compare all others to Max
                for(int i=0;i<Patterns_to_Compare.size();i++)
                {
                    if(i!=Base_Key)
                    {
                        std::vector<int> Compare_Key = Patterns_to_Compare[i].Get_Max_Sequence_Key();
                        std::vector<int> Compare_Map = Patterns_to_Compare[i].Get_Max_Sequence_Map();
                        double Compare_Duration = Patterns_to_Compare[i].Get_Max_Sequence_Relative_Duration();

                        if(Compare_Key.size()>0) //Is there a key to compare to?
                        {

                            //Change Compare_Key & Map to be the same as Max_Key & Map
                            //if sizes are equal, they are not already the same, and Max_Key is close in frequency to Compare_Key in the Comparison series
//                            if(!Key_is_Equal(Max_Key,Compare_Key))
//                                if(Max_Key.size()==Compare_Key.size())
//                                {
//                                    //std::vector<int> Test_Map = Patterns_to_Compare[i].Get_Matching_Sequence_Map(Max_Key);
//                                    std::pair<double,std::vector<int>> Test_Pattern = Patterns_to_Compare[i].Get_Matching_Sequence(Max_Key);

//                                    std::vector<int> Test_Map = Test_Pattern.second;
//                                    double Test_Duration = Test_Pattern.first*(double)Test_Map.size();

//                                    if(std::abs(Test_Duration-Compare_Duration)<=Algorithm_Parameters.Max_Cycle_Switch_Threshold)
//                                    {
//                                        Compare_Key = Max_Key;
//                                        Compare_Map = Test_Map;
//                                        Compare_Duration = Test_Duration;
//                                    }
//                                }
                            //End change Compare_Key & Map

                            //Find & Compare Cycle_Details
                            Calculate_Cycle_Descriptive_Stats(i,
                                                              Time_series,
                                                              Raw_Cycle_Indicies[i],
                                                              Sequences_to_Compare[i].Get_Sequence(), //CHANGED
                                                              Sequences_to_Compare[i].Get_Key(),
                                                              Compare_Key,
                                                              Compare_Duration,
                                                              Max_Cycle_Details[i]);

                            Calculate_Cycle_Desc_Difference_Stats(Max_Cycle_Details[Base_Key],Max_Cycle_Details[i]);
                            //End find & compare cycle details

                            //Find smallest difference in timing from Max for each Comparison cycle

                            //double Max_Period = Patterns_to_Compare[Base_Key].Get_Max_Sequence_Duration();

                            std::vector<double> Phase_Diffs;

                                //std::vector<std::vector<double>> Cycle_Times = Extrema_to_Compare.Get_Cycle_Times();

                            std::vector<double> Compare_Cycle_Times;

                            for(int j=0;j<Compare_Map.size();j++)
                                for(int k=0;k<Compare_Key.size();k++)
                                    Compare_Cycle_Times.push_back(Cycle_Times[i][Compare_Map[j]+k]);

                            int Max_Cycle_Start = 0;

                            for(int j=0;j<Compare_Cycle_Times.size();j++)
                            {
                                double Min_Difference = std::abs(Max_Cycle_Times[Max_Cycle_Start]-Compare_Cycle_Times[j]);

                                //double Min_Difference = Compare_Cycle_Times[j]-Max_Cycle_Times[Max_Cycle_Start];

                                for(int m=Max_Cycle_Start+1;m<Max_Cycle_Times.size();m++)
                                {
                                    double Current_Difference = std::abs(Max_Cycle_Times[m]-Compare_Cycle_Times[j]);
                                    //double Current_Difference = Compare_Cycle_Times[j]-Max_Cycle_Times[m];

                                    if(Min_Difference-Current_Difference>Algorithm_Parameters.Difference_Error_Tolerance)
                                        Min_Difference = Current_Difference;
                                    else
                                    {
                                        Max_Cycle_Start = m-1;

                                        Phase_Diffs.push_back((Compare_Cycle_Times[j]-Max_Cycle_Times[m-1]));

                                        break;
                                    }
                                }
                            }
                            //End find smallest difference

                            boost_matrix Phase_Diff_Mat(boost::extents[Phase_Diffs.size()][1]);

                            for(int j=0;j<Phase_Diff_Mat.shape()[0];j++)
                                Phase_Diff_Mat[j][0] = Phase_Diffs[j];

                            //Find patterns and stats of phase differences
                            Calculate_Cycle_Phase_Difference_Stats(Phase_Diff_Mat,
                                                                   Algorithm_Parameters,
                                                                   Simulation_Resolution,
                                                                   Max_Cycle_Details[i]);
                        }
                    }
                }
            }
        }
    }
}

//Methods called by constructor

bool Find_Cycles::Key_is_Equal(std::vector<int> Compare_1, std::vector<int> Compare_2)
{
    if(Compare_1.size()==Compare_2.size())
    {
        bool all_elements_matched = true;
        for(int i=0;i<Compare_1.size();i++)
        {
            if(Compare_1[i]!=Compare_2[i])
            {
                all_elements_matched=false;
                break;
            }
        }

        if(all_elements_matched)
            return true;
    }
    return false;
}

void Find_Cycles::Calculate_Fixed_Point_Stats(const int column,
                                              const std::pair<boost_matrix_1d_view, boost_matrix_2d_view> &Time_series,
                                              Cycle_Details &Output)
{
    Output.Extrema_Variability = 0;
    Output.Center_Fit = 0;

    Output.Relative_Duration = 1;

    Output.Maxima_Timing.clear();
    Output.Maxima.clear();
    Output.Minima.clear();
    Output.Period = 0;

    Output.Maxima_Timing.push_back(0);
    Output.Maxima.push_back(Time_series.second[0][column]);
    Output.Minima.push_back(Time_series.second[0][column]);

    Output.Max_Amplitude = 0;

    Output.Example_Time_series.clear();

    Output.Example_Time_series.push_back(std::pair<double,double>(0,Time_series.second[0][column]));

    Output.Average = Time_series.second[0][column];
    Output.SD = 0;
    Output.Skewness = 0;

    Output.Max_Maxima_Diff = 0;
    Output.Min_Minima_Diff = 0;
    Output.Period_Diff = 0;
    Output.Max_Amplitude_Diff = 0;

    Output.Average_Diff = 0;
    Output.SD_Diff = 0;
    Output.Skewness_Diff = 0;

    Output.Common_Phase_Diff.clear();
    Output.Common_Phase_Diff.push_back(0);

    Output.Common_Phase_Diff_Average = 0;
    Output.Common_Phase_Diff_SD = 0;

    Output.Common_Phase_Diff_Frequency = 0;

    Output.Average_Phase_Diff = 0;
    Output.SD_Phase_Diff = 0;
}

void Find_Cycles::Calculate_Cycle_Descriptive_Stats(const int column,
                                                    const std::pair<boost_matrix_1d_view, boost_matrix_2d_view> &Time_series,
                                                    const std::vector<std::pair<int,int>> &Raw_Cycle_Indicies,
                                                    const std::vector<std::vector<int>> &Element_Sequence,
                                                    const boost_matrix &Element_Key,
                                                    const std::vector<int> &Cycle_Sequence,
                                                    const double Cycle_Duration,
                                                    Cycle_Details &Output)
{

    //Assumed Element_Key is Duration, Maxima 1, Maxima 2, Minima

    Output.Relative_Duration = Cycle_Duration;

    Output.Maxima_Timing.clear();
    Output.Maxima.clear();
    Output.Minima.clear();
    Output.Period = 0;

    for(int i=0;i<Cycle_Sequence.size();i++)
    {
        Output.Maxima_Timing.push_back(Element_Key[Cycle_Sequence[i]][0]);
        Output.Maxima.push_back(Element_Key[Cycle_Sequence[i]][1]);
        Output.Minima.push_back(Element_Key[Cycle_Sequence[i]][3]);

        Output.Period+=Element_Key[Cycle_Sequence[i]][0];
    }

    Output.Max_Amplitude = *std::max_element(Output.Maxima.begin(),Output.Maxima.end())-*std::min_element(Output.Minima.begin(),Output.Minima.end());

    Output.Example_Time_series.clear();
    for(int i=0;i<Element_Sequence.size();i++)
    {
        bool match = true;
        for(int j=0;j<Cycle_Sequence.size();j++)
        {
            bool fuzzy_match = false;
            for(int k=0;k<Element_Sequence[i+j].size();k++)
            {
                if(Element_Sequence[i+j][k]==Cycle_Sequence[j])
                {
                    fuzzy_match = true;
                    break;
                }
            }

            if(!fuzzy_match)
            {
                match=false;
                break;
            }
        }
        if(match)
        {
            for(int j=0;j<Cycle_Sequence.size();j++)
                for(int k=Raw_Cycle_Indicies[i+j].first;k<Raw_Cycle_Indicies[i+j].second;k++)
                    Output.Example_Time_series.push_back(std::pair<double,double>(Time_series.first[k],Time_series.second[k][column]));

            break;
        }
    }

    if(Output.Example_Time_series.size()>0)
    {

    Output.Average = 0;

    for(int i=0;i<Output.Example_Time_series.size();i++)
        Output.Average += Output.Example_Time_series[i].second;

    Output.Average /= (double)Output.Example_Time_series.size();

    Output.SD = 0;
    Output.Skewness = 0;

    for(int i=0;i<Output.Example_Time_series.size();i++)
    {
        Output.SD += pow(Output.Example_Time_series[i].second-Output.Average,2);
        Output.Skewness += pow(Output.Example_Time_series[i].second-Output.Average,3);
    }

    Output.SD /= (double)Output.Example_Time_series.size();
    Output.SD = pow(Output.SD,.5);

    Output.Skewness /= (double)Output.Example_Time_series.size();
    Output.Skewness /= pow(Output.SD,3);
    }
    else
    {
        Output.Average = 0; //Something fucky is happening
    }
}

void Find_Cycles::Calculate_Cycle_Desc_Difference_Stats(const Cycle_Details &Base, Cycle_Details &Output)
{
    Output.Max_Maxima_Diff = *std::max_element(Output.Maxima.begin(),Output.Maxima.end()) - *std::max_element(Base.Maxima.begin(),Base.Maxima.end());
    Output.Min_Minima_Diff = *std::min_element(Output.Minima.begin(),Output.Minima.end()) - *std::min_element(Base.Minima.begin(),Base.Minima.end());
    Output.Period_Diff = Output.Period - Base.Period;
    Output.Max_Amplitude_Diff = Output.Max_Amplitude - Base.Max_Amplitude;

    Output.Average_Diff = Output.Average - Base.Average;
    Output.SD_Diff = Output.SD-Base.SD;
    Output.Skewness_Diff = Output.Skewness-Base.Skewness;
}

void Find_Cycles::Calculate_Cycle_Phase_Difference_Stats(boost_matrix &Phase_Diffs,
                                                         const Cycle_Detection_Parameters &Algorithm_Parameters,
                                                         const double Simulation_Resolution,
                                                         Cycle_Details &Output)
{
    Output.Common_Phase_Diff.clear();
    Output.Common_Phase_Diff_Average = 0;
    Output.Common_Phase_Diff_SD = 0;
    Output.Common_Phase_Diff_Frequency = 0;

    if(Phase_Diffs.shape()[0]>0)
    {
        std::vector<bool> Scale_SD_by_Mean(1,false);
        std::vector<double> Else_Scale_SD_by(1,10*Simulation_Resolution);

        Sequence Process_Phases(Phase_Diffs,
                                    Scale_SD_by_Mean,
                                    Else_Scale_SD_by,
                                    Algorithm_Parameters.Sequence_Center_Target_CV,
                                    Algorithm_Parameters.Sequence_Threshold_to_Increase_Centers,
                                    Algorithm_Parameters.Sequence_Convergence_Criteria,
                                    Algorithm_Parameters.Sequence_Recenter_Iterations,
                                    Algorithm_Parameters.Sequence_New_Center_Iterations,
                                    Algorithm_Parameters.Sequence_Min_Centers,
                                    Algorithm_Parameters.Sequence_Max_Centers);

        Pattern_Sequence Phase_Pattern(Process_Phases.Get_Sequence(),
                Algorithm_Parameters.Min_Occs_for_Common_Pattern,
                Algorithm_Parameters.Min_Prop_for_Common_Pattern,
                Algorithm_Parameters.Max_Cycle_Switch_Threshold,
                Phase_Diffs.shape()[0]);

        std::vector<int> Phase_Pattern_Map = Phase_Pattern.Get_Max_Sequence_Map();

        if(Phase_Pattern_Map.size()>0)
        {
            std::vector<int> Phase_Pattern_Key = Phase_Pattern.Get_Max_Sequence_Key();

            for(int i=0;i<Phase_Pattern_Key.size();i++)
            {
                Output.Common_Phase_Diff.push_back(Process_Phases.Get_Key()[Phase_Pattern_Key[i]][0]);
                Output.Common_Phase_Diff_Average+=Process_Phases.Get_Key()[Phase_Pattern_Key[i]][0];
            }

            Output.Common_Phase_Diff_Average/=(double)Phase_Pattern_Key.size();

            Output.Common_Phase_Diff_SD = 0;

            for(int i=0;i<Phase_Pattern_Key.size();i++)
                Output.Common_Phase_Diff_SD+=pow(Process_Phases.Get_Key()[Phase_Pattern_Key[i]][0]-Output.Common_Phase_Diff_Average,2);

            Output.Common_Phase_Diff_SD/=(double)Phase_Pattern_Key.size();
            Output.Common_Phase_Diff_SD=pow(Output.Common_Phase_Diff_SD,.5);

            Output.Common_Phase_Diff_Frequency = (Phase_Pattern.Get_Max_Sequence_Relative_Duration());
        }

        Output.Average_Phase_Diff = 0;

        for(int j=0;j<Phase_Diffs.shape()[0];j++)
        {
            Output.Average_Phase_Diff+=Phase_Diffs[j][0];
        }

        Output.Average_Phase_Diff/=(double)Phase_Diffs.shape()[0];

        Output.SD_Phase_Diff = 0;

        for(int j=0;j<Phase_Diffs.shape()[0];j++)
        {
            Output.SD_Phase_Diff+=pow(Phase_Diffs[j][0]-Output.Average_Phase_Diff,2);
        }

        Output.SD_Phase_Diff/=(double)Phase_Diffs.shape()[0];
        Output.SD_Phase_Diff=pow(Output.SD_Phase_Diff,.5);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////

//Convenience function

std::pair<boost_matrix_1d_view,boost_matrix_2d_view> Subset_Time_Series(boost_matrix &Matrix, std::pair<int, int> Test_Interval, int Species, int nSpecies)
{
    typedef boost::multi_array_types::index_range range_t;
    boost_matrix::index_gen indices;

    if(Test_Interval.first<Matrix.shape()[0]&&Test_Interval.second<Matrix.shape()[0])
    {

        boost_matrix_1d_view Times = Matrix[indices[range_t().start(Test_Interval.first).finish(Test_Interval.second)][0]];
        boost_matrix_2d_view Species_Data = Matrix[indices[range_t().start(Test_Interval.first).finish(Test_Interval.second)][range_t().start(Species+1).stride(nSpecies)]];

        return std::pair<boost_matrix_1d_view,boost_matrix_2d_view>(Times,Species_Data);
    }
    else
    {
        boost_matrix_1d_view Times = Matrix[indices[range_t().start(0).finish(Matrix.shape()[0]-1)][0]];
        boost_matrix_2d_view Species_Data = Matrix[indices[range_t().start(0).finish(Matrix.shape()[0]-1)][range_t().start(Species+1).stride(nSpecies)]];

        return std::pair<boost_matrix_1d_view,boost_matrix_2d_view>(Times,Species_Data);
    }
}
