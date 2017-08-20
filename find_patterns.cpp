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
#include "find_patterns.h"

Pattern_Sequence::Pattern_Sequence(const std::vector<std::vector<int>> &Element_Sequence,
                                   const int Min_Occurences,
                                   const double Min_Duration,
                                   const double Max_Switch_Threshold,
                                   const int Max_Pattern_Length)
{
    
    double Common_Duration = .5;
    
    //Method makes assumptions about composition of Element Key: must be Dist, first max, second max, min

    //std::vector<double> Sequence_Weighted_Length;
    std::vector<std::vector<int>> Sequence_Pattern_Key;
    std::vector<std::vector<int>> Sequence_Map;
    std::vector<std::vector<int>> Sequence_Size_Map;

    std::vector<std::vector<int>> Sequence_Element_Composition;

    //Total_Weighted_Duration = 0;

    Sequence_Size_Map.push_back(std::vector<int>());
    Common_Sequences_Size_Map.push_back(std::vector<int>());

    std::vector<int> Observed_Elements;

    for(int i=0;i<Element_Sequence.size();i++)
        for(int j=0;j<Element_Sequence[i].size();j++)
        {
            //Total_Weighted_Duration += Element_Weight[Element_Sequence[i]];

            if(Element_Sequence[i][j]>=Sequence_Pattern_Key.size())
            {
                //Sequence_Weighted_Length.resize(Element_Sequence[i]+1);
                Sequence_Pattern_Key.resize(Element_Sequence[i][j]+1);
                Sequence_Map.resize(Element_Sequence[i][j]+1);
            }

            if(Sequence_Pattern_Key[Element_Sequence[i][j]].size()==0)
            {
                Observed_Elements.push_back(Element_Sequence[i][j]);

                //Sequence_Weighted_Length[Element_Sequence[i]]=Element_Weight[Element_Sequence[i]];
                Sequence_Pattern_Key[Element_Sequence[i][j]].push_back(Element_Sequence[i][j]);
                Sequence_Size_Map.back().push_back(Element_Sequence[i][j]);
            }

            Sequence_Map[Element_Sequence[i][j]].push_back(i);
        }

    for(int i=0;i<Sequence_Pattern_Key.size();i++)
    {
        Sequence_Element_Composition.push_back(std::vector<int>(Observed_Elements.size(),0));

        for(int j=0;j<Sequence_Pattern_Key[i].size();j++)
            for(int k=0;k<Observed_Elements.size();k++)
                if(Sequence_Pattern_Key[i][j]==Observed_Elements[k])
                    Sequence_Element_Composition[i][k]++;
    }

//    for(int i=0;i<Observed_Elements.size();i++)
//        Sequence_Weighted_Length[Observed_Elements[i]];

    for(int i=0;i<Sequence_Map.size();i++)
    {
        //if(((double)Sequence_Map[i].size()*Sequence_Weighted_Length[i]/Total_Weighted_Duration)>=Common_Duration)
        if(((double)Sequence_Map[i].size()/(double)Element_Sequence.size())>=Common_Duration)
        {
            //Common_Sequences_Weighted_Length.push_back(Sequence_Weighted_Length[i]);
            Common_Sequences_Key.push_back(Sequence_Pattern_Key[i]);
            Common_Sequences_Map.push_back(Sequence_Map[i]);
            Common_Sequences_Size_Map.back().push_back(Common_Sequences_Key.size()-1);

            Common_Sequences_Edge_Corrected_Total_Length.push_back((double)Element_Sequence.size());
            //Total_Weighted_Duration_Common_Sequences_Edge_Corrected.push_back(Total_Weighted_Duration);
        }
    }

    std::vector<std::vector<int>> Full_Sequence_Map = Sequence_Map;

    int k_lim = 0;

    if(Element_Sequence.size()<Max_Pattern_Length)
        k_lim = Element_Sequence.size();
    else
        k_lim = Max_Pattern_Length;

    for(int k=2;k<=k_lim;k++)
    {
        Sequence_Size_Map.push_back(std::vector<int>());
        Common_Sequences_Size_Map.push_back(std::vector<int>());

        for(int i=0;i<Sequence_Size_Map[k-2].size();i++)
            for(int j=0;j<Observed_Elements.size();j++)
            {
                bool matched = false;

                //double Temp_Weighted_Length;
                std::vector<int> Temp_Full_Map;
                std::vector<int> Temp_Sequence_Map;

                std::vector<int> Temp_Pattern_Key;
                std::vector<int> Temp_Element_Composition;

                for(int m=0;m<Full_Sequence_Map[Sequence_Size_Map[k-2][i]].size();m++)
                    if(Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]+k-1<Element_Sequence.size())
                        for(int n=0;n<Element_Sequence[Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]+k-1].size();n++)
                        {
                        if(Element_Sequence[Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]+k-1][n]==Observed_Elements[j])
                            if(matched)
                            {
                                if(Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]>Temp_Sequence_Map.back()+k-1)
                                    Temp_Sequence_Map.push_back(Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]);

                                Temp_Full_Map.push_back(Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]);
                            }
                            else
                            {
                                matched=true;

                                Temp_Full_Map.push_back(Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]);
                                Temp_Sequence_Map.push_back(Full_Sequence_Map[Sequence_Size_Map[k-2][i]][m]);

                                Temp_Pattern_Key = Sequence_Pattern_Key[Sequence_Size_Map[k-2][i]];
                                Temp_Pattern_Key.push_back(Sequence_Pattern_Key[Observed_Elements[j]][0]);

                                Temp_Element_Composition = Sequence_Element_Composition[Sequence_Size_Map[k-2][i]];
                                Temp_Element_Composition[j]++;

                                //Temp_Weighted_Length = Sequence_Weighted_Length[Sequence_Size_Map[k-2][i]];
                                //Temp_Weighted_Length += Sequence_Weighted_Length[Observed_Elements[j]];
                            }
                        }

                //Condition to be added to pool of known fragments
                if(Temp_Sequence_Map.size()>=Min_Occurences&&(double)Temp_Sequence_Map.size()>=Min_Duration*(double)Element_Sequence.size())
                {
                    bool Pattern_is_Unique = true;

                    for(int m = 0; m < Sequence_Size_Map.back().size();m++)
                    {
                        bool Same_Element_Composition = true;

                        for(int n=0; n<Observed_Elements.size();n++)
                            if(Temp_Element_Composition[n]!=Sequence_Element_Composition[Sequence_Size_Map.back()[m]][n])
                            {
                                Same_Element_Composition = false;
                                break;
                            }

                        if(Same_Element_Composition)
                            if(Patterns_are_Automorphic(Temp_Pattern_Key,Sequence_Pattern_Key[Sequence_Size_Map.back()[m]]))
                            {
                                Pattern_is_Unique = false;
                                break;
                            }
                    }

                    if(Pattern_is_Unique)
                    {
                        Full_Sequence_Map.push_back(Temp_Full_Map);
                        Sequence_Map.push_back(Temp_Sequence_Map);
                        Sequence_Pattern_Key.push_back(Temp_Pattern_Key);
                        Sequence_Element_Composition.push_back(Temp_Element_Composition);

                        //Sequence_Weighted_Length.push_back(Temp_Weighted_Length);

                        Sequence_Size_Map.back().push_back(Sequence_Pattern_Key.size()-1);

                        //Correct duration

                        double Total_Duration_Edge_Corrected = (double)Element_Sequence.size();

                        //Remove front duration of sequence from calc of relative frequency if it's shorter than the sequence and matches the end of it

                        if(Temp_Sequence_Map[0]<Temp_Pattern_Key.size())
                        {
                            bool Edge_Matches = true;
                            double Duration_to_Remove = 0;

                            for(int m = 0; m < Temp_Sequence_Map[0];m++)
                            {
                                bool fuzzy_match = false;
                                for(int n=0; n<Element_Sequence[Temp_Sequence_Map[0]-m-1].size();n++)
                                    if(Element_Sequence[Temp_Sequence_Map[0]-m-1][n]==Temp_Pattern_Key[Temp_Pattern_Key.size()-m-1])
                                    {
                                        fuzzy_match = true;
                                        break;
                                    }

                                if(!fuzzy_match)
                                {
                                    Edge_Matches = false;
                                    break;
                                }

                                //Duration_to_Remove += Element_Weight[Element_Sequence[Temp_Sequence_Map[0]-m-1]];
                                Duration_to_Remove++;
                            }

                            if(Edge_Matches)
                                Total_Duration_Edge_Corrected -=Duration_to_Remove;
                        }

                        //Remove back duration of sequence from calc of relative frequency if it's shorter than the sequence and matches the beginning of it

                        if(Element_Sequence.size()-(Temp_Sequence_Map.back()+Temp_Pattern_Key.size())<Temp_Pattern_Key.size())
                        {
                            bool Edge_Matches = true;
                            double Duration_to_Remove = 0;

                            for(int m = 0; m < Element_Sequence.size()-(Temp_Sequence_Map.back()+Temp_Pattern_Key.size());m++)
                            {
                                bool fuzzy_match = false;
                                for(int n=0; n<Element_Sequence[Temp_Sequence_Map.back()+Temp_Pattern_Key.size()+m].size();n++)
                                    if(Element_Sequence[Temp_Sequence_Map.back()+Temp_Pattern_Key.size()+m][n]==Temp_Pattern_Key[m])
                                    {
                                        fuzzy_match = true;
                                        break;
                                    }

                                if(!fuzzy_match)
                                {
                                    Edge_Matches = false;
                                    break;
                                }

                                //Duration_to_Remove += Element_Weight[Element_Sequence[Temp_Sequence_Map.back()+Temp_Pattern_Key.size()+m]];
                                Duration_to_Remove++;
                            }

                            if(Edge_Matches)
                                Total_Duration_Edge_Corrected -=Duration_to_Remove;
                        }

                        //

                        //double Relative_Duration = ((double)Temp_Sequence_Map.size()*Temp_Weighted_Length/Total_Duration_Edge_Corrected);
                        double Relative_Duration = ((double)Temp_Sequence_Map.size()*(double)Temp_Pattern_Key.size())/Total_Duration_Edge_Corrected;

                        if(Relative_Duration >= Common_Duration)
                        { 
                            //Common_Sequences_Weighted_Length.push_back(Temp_Weighted_Length);
                            Common_Sequences_Key.push_back(Temp_Pattern_Key);
                            Common_Sequences_Map.push_back(Temp_Sequence_Map);
                            Common_Sequences_Size_Map.back().push_back(Common_Sequences_Key.size()-1);
                            Common_Sequences_Edge_Corrected_Total_Length.push_back(Total_Duration_Edge_Corrected);

                            //Total_Weighted_Duration_Common_Sequences_Edge_Corrected.push_back(Total_Duration_Edge_Corrected);

                            if(Relative_Duration==1) //break out, match is perfect
                                goto Full_Match_Found;

                        }
                    }
                }
            }

        if(Sequence_Size_Map[k-1].size()==0)
            break;
    }

//    if(Common_Sequences_Key.size()>0)
//    {
//        std::vector<double> Sequence_Duration;
//        for(int i=0;i<Common_Sequences_Key.size();i++)
//        {
//            Sequence_Duration.push_back(Common_Sequences_Weighted_Length[i]*(double)Common_Sequences_Map[i].size());
//        }

        Full_Match_Found:

        double Best_Duration = 0;
        double Best_Sequence = 0;

        for(int i=0;i<Common_Sequences_Size_Map.size();i++)
        {
            if(Common_Sequences_Size_Map[i].size()>0)
            {
                std::vector<double> Sequence_Duration;

                for(int j=0;j<Common_Sequences_Size_Map[i].size();j++)
                {
                    Sequence_Duration.push_back(((double)Common_Sequences_Map[Common_Sequences_Size_Map[i][j]].size()*(double)Common_Sequences_Key[Common_Sequences_Size_Map[i][j]].size())/
                                                Common_Sequences_Edge_Corrected_Total_Length[Common_Sequences_Size_Map[i][j]]);

                    //Sequence_Duration.push_back(Common_Sequences_Weighted_Length[Common_Sequences_Size_Map[i][j]]
                    //        /Total_Weighted_Duration_Common_Sequences_Edge_Corrected[Common_Sequences_Size_Map[i][j]]
                    //        *(double)Common_Sequences_Map[Common_Sequences_Size_Map[i][j]].size());
                }

                std::vector<double>::iterator Max_Duration = std::max_element(Sequence_Duration.begin(),Sequence_Duration.end());

                int Size_Best_Sequence = std::distance(Sequence_Duration.begin(),Max_Duration);

                if(*Max_Duration-Best_Duration>Max_Switch_Threshold)
                {
                    Best_Duration = *Max_Duration;
                    Best_Sequence = Common_Sequences_Size_Map[i][Size_Best_Sequence];
                }
            }
        }

//        std::vector<double>::iterator Max_Duration = std::max_element(Sequence_Duration.begin(),Sequence_Duration.end());
//        Max_Sequence = std::distance(Sequence_Duration.begin(),Max_Duration);
//    }
//    else
//    {
//        Max_Sequence = 0;
//    }

        Max_Sequence = Best_Sequence;
}

std::vector<int> Pattern_Sequence::Get_Max_Sequence_Key() const
{
    if(Common_Sequences_Key.size()>Max_Sequence)
        return Common_Sequences_Key[Max_Sequence];
    else
        return std::vector<int>();
}

std::vector<int> Pattern_Sequence::Get_Max_Sequence_Map() const
{
    if(Common_Sequences_Map.size()>Max_Sequence)
        return Common_Sequences_Map[Max_Sequence];
    else
        return std::vector<int>();
}

double Pattern_Sequence::Get_Max_Sequence_Duration() const
{
    if(Common_Sequences_Map.size()>Max_Sequence)
        return (double)Common_Sequences_Map[Max_Sequence].size()*(double)Common_Sequences_Key[Max_Sequence].size();
    else
        return 0;
}

double Pattern_Sequence::Get_Max_Sequence_Relative_Duration() const
{
    if(Common_Sequences_Map.size()>Max_Sequence)
        return ((double)Common_Sequences_Map[Max_Sequence].size()*(double)Common_Sequences_Key[Max_Sequence].size())
                /Common_Sequences_Edge_Corrected_Total_Length[Max_Sequence];
    else
        return 0;

}

//double Pattern_Sequence::Get_Total_Duration() const
//{
//    return Total_Weighted_Duration;
//}


std::vector<int> Pattern_Sequence::Get_Matching_Sequence_Map(std::vector<int> Key_to_Match) const
{
    if(Common_Sequences_Size_Map.size()>=Key_to_Match.size())
        for(int i=0;i<Common_Sequences_Size_Map[Key_to_Match.size()-1].size();i++)
        {
            bool all_elements_matched = true;
            for(int j=0;j<Key_to_Match.size();j++)
            {
                if(Common_Sequences_Key[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]][j]!=Key_to_Match[j])
                {
                    all_elements_matched=false;
                    break;
                }
            }

            if(all_elements_matched)
                return Common_Sequences_Map[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]];
        }

    return std::vector<int>();
}

std::pair<double,std::vector<int>> Pattern_Sequence::Get_Matching_Sequence(std::vector<int> Key_to_Match) const //Element 1 is relative duration, Element 2 is Map
{
    if(Common_Sequences_Size_Map.size()>=Key_to_Match.size())
        for(int i=0;i<Common_Sequences_Size_Map[Key_to_Match.size()-1].size();i++)
        {
            bool all_elements_matched = true;
            for(int j=0;j<Key_to_Match.size();j++)
            {
                if(Common_Sequences_Key[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]][j]!=Key_to_Match[j])
                {
                    all_elements_matched=false;
                    break;
                }
            }

            if(all_elements_matched)
                return std::pair<double,std::vector<int>> (((double)Common_Sequences_Map[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]].size()
                                                          *(double)Common_Sequences_Key[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]].size())
                                                            /Common_Sequences_Edge_Corrected_Total_Length[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]]
                                                           ,Common_Sequences_Map[Common_Sequences_Size_Map[Key_to_Match.size()-1][i]]);
        }

    return std::pair<double,std::vector<int>>();
}

bool Patterns_are_Automorphic(std::vector<int> Sequence_1,std::vector<int> Sequence_2)
{
    if(Sequence_1.size()!=Sequence_2.size())
        return false;

    //bool is_Automorphism = false;

    for(int i=0;i<Sequence_2.size();i++)
    {
        bool is_Matched = true;
        for(int j=0;j<Sequence_2.size();j++)
        {
            int test_pt = i+j;

            if(test_pt >= Sequence_2.size())
                test_pt -= Sequence_2.size();

            if(Sequence_1[j]!=Sequence_2[test_pt])
            {
                is_Matched = false;
                break;
            }
        }

        if(is_Matched)
            return true;
    }

    return false;
}
