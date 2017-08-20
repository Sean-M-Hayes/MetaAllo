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
#ifndef FIND_PATTERNS_H
#define FIND_PATTERNS_H

#include "forward_declarations.h"
#include <vector>

class Pattern_Sequence{

private:
    //std::vector<double> Common_Sequences_Weighted_Length; //

    std::vector<std::vector<int>> Common_Sequences_Key;
    std::vector<std::vector<int>> Common_Sequences_Map;

    std::vector<std::vector<int>> Common_Sequences_Size_Map;

    std::vector<double> Common_Sequences_Edge_Corrected_Total_Length;

    int Max_Sequence;
    //double Total_Weighted_Duration; //
    //std::vector<double> Total_Weighted_Duration_Common_Sequences_Edge_Corrected; //


public:
    Pattern_Sequence(const std::vector<std::vector<int> > &Element_Sequence,
                     const int Min_Occurences,
                     const double Min_Duration,
                     const double Max_Switch_Threshold,
                     const int Max_Pattern_Length);

    std::vector<int> Get_Max_Sequence_Key() const;

    std::vector<int> Get_Max_Sequence_Map() const;

    double Get_Max_Sequence_Duration() const;

    double Get_Max_Sequence_Relative_Duration() const;

    //double Get_Total_Duration() const;

    std::pair<double,std::vector<int>> Get_Matching_Sequence(std::vector<int> Key_to_Match) const;

    std::vector<int> Get_Matching_Sequence_Map(std::vector<int> Key_to_Match) const;

};

bool Patterns_are_Automorphic(std::vector<int> Sequence_1,std::vector<int> Sequence_2);

#endif // FIND_PATTERNS_H
