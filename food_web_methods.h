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
#ifndef FOOD_WEB_METHODS_H
#define FOOD_WEB_METHODS_H

//

#include "forward_declarations.h"

#include <vector>

//typedef boost::multi_array<double, 2> boost_matrix;

struct Food_Web_Metadata{

    boost_matrix Food_Web;

    double H;
    double AMI;
    double Hc;

    double Connectance;

    //std::vector< double > H_in;
    //std::vector< double > H_out;
    std::vector< double > N_in;
    std::vector< double > N_out;

    std::vector< double > is_Producer;
    std::vector< double > is_Consumer;

    //std::vector< double > Trophic_Position;
    std::vector< double > Min_Trophic_Position;
    std::vector< double > Avg_Trophic_Position;
    std::vector< double > Max_Trophic_Position;
    std::vector<std::vector<int>> Trophic_Positions;

    std::vector< double > Generality;
    std::vector< double > Positional_Index;

    Food_Web_Metadata();
    Food_Web_Metadata(const boost_matrix &Food_Web);

    ~Food_Web_Metadata(){}

    void operator()(const Food_Web_Metadata &New_Values);

    bool is_Equal(const Food_Web_Metadata &Compare_To);

};

double Log2( double n );

boost_matrix Niche_Model(int nSpecies,double Connectance);

Food_Web_Metadata Calculate_Food_Web_Metadata (const boost_matrix &Food_Web);

std::vector<std::vector<int>> Calculate_Food_Chains (const boost_matrix &Food_Web);



#endif // FOOD_WEB_METHODS_H
