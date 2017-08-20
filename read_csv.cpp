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
#include "read_csv.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/tokenizer.hpp>

boost_matrix Read_CSV_to_Double(std::string File_Path)
{
    std::ifstream file(File_Path);

    std::vector<std::vector<std::string>> vec;
    std::string line;

    typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;

    while (getline(file,line))
    {
        Tokenizer tok(line);

        vec.push_back(std::vector<std::string>());
        vec.back().assign(tok.begin(),tok.end());
    }

    std::vector<std::vector<double>> vec_double;

    for(int i=0;i<vec.size();i++)
    {
        std::vector<double> to_Build;

        for(int i2=0;i2<vec[i].size();i2++)
        {
            try
            {
                to_Build.push_back(std::stod(vec[i][i2]));
            }
            catch (std::invalid_argument ex)
            {
                break;
            }
        }

        if(to_Build.size()==vec[i].size()-1)
        {
            vec_double.push_back(to_Build);
        }
    }

    boost_matrix to_Output(boost::extents[vec_double.size()][vec_double[0].size()]);

    for(int i=0;i<to_Output.shape()[0];i++)
        for(int i2=0;i2<to_Output.shape()[1];i2++)
            to_Output[i][i2]=vec_double[i][i2];

    return to_Output;
}

bool Read_CSV_Matrices(std::string File_Path, std::vector<boost_matrix>& Output, int force_row_size, int force_col_size)
{
    Output.clear();

    //Assumes elements 0 & 1 are identifier

    std::ifstream file(File_Path);

    std::vector<std::vector<std::string>> vec;
    std::string line;

    typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;

    while (getline(file,line))
    {
        Tokenizer tok(line);

        vec.push_back(std::vector<std::string>());
        vec.back().assign(tok.begin(),tok.end());
    }

    std::vector<std::vector<double>> vec_double;
    std::vector<double> matrix_row_size;

    double row_length = 0;

    for(int i=0;i<vec.size();i++)
    {
        std::vector<double> to_Build;

        for(int i2=0;i2<vec[i].size();i2++)
        {
            try
            {
                to_Build.push_back(std::stod(vec[i][i2]));
            }
            catch (std::invalid_argument ex)
            {
                break;
            }
        }

        if(to_Build.size()>2)
        {
            if(vec_double.size()>0)
                if(vec_double.back()[0]!=to_Build[0]||vec_double.back()[1]!=to_Build[1])
                {
                    matrix_row_size.push_back(row_length);
                    row_length=0;
                }

            vec_double.push_back(to_Build);

            row_length++;
        }
    }

    matrix_row_size.push_back(row_length);

    int vec_row = 0;

    for(int i=0;i<matrix_row_size.size();i++)
    {
        boost_matrix to_Build(boost::extents[force_row_size][force_col_size]);

        for(int j=0;j<to_Build.shape()[0];j++)
        {
            for(int k=0;k<to_Build.shape()[1];k++)
            {

                if(vec_double[vec_row].size()>k+2&&matrix_row_size[i]>j)
                    to_Build[j][k]=vec_double[vec_row][k+2];
                else
                    to_Build[j][k]=0;
            }

            if(matrix_row_size[i]>=j)
                vec_row++;
        }

        Output.push_back(to_Build);
    }

    if(Output.size()>0)
        return true;
    else
        return false;

}
