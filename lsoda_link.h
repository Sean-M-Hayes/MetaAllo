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
#ifndef LSODA_LINK_H
#define LSODA_LINK_H

//#include <boost/multi_array.hpp>

//typedef boost::multi_array<double, 2> boost_matrix;

#include "forward_declarations.h"

extern "C"
{
    void dlsoda_(void (*F)(int*,double*,double[],double[]),
                int *NEQ,
                double *Y,
                double *T,
                double *TOUT,
                int *ITOL,
                double *RTOL,
                double *ATOL,
                int *ITASK,
                int *ISTATE,
                int *IOPT,
                double *RWORK,
                int *LRW,
                int *IWORK,
                int *LIW,
                void (*JAC)(),
                int *JT);


}

boost_matrix LSODA_Solve( std::vector< double > x ,
                    void (*Simulation_Function)(int*,double*,double[],double[]),
                    void (*Jacobian_Function)(),
                    double T_Starting,
                    double Simulation_Length,
                    double Simulation_Resolution,
                    double Extinction_Threshold,
                    int nSpecies );

bool LSODA_Solve( std::vector< double > x ,
                    void (*Simulation_Function)(int*,double*,double[],double[]),
                    void (*Jacobian_Function)(),
                    double T_Starting,
                    double Simulation_Length,
                    double Simulation_Resolution,
                    double Extinction_Threshold,
                    int nSpecies,
                    boost_matrix &Output);



#endif // LSODA_LINK_H

