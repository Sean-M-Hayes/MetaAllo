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
#include "lsoda_link.h"

boost_matrix LSODA_Solve(std::vector< double > x ,
                    void (*Simulation_Function)(int*,double*,double[],double[]),
                    void (*Jacobian_Function)(), double T_Starting,
                    double Simulation_Length,
                    double Simulation_Resolution,
                    double Extinction_Threshold,
                    int nSpecies)
{

    //double Extinction_Threshold = 1e-10;

    int NEQ = x.size();
    double Y[NEQ];

    for(int i = 0;i<NEQ;i++)
        Y[i]=x[i];

    bool Species_Check[nSpecies];

    std::vector< double > obs_times;

    std::vector< double > temp_store(NEQ,0);

    std::vector< std::vector<double> > obs_store;

    //bool Feas_test = true;

    int ITOL = 1;
    double ATOL = 1e-6;

    int ITASK = 1;
    int ISTATE = 1;
    int IOPT = 0;

    int JT = 2;

    double RTOL = 1e-6;

    int LIW = 20+NEQ;
    int IWORK[LIW];

    int LRW = 22+NEQ*std::max(16,NEQ+9);
    double RWORK[LRW];

    //Optional Solver Inputs (IOPT must = 1)
    //0 = use default value

//    //IXPR
//    IWORK[4]=0;
//    //MXSTEP
//    IWORK[5]=500;
//    //MXHNIL
//    IWORK[6]=0;
//    //MXORDN
//    IWORK[7]=0;
//    //MXORDS
//    IWORK[8]=0;

//    IWORK[9]=0;

//    //H0
//    RWORK[4]=0;
//    //HMAX
//    RWORK[5]=0;
//    //HMIN
//    RWORK[6]=0;

//    RWORK[7]=0;
//    RWORK[8]=0;
//    RWORK[9]=0;

    //Starts at t=0, calls dlsoda_ for observations at t=Simulation_Resolution,Simulation_Resolution*2,...,Simulation_Length

    double T = 0;

    for(double TOUT = Simulation_Resolution;TOUT<=Simulation_Length;TOUT=TOUT+Simulation_Resolution)
    {
        dlsoda_(Simulation_Function,
                &NEQ,
                Y,
                &T,
                &TOUT,
                &ITOL,
                &RTOL,
                &ATOL,
                &ITASK,
                &ISTATE,
                &IOPT,
                RWORK,
                &LRW,
                IWORK,
                &LIW,
                Jacobian_Function,
                &JT);

        obs_times.push_back(T);

        int n_Non_Feas = 0;

        for(int i=0;i<nSpecies;i++)
        {
            temp_store[i] = Y[i];

            if(Y[i]<=Extinction_Threshold)
            {
                n_Non_Feas++;
                Species_Check[i]=false;
            }
            else
                Species_Check[i]=true;
        }

        for(int i = nSpecies;i<NEQ;i++)
        {
            temp_store[i] = Y[i];

            if(Y[i]>Extinction_Threshold&&!Species_Check[i%nSpecies])
            {
                n_Non_Feas--;
                Species_Check[i%nSpecies]=true;
            }
        }

        obs_store.push_back(temp_store);

        if(ISTATE<0||n_Non_Feas>0)
            break;
    }

    boost_matrix sim_ana (boost::extents[obs_times.size()][NEQ+1]);

    for(int i = 0;i<obs_times.size();i++)
    {
        sim_ana[i][0] = obs_times[i]+T_Starting;

        for(int i2 = 0;i2<NEQ;i2++)
            sim_ana[i][i2+1] = obs_store[i][i2];
    }

    return sim_ana;
}

bool LSODA_Solve(std::vector< double > x ,
                    void (*Simulation_Function)(int*,double*,double[],double[]),
                    void (*Jacobian_Function)(),
                    double T_Starting,
                    double Simulation_Length,
                    double Simulation_Resolution,
                    double Extinction_Threshold,
                    int nSpecies,
                    boost_matrix &Output)
{

    //double Extinction_Threshold = 1e-10;

    int NEQ = x.size();
    double Y[NEQ];

    for(int i = 0;i<NEQ;i++)
        Y[i]=x[i];

    bool Species_Check[nSpecies];

    std::vector< double > obs_times;

    std::vector< double > temp_store(NEQ,0);

    std::vector< std::vector<double> > obs_store;

    //bool Feas_test = true;

    int ITOL = 1;
    double ATOL = 1e-6;

    int ITASK = 1;
    int ISTATE = 1;
    int IOPT = 0;

    int JT = 2;

    double RTOL = 1e-6;

    int LIW = 20+NEQ;
    int IWORK[LIW];

    int LRW = 22+NEQ*std::max(16,NEQ+9);
    double RWORK[LRW];

    //Optional Solver Inputs (IOPT must = 1)
    //0 = use default value

//    //IXPR
//    IWORK[4]=0;
//    //MXSTEP
//    IWORK[5]=500;
//    //MXHNIL
//    IWORK[6]=0;
//    //MXORDN
//    IWORK[7]=0;
//    //MXORDS
//    IWORK[8]=0;

//    IWORK[9]=0;

//    //H0
//    RWORK[4]=0;
//    //HMAX
//    RWORK[5]=0;
//    //HMIN
//    RWORK[6]=0;

//    RWORK[7]=0;
//    RWORK[8]=0;
//    RWORK[9]=0;

    //Starts at t=0, calls dlsoda_ for observations at t=Simulation_Resolution,Simulation_Resolution*2,...,Simulation_Length

    bool Feasible = true;
    double T = 0;

    for(double TOUT = Simulation_Resolution;TOUT<=Simulation_Length;TOUT=TOUT+Simulation_Resolution)
    {
        dlsoda_(Simulation_Function,
                &NEQ,
                Y,
                &T,
                &TOUT,
                &ITOL,
                &RTOL,
                &ATOL,
                &ITASK,
                &ISTATE,
                &IOPT,
                RWORK,
                &LRW,
                IWORK,
                &LIW,
                Jacobian_Function,
                &JT);

        obs_times.push_back(T);

        int n_Non_Feas = 0;

        for(int i=0;i<nSpecies;i++)
        {
            temp_store[i] = Y[i];

            if(Y[i]<=Extinction_Threshold)
            {
                n_Non_Feas++;
                Species_Check[i]=false;
            }
            else
                Species_Check[i]=true;
        }

        for(int i = nSpecies;i<NEQ;i++)
        {
            temp_store[i] = Y[i];

            if(Y[i]>Extinction_Threshold&&!Species_Check[i%nSpecies])
            {
                n_Non_Feas--;
                Species_Check[i%nSpecies]=true;
            }
        }

        obs_store.push_back(temp_store);

        if(ISTATE<0||n_Non_Feas>0)
        {
            Feasible = false;
            break;
        }
    }

    //boost_matrix sim_ana (boost::extents[obs_times.size()][NEQ+1]);

    Output.resize(boost::extents[obs_times.size()][NEQ+1]);

    for(int i = 0;i<obs_times.size();i++)
    {
        Output[i][0] = obs_times[i]+T_Starting;

        for(int i2 = 0;i2<NEQ;i2++)
            Output[i][i2+1] = obs_store[i][i2];
    }

    return Feasible;
}
