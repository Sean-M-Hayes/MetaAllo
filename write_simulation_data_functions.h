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
#ifndef WRITE_SIMULATION_DATA_FUNCTIONS_H
#define WRITE_SIMULATION_DATA_FUNCTIONS_H

#include "allometric_metacommunity.h"

void Write_Timeseries(std::string Write_Path, int nSpecies, const boost_matrix &Time_Series);

void Start_Time_series(std::string Write_File, int nSpecies, int nPatches);

void Write_or_Append_Time_series(std::string Write_File, int nSpecies, const boost_matrix &Time_Series);

void Write_Food_Web_Summary_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output, int n_Food_Webs);

void Write_Spatial_Structure_Summary_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output, int n_Spatial_Structures);

void Write_Species_Output(std::string Write_Path, bool Feasible_Only, const std::vector<Sim_Set_Output> &Output);

void Write_Community_Output(std::string Write_Path, bool Feasible_Only, const std::vector<Sim_Set_Output> &Output);

void Write_Community_Parameters(std::string Write_Path, bool Feasible_Only, const std::vector<Sim_Set_Output> &Output, int n_Food_Webs, int n_Spatial_Structures);

void Write_Set_Summary_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output);

void Write_Patch_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output);

void Write_Population_Output(std::string Write_Path, const std::vector<Sim_Set_Output> &Output);

void Write_Population_Cycle_Timeseries(std::string Write_Path, double Simulation_Resolution, const Dynamics &Output);

void Write_Population_Cycle_Timeseries(std::string Write_Path, double Simulation_Resolution, const std::vector<Sim_Set_Output> &Output);

#endif // WRITE_SIMULATION_DATA_FUNCTIONS_H
