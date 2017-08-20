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
#ifndef INPUT_READING_FUNCTIONS
#define INPUT_READING_FUNCTIONS

#include "forward_declarations.h"

#include <vector>
#include <QString>
#include <QLineEdit>
#include <QSpinBox>
#include <QTableWidget>

//Loads output from QObject into provided container, returns bool specifying if transfer was successful or not

bool Collect_Input(QSpinBox *Input_Box, int &Output_Int);

bool Collect_Input(QLineEdit *Input_Line, std::vector<int> &Output_Vector);

bool Collect_Input(QLineEdit *Input_Line, std::vector<double> &Output_Vector);

bool Collect_Input(QLineEdit *Input_Line, double &Output_Double);

bool Collect_Input(QLineEdit *Input_Line, int &Output_Int);

bool Collect_Input(QLineEdit *Input_Line, QString &Output_String);

bool Collect_Input(QLineEdit *Input_Line, std::string &Output_String);

bool Collect_Input_Table_Column(int Column_Index, QTableWidget* Input_Table, std::vector<double> &Output_Vector);

bool Collect_Input_Table_Row (int Row_Index, QTableWidget* Input_Table, std::vector<double> &Output_Vector);

bool Collect_Input_Table (QTableWidget* Input_Table, boost_matrix &Output_Matrix);


#endif // INPUT_READING_FUNCTIONS

