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

#include "input_reading_functions.h"

bool Collect_Input(QSpinBox *Input_Box, int &Output_Int)
{
    QString read = Input_Box->text();
    if(read.isEmpty())
        return false;
    else
    {
        Output_Int=read.toInt();
        return true;
    }
//        if(read.toInt()>0)
//        {
//            Output_Int=read.toInt();
//            return true;
//        }
//        else
//            return false;
}

bool Collect_Input(QLineEdit *Input_Line, std::vector<int> &Output_Vector)
{
    QString read = Input_Line->text();
    QRegExp rx(",");
    QStringList read_list = read.split(rx, QString::SkipEmptyParts);

    for(int i=0;i<read_list.size();i++)
    {
        int test = read_list.at(i).toInt();
            Output_Vector.push_back(test);
    }

    if(Output_Vector.size()>0)
        return true;
    else
        return false;

}

bool Collect_Input(QLineEdit *Input_Line, std::vector<double> &Output_Vector)
{
    QString read = Input_Line->text();
    QRegExp rx(",");
    QStringList read_list = read.split(rx, QString::SkipEmptyParts);

    for(int i=0;i<read_list.size();i++)
    {
        double test = read_list.at(i).toDouble();
            Output_Vector.push_back(test);
    }

    if(Output_Vector.size()>0)
        return true;
    else
        return false;

}

bool Collect_Input(QLineEdit *Input_Line, double &Output_Double)
{
    QString read = Input_Line->text();

    if(read.isEmpty())
        return false;
    else
    {
        Output_Double = read.toDouble();

        return true;
    }

//        if(read.toDouble()>=0)
//        {
//            Output_Double = read.toDouble();

//            return true;
//        }
//        else
//            return false;
}

bool Collect_Input(QLineEdit *Input_Line, int &Output_Int)
{
    QString read = Input_Line->text();

    if(read.isEmpty())
        return false;
    else
    {
        Output_Int = read.toInt();
        return true;
    }
//        if(read.toInt()>0)
//        {
//            Output_Int = read.toInt();
//            return true;
//        }
//        else
//            return false;
}

bool Collect_Input(QLineEdit *Input_Line, QString &Output_String)
{
    QString read = Input_Line->text();

    if(read.isEmpty())
        return false;
    else
    {
        Output_String = read;
        return true;
    }
}

bool Collect_Input(QLineEdit *Input_Line, std::string &Output_String)
{
    QString read = Input_Line->text();

    if(read.isEmpty())
        return false;
    else
    {
        Output_String = read.toStdString();
        return true;
    }
}

bool Collect_Input_Table_Column(int Column_Index, QTableWidget* Input_Table, std::vector<double> &Output_Vector)
{
    Output_Vector.clear();

    if(Input_Table->rowCount()>0)
    {
        for(int i=0;i<Input_Table->rowCount();i++)
        {
            QTableWidgetItem * Input_Test = Input_Table->item(i,Column_Index);
            if(Input_Test==0)
            {
                Output_Vector.push_back(0);
            }
            else
            {
                Output_Vector.push_back(Input_Test->text().toDouble());
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool Collect_Input_Table_Row (int Row_Index, QTableWidget* Input_Table, std::vector<double> &Output_Vector)
{
    Output_Vector.clear();

    if(Input_Table->columnCount()>0)
    {
        for(int i=0;i<Input_Table->columnCount();i++)
        {
            QTableWidgetItem * Input_Test = Input_Table->item(Row_Index,i);
            if(Input_Test==0)
            {
                Output_Vector.push_back(0);
            }
            else
            {
                Output_Vector.push_back(Input_Test->text().toDouble());
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}

bool Collect_Input_Table (QTableWidget* Input_Table, boost_matrix &Output_Matrix)
{
    if(Input_Table->rowCount()>0&&Input_Table->columnCount()>0)
    {
        Output_Matrix.resize(boost::extents[Input_Table->rowCount()][Input_Table->columnCount()]);

        for(int i=0;i<Input_Table->rowCount();i++)
        {
            for(int i2=0;i2<Input_Table->columnCount();i2++)
            {
                QTableWidgetItem * Input_Test = Input_Table->item(i,i2);

                if(Input_Test==0)
                {
                    Output_Matrix[i][i2]=0;
                }
                else
                {
                    Output_Matrix[i][i2] = Input_Test->text().toDouble();
                }
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}
