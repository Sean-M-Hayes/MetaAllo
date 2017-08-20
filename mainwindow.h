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
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTableWidget>

//#include "allometric_metacommunity.h"
#include "parameter_generator_functions.h"
// <string> & "food_web_methods.h" from "parameter_generator_functions.h"
// <boost/multi_array.hpp> & <vector> from "food_web_methods.h"
// boost_matrix definition from "forward_declarations.h";

//#include <string>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_Simulate_clicked();

    void on_Species_Input_valueChanged(int arg1);

    void on_Patches_Input_valueChanged(int arg1);

    void on_Scale_by_Size_Iterates_Input_valueChanged(int arg1);

    void on_actionSave_triggered();

    void on_actionLoad_triggered();

    void on_Write_Input_Prompt_clicked();

    //void on_Food_Web_Read_Input_Prompt_clicked();

    void on_Food_Web_Read_Input_Prompt_clicked();

    void on_Multi_Food_Web_Spin_valueChanged(int arg1);

    void on_Food_Web_Multi_Input_itemChanged(QTableWidgetItem *item);

    void on_Food_Web_Add_clicked();

    void on_Food_Web_Remove_clicked();

    void on_Spatial_Structure_Read_Input_Prompt_clicked();

    void on_Multi_Spatial_Structure_Spin_valueChanged(int arg1);

    void on_Spatial_Structure_Multi_Input_itemChanged(QTableWidgetItem *item);

    void on_Spatial_Structure_Add_clicked();

    void on_Spatial_Structure_Remove_clicked();

private:
    Ui::MainWindow *ui;

    Parameter_Generators Input_Cache;

    std::vector<boost_matrix> Multi_Food_Web;
    std::vector<boost_matrix> Multi_Spatial_Structure;

    void Food_Web_Method();
    void Spatial_Structure_Method();
    void Body_Mass_Method();
    void Metabolism_Method();
    void Functional_Response_Method();
    void Dispersal_Method();
    void Carrying_Capacity_Method();
    void Initial_Abundance_Method();
    void Extinction_Threshold_Method();

};

#endif // MAINWINDOW_H
