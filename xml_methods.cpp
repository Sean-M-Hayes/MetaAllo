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
#include "xml_methods.h"

void Read_XML_Matrix(QXmlStreamReader* xmlReader, QTableWidget* Output_Table)
{
    int i = 0;
    while(!xmlReader->hasError()&&i<Output_Table->rowCount())
    {
        xmlReader->readNext();
        if(xmlReader->isStartElement()&&xmlReader->name()=="row")
        {
            int i2=0;
            while(!xmlReader->hasError()&&i2<Output_Table->columnCount())
            {
                xmlReader->readNext();
                if(xmlReader->isStartElement()&&xmlReader->name()=="col")
                {
                    Output_Table->setItem(i,i2,new QTableWidgetItem(xmlReader->readElementText()));

                    i2++;
                }
                if(xmlReader->isEndElement()&&xmlReader->name()!="col")
                    break;
            }

            i++;
        }
        if(xmlReader->isEndElement()&&!(xmlReader->name()=="row"||xmlReader->name()=="col"))
            break;
    }
}

void Read_XML_Matrix(QXmlStreamReader* xmlReader, boost_matrix &Output)
{
    int i = 0;
    while(!xmlReader->hasError())
    {
        xmlReader->readNext();
        if(xmlReader->isStartElement()&&xmlReader->name()=="row")
        {
            if(i>=Output.shape()[0])
                Output.resize(boost::extents[i+1][Output.shape()[1]]);

            int j=0;
            while(!xmlReader->hasError())
            {
                xmlReader->readNext();
                if(xmlReader->isStartElement()&&xmlReader->name()=="col")
                {
                    if(j>=Output.shape()[1])
                        Output.resize(boost::extents[Output.shape()[0]][j+1]);

                    Output[i][j] = xmlReader->readElementText().toDouble();

                    j++;
                }
                if(xmlReader->isEndElement()&&xmlReader->name()!="col")
                    break;
            }

            i++;
        }
        if(xmlReader->isEndElement()&&!(xmlReader->name()=="row"||xmlReader->name()=="col"))
            break;
    }
}

void Read_XML_Column(QXmlStreamReader* xmlReader, QTableWidget* Output_Table,int Column)
{
    int i = 0;
    while(!xmlReader->hasError()&&i<Output_Table->rowCount())
    {
        xmlReader->readNext();
        if(xmlReader->isStartElement()&&xmlReader->name()=="row")
        {
            Output_Table->setItem(i,Column,new QTableWidgetItem(xmlReader->readElementText()));

            i++;
        }
        if(xmlReader->isEndElement()&&xmlReader->name()!="row")
            break;
    }
}
