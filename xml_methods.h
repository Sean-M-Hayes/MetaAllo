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
#ifndef XML_METHODS
#define XML_METHODS

#include <QXmlStreamWriter>
#include <QTableWidget>
#include "forward_declarations.h"

void Read_XML_Matrix(QXmlStreamReader *xmlReader, QTableWidget *Output_Table);

void Read_XML_Matrix(QXmlStreamReader* xmlReader, boost_matrix &Output);

void Read_XML_Column(QXmlStreamReader* xmlReader, QTableWidget* Output_Table,int Column);

#endif // XML_METHODS

