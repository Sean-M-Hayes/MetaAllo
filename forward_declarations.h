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
#ifndef FORWARD_DECLARATIONS
#define FORWARD_DECLARATIONS

#include <boost/multi_array.hpp>

typedef boost::multi_array<double, 2> boost_matrix;
typedef boost_matrix::array_view<2>::type boost_matrix_2d_view;
typedef boost_matrix::array_view<1>::type boost_matrix_1d_view;

#endif // FORWARD_DECLARATIONS

