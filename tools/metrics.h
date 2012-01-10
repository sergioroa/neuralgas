/*
 *   This file is part of NeuralGas.
 *
 *   NeuralGas is free software: you can redistribute it and/or modify it
 *   under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   NeuralGas is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with NeuralGas.  If not, see <http://www.gnu.org/licenses/>.
 */

/** 
* \file metrics.h
* \author Sergio Roa
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/
#ifndef METRICS_H
#define METRICS_H

#include <Graphs/Vector.h>
#include <math.h>

namespace neuralgas
{

/**
 * \enum _meandistance_mode
 * \brief for using different mechanisms for calculating mean distances
 */
enum _meandistance_mode {harmonic, /**< calculate harmonic mean distances  */
			 arithmetic /**< calculate arithmetic mean distances */
};

//! \brief The euclidean distance function
/*! 
  
  \param x first vector
  \param y second vector
  
  \return distance
*/template<typename T, typename S>
inline T euclidean (const Vector<T>& x, const Vector<T>& y)
{
	T result;
	T value;
	Vector<T> z = x - y;
	
	//result =  _zero;
	result = 0;
	unsigned int tsize = z.size();
	
	for (unsigned int i=0; i < tsize; i++)
	{
		value  = (z[i]*z[i]);
		result+=  value;
	}
	return T(sqrt(result));
}

} //namespace neuralgas

#endif
