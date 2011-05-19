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

namespace neuralgas
{

template<typename T, typename S>
T euclidean (const Vector<T>& x, const Vector<T>& y)
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
