/** 
* \file MackeyGlass.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef MACKEYGLASS_H
#define MACKEYGLASS_H

#include <math.h>
#include "DataGenerator.h"

namespace neuralgas {

//! \class MackeyGlass
/*! \brief MackeyGlass dataset generated from a chaotic function
 */
class MackeyGlass : public DataGenerator<double>
{
      public:
             MackeyGlass();
             ~MackeyGlass();
             void setPower(const int&);
             void setBoundary(const double&);
             void setPastTimeSteps(const int&);
             void generate(const int& number){DataGenerator<double>::generate(number);}
             virtual  Vector<double>*  next();
      private:
              void  alloc(const int&);
              virtual  Vector<double>* generate();
              int   _pastTimeSteps;
              int   _power;
              double* _memory;
              int    _time;
};

MackeyGlass::MackeyGlass() : DataGenerator<double>(1)
{
 setPastTimeSteps(0);
 setBoundary(0.0);
 setPower(1);
 _time=0;
}

MackeyGlass::~MackeyGlass()
{
 if (_memory != NULL)
    delete _memory;

}

void MackeyGlass::setPower(const int& power)
{_power=power;}

void MackeyGlass::setBoundary(const double& boundary)
{_memory[0]=boundary;}

void MackeyGlass::setPastTimeSteps(const int& pastTimeSteps)
{    
     alloc(pastTimeSteps+1);
     _pastTimeSteps=pastTimeSteps;
}

Vector<double>* MackeyGlass::next()
{return generate();}

void MackeyGlass::alloc(const int& size)
{
     if (_memory!=NULL)
        delete _memory;
     _memory = new double[size];
}

Vector<double>* MackeyGlass::generate()
{
 double dx = 0.0;
 
 int frame;

 if ( _time > _pastTimeSteps )
 {
    frame = ( _time - _pastTimeSteps ) % (_pastTimeSteps + 1 );
    dx = -0.1 * _memory[ _time % (_pastTimeSteps + 1 )  ]+(0.2 * _memory[ frame ])/(1+pow( _memory[ frame ],_power));
 }
 else
 {
    frame = _time ;
    dx = -0.1 * _memory[ _time ]+(0.2 * _memory[0])/(1+pow( _memory[0], _power));
 }
 _memory[  (_time + 1) % (_pastTimeSteps + 1 ) ] = _memory[  _time  % (_pastTimeSteps + 1 ) ] + dx;
 
 _time++;
 
 Vector<double> * v = new Vector<double>(1);
 (*v)[0] = _memory[frame];
 return v;
}

} // namespace neuralgas

#endif
