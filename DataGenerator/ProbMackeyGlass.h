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
* \file ProbMackeyGlass.h
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef PROBMACKEYGLASS_H
#define PROBMACKEYGLASS_H

#include <math.h>
#include "DataGenerator.h"
#include "MackeyGlass.h"

namespace neuralgas {

//! \class ProbMackeyGlass
/*! \brief A MackeyGlass dataset generated from 2 MackeyGlass sources
 */
class ProbMackeyGlass : public DataGenerator<double>
{
      public:
             ProbMackeyGlass();
             void random(const bool&, const int&);
             void generate(const int& );
             virtual  Vector<double>*  next();
  
      private:
              virtual  Vector<double>* generate();
              std::vector<MackeyGlass*> v_mg;
              int _number;
              bool _random;
              short _r_index;
};

ProbMackeyGlass::ProbMackeyGlass():DataGenerator<double>(1)
{
 _r_index=0;
 _random=true;
 _number=0;
 MackeyGlass* tmp1=new MackeyGlass;
 MackeyGlass* tmp2=new MackeyGlass;
 v_mg.push_back(tmp1);
 v_mg.push_back(tmp2);
 v_mg[0]->setPastTimeSteps(17);
 v_mg[0]->setBoundary(0.4);
 v_mg[0]->setPower(10);
 v_mg[1]->setPastTimeSteps(5);
 v_mg[1]->setBoundary(0.1);
 v_mg[1]->setPower(8);

}

void ProbMackeyGlass::random(const bool& random, const int& number=0)
{
     _random=random;
     _number=number;
}

Vector<double>* ProbMackeyGlass::generate()
{
 return v_mg[_r_index]->next();
}

Vector<double>* ProbMackeyGlass::next()
{
 return generate();
}
 
void ProbMackeyGlass::generate(const int& number)
{
     if (_random)
     {
         for(int i=0; i < number; i++)
         {
          int ran = rand()%3;
          if (ran == 1) 
             _r_index = 1;
          else
              _r_index=0;
          this->_data->push_back(generate());
         }                 
     }
     else
     {
         _r_index=0;
         for(int i=0; i < _number; i++)
         {
          this->_data->push_back(generate());
         }
         _r_index=1;
         for(int i=_number; i < number; i++)
         {
          this->_data->push_back(generate());
         }
     }
}

} // namespace neuralgas

#endif
