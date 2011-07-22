/** 
* \class ProbMackeyGlass
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef PROBMACKEYGLASS_H
#define PROBMACKEYGLASS_H

#include <math.h>
#include "DataGenerator.h"
#include "MackeyGlass.h"

namespace neuralgas {

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
