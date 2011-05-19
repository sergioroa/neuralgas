/** 
* \class NoisyAutomata
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef NOISYAUTOMATA_H
#define NOISYAUTOMATA_H

#include "DataGenerator.h"
#include <tools/math_helpers.h>

namespace neuralgas {

class NoisyAutomata : public DataGenerator<float>
{
      public:
             NoisyAutomata();
             void generate(const int& number){DataGenerator<float>::generate(number);}
             virtual void reset();
             void setSigma(const float&);
             void setTransProb(const float&);
             void setState(const int&);
      private:
             virtual  Vector<float>* generate();
             virtual  Vector<float>* next(){return generate();}
             // Vector<float>* normal_distribution(const float&, const float&);
             // float    approx_gauss(const float&);
             // void     box_muller(float&,float&);
             short    _state;
             float    _sigma;
             float    _transProb;
      
};

NoisyAutomata::NoisyAutomata() : DataGenerator<float>(2)
{
 srand( (unsigned)time( NULL ) );
 _state=0;
 _transProb=0.5;
 _sigma=1;
}

void NoisyAutomata::reset()
{
 DataGenerator<float>::reset();
 _state = 0;
 _transProb=0.5;
 _sigma=1;

}

void NoisyAutomata::setState(const int& state)
{
 _state = state;
}

void NoisyAutomata::setSigma(const float& sigma)
{
 _sigma=sigma;
}

void NoisyAutomata::setTransProb(const float& transProb)
{
 _transProb=transProb;
}

Vector<float>* NoisyAutomata::generate()
{
 float prob  =float (rand() % 1000) /1000;
 // states are coded as follows
 // 0 : ba
 // 1 : ab
 // 2 : ca
 // 3 : ac

 if ( _state==0)
 { 
   if (prob <= _transProb)
      {
      
            _state=1;
            return normal_distribution(1,0,_sigma,_sigma);
      }
   else 
      {
            
            _state=3;
            
            return normal_distribution(0,1,_sigma,_sigma);
      }
 }
 else if( _state==1)
 {
   _state=0;
   return normal_distribution(0,0,_sigma,_sigma);
 }
 else if (_state==2)
 {
   if (prob <= _transProb)
      {
            _state=3;
            return normal_distribution(0,1,_sigma,_sigma);
      }
   else 
      {

            _state=1;
            return normal_distribution(1,0,_sigma,_sigma);
      } 
 }
 else if(_state==3)
 {
    _state=2;
    return normal_distribution(0,0,_sigma,_sigma);
 }
 return 0;
}

} // namespace neuralgas

#endif
