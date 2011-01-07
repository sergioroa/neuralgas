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

#define PI 3.14159

#include <math.h>
#include "DataGenerator.h"

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
             virtual  Vector<float>* next(){generate();}
             Vector<float>* normal_distribution(const float&, const float&);
             float    approx_gauss(const float&);
             void     box_muller(float&,float&);
             short    _state;
             float    _sigma;
             float    _transProb;
             float    _sqrt2;
             float    _sqrtPI;
      
};

NoisyAutomata::NoisyAutomata() : DataGenerator<float>(2)
{
 srand( (unsigned)time( NULL ) );
 _state=0;
 _transProb=0.5;
 _sigma=1;
 _sqrt2 = sqrt(2);
 _sqrtPI = sqrt(M_PI);
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

float NoisyAutomata::approx_gauss(const float& value)
{
      float dz = 10E-5;
      float integral = 0.0;
      float z = 0.0;
    
      float tmpvalue=(value >0 )? value : -value; 
 
      while ( integral < tmpvalue  )
      {
            integral+=exp(- pow(z,2))*dz;
            z += dz;
      } 
      return ( value > 0) ? z : -z;
}

void NoisyAutomata::box_muller(float& a1, float& a2)
{
     float rho = sqrt( -2* (log(a1)));
     float phi = 2*M_PI*a2;
     a1 = rho * cos(phi);
     a2 = rho * sin(phi);
}

Vector<float>* NoisyAutomata::normal_distribution(const float& m1, const float& m2)
{
 Vector<float>* v=new Vector<float>(2);

 float a;
 float a1,a2;
 do
 {
  a1 = float(rand()%10000) / 10000;
 }while(a1==0.0);
 a2 = float(rand()%10000) / 10000;
 
 box_muller(a1,a2); 
 (*v)[0] = a1  * _sigma + m1; 
 (*v)[1] = a2  * _sigma + m2; 

 /*do
 {
       a = _sqrtPI * ( float(rand()%10000) / 10000) - _sqrtPI / 2;
 }
 while(a>0.88 or a <-0.88); 

 (*v)[0] = _sqrt2 * approx_gauss(a)  * _sigma + m1; 
 
 do
 {
       a = _sqrtPI * (float(rand()%10000) / 10000) - _sqrtPI / 2;
 }
 while(a>0.88 or a <-0.88); 
 
 (*v)[1] = _sqrt2 * approx_gauss(a) * _sigma + m2; 
 */
 return v;
             
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
            return normal_distribution(1,0);
      }
   else 
      {
            
            _state=3;
            
            return normal_distribution(0,1);
      }
 }
 else if( _state==1)
 {
   _state=0;
   return normal_distribution(0,0);
 }
 else if (_state==2)
 {
   if (prob <= _transProb)
      {
            _state=3;
            return normal_distribution(0,1);
      }
   else 
      {

            _state=1;
            return normal_distribution(1,0);
      } 
 }
 else if(_state==3)
 {
    _state=2;
    return normal_distribution(0,0);
 }
 
}

} // namespace neuralgas

#endif
