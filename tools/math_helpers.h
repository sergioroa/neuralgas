/** 
* \file math_helpers.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef MATHHELPERS_H
#define MATHHELPERS_H

#include <math.h>
#include <Graphs/Vector.h>

namespace neuralgas {

float approx_gauss(const float& value)
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

void box_muller(float& a1, float& a2)
{
     float rho = sqrt( -2* (log(a1)));
     float phi = 2*M_PI*a2;
     a1 = rho * cos(phi);
     a2 = rho * sin(phi);
}


Vector<float>* normal_distribution(const float& m1, const float& m2, const float& _sigma1, const float& _sigma2)
{
 Vector<float>* v=new Vector<float>(2);

 //float a;
 float a1,a2;
 do
 {
  a1 = float(rand()%10000) / 10000;
 }while(a1==0.0);
 a2 = float(rand()%10000) / 10000;
 
 box_muller(a1,a2); 
 (*v)[0] = a1  * _sigma1 + m1; 
 (*v)[1] = a2  * _sigma2 + m2; 

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




} // namespace neuralgas

#endif // MATHHELPERS_H
