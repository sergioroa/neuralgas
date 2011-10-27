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

//! \brief Approximated Gauss by numerical integration
/*! 
  \param value variable
  \return result
*/
double approx_gauss(const double& value)
{
      double dz = 10E-5;
      double integral = 0.0;
      double z = 0.0;
    
      double tmpvalue=(value >0 )? value : -value; 
 
      while ( integral < tmpvalue  )
      {
            integral+=exp(- pow(z,2))*dz;
            z += dz;
      } 
      return ( value > 0) ? z : -z;
}

//! \brief Box Muller method for obtaining approx. Gauss distribution
/*! 
  
  \param a1 first uniformly distributed random variable
  \param a2 second uniformly distributed random variable
*/void box_muller(double& a1, double& a2)
{
     double rho = sqrt( -2* (log(a1)));
     double phi = 2*M_PI*a2;
     a1 = rho * cos(phi);
     a2 = rho * sin(phi);
}


//! \brief Generation of bi-dimensional normal distribution by applying Box Muller method
/*! 
  
  \param m1 first mean
  \param m2 second mean
  \param _sigma1 first std dev
  \param _sigma2 second std dev
  
  \return normal distributed bi-dimensional vector
*/Vector<double>* normal_distribution(const double& m1, const double& m2, const double& _sigma1, const double& _sigma2)
{
 Vector<double>* v=new Vector<double>(2);

 //double a;
 double a1,a2;
 do
 {
  a1 = double(rand()%10000) / 10000;
 }while(a1==0.0);
 a2 = double(rand()%10000) / 10000;
 
 box_muller(a1,a2); 
 (*v)[0] = a1  * _sigma1 + m1; 
 (*v)[1] = a2  * _sigma2 + m2; 

 /*do
 {
       a = _sqrtPI * ( double(rand()%10000) / 10000) - _sqrtPI / 2;
 }
 while(a>0.88 or a <-0.88); 

 (*v)[0] = _sqrt2 * approx_gauss(a)  * _sigma + m1; 
 
 do
 {
       a = _sqrtPI * (double(rand()%10000) / 10000) - _sqrtPI / 2;
 }
 while(a>0.88 or a <-0.88); 
 
 (*v)[1] = _sqrt2 * approx_gauss(a) * _sigma + m2; 
 */
 return v;
             
}             




} // namespace neuralgas

#endif // MATHHELPERS_H
