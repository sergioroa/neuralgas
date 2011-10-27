/** 
* \file NoisyAutomata.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef NOISYAUTOMATA_H
#define NOISYAUTOMATA_H

#include "DataGenerator.h"
#include <tools/math_helpers.h>
#include <fstream>

namespace neuralgas {

//! \enum noisyautomata_encoding
/*! \brief these are the possible modes for storing a dataset in CrySSMEx format
 */
enum noisyautomata_encoding { vectorial, symbolic };


//! \class NoisyAutomata
/*! \brief Generation of noisy automata datasets from 3 Gaussian distributions as
     proposed by Hammer et al.
 */
class NoisyAutomata : public DataGenerator<double>
{
      public:
             NoisyAutomata();
             void generate(const int& number);
             virtual void reset();
             void setSigma(const double&);
             void setTransProb(const double&);
             void setState(const int&);
             void openCrySSMExFile (const char*);
             void setInputFormat (unsigned int);
             void setOutputFormat (unsigned int);
      private:
             virtual  Vector<double>* generate();
             virtual  Vector<double>* next(){return generate();}
             short    _state;
	     /// std dev
             double    _sigma;
             double    _transProb;
             std::ofstream crySSMExFile;
             unsigned int input_format;
             unsigned int output_format;
      
};

NoisyAutomata::NoisyAutomata() : DataGenerator<double>(2)
{
 srand( (unsigned)time( NULL ) );
 _state=0;
 _transProb=0.5;
 _sigma=1;
 input_format = symbolic;
 output_format = symbolic;
}

void NoisyAutomata::reset()
{
 DataGenerator<double>::reset();
 _state = 0;
 _transProb=0.5;
 _sigma=1;

}

void NoisyAutomata::setState(const int& state)
{
 _state = state;
}

void NoisyAutomata::setSigma(const double& sigma)
{
 _sigma=sigma;
}

void NoisyAutomata::setTransProb(const double& transProb)
{
 _transProb=transProb;
}

void NoisyAutomata::setInputFormat (unsigned int format)
{
  input_format = format;
}

void NoisyAutomata::setOutputFormat (unsigned int format)
{
  output_format = format;
}

void NoisyAutomata::openCrySSMExFile (const char* filename = "cryssmexdata.cry")
{
  crySSMExFile.open(filename);
}


void NoisyAutomata::generate(const int& number)
{
  if (crySSMExFile.is_open())
  {
    crySSMExFile << "# input dim" << std::endl;
    if (input_format == symbolic)
      crySSMExFile << 0 << std::endl;
    else if (input_format == vectorial)
      crySSMExFile << _dim << std::endl;

    crySSMExFile << "# state dim" << std::endl;
    crySSMExFile << _dim << std::endl;
    crySSMExFile << "# output dim" << std::endl;
    if (output_format == symbolic)
      crySSMExFile << 0 << std::endl;
    else if (output_format == vectorial)
      crySSMExFile << _dim << std::endl;

    if (input_format == symbolic)
    {
      crySSMExFile << "# nr input symbols" << std::endl << 3 << std::endl;
      crySSMExFile << "# examples" << std::endl << "a b c" << std::endl;
    }
    else if (input_format == vectorial)
    {
      // crySSMExFile << "# nr input symbols" << std::endl << 0.0 << std::endl;
      crySSMExFile << "# nr input symbols" << std::endl << 3 << std::endl;
      crySSMExFile << "# examples" << std::endl;
      crySSMExFile << "0.0 0.0 a" << std::endl;
      crySSMExFile << "1.0 0.0 b" << std::endl;
      crySSMExFile << "0.0 1.0 c" << std::endl;
    }
    if (output_format == symbolic)
    {
      crySSMExFile << "# nr output symbols" << std::endl;
      crySSMExFile << 3 << std::endl;
      crySSMExFile << "# examples" << std::endl << "a b c" << std::endl;
    }
    else if (output_format == vectorial)
    {
      crySSMExFile << "# nr output symbols" << std::endl;
      // crySSMExFile << 0.0 << std::endl;
      crySSMExFile << 4 << std::endl;
      crySSMExFile << "1.0 0.0 ab" << std::endl;
      crySSMExFile << "-1.0 0.0 ba" << std::endl;
      crySSMExFile << "0.0 1.0 ac" << std::endl;
      crySSMExFile << "0.0 -1.0 ca" << std::endl;
    }
    
    
  }
  DataGenerator<double>::generate(number);
  if (crySSMExFile.is_open())
    crySSMExFile.close();
}

void write_vector (std::ofstream& file, const Vector<double>& v)
{
  for (unsigned int i=0; i<v.size(); i++)
    file << v[i] << "  ";
}

Vector<double>* NoisyAutomata::generate()
{
 double prob  =double (rand() % 1000) /1000;
 // states are coded as follows
 // 0 : ba
 // 1 : ab
 // 2 : ca
 // 3 : ac

 Vector<double>* previous_item;
 Vector<double>* current_item;

 if (_data->size() == 0)
 {
   previous_item = normal_distribution(0,0,_sigma,_sigma);
   return previous_item; 
 }
 else
   previous_item = _data->back();

 if ( _state==0)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, *previous_item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "a  ";
       
   if (prob < _transProb)
   {
     
     _state=1;
     current_item = normal_distribution(1,0,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, *current_item);
       if (output_format == symbolic)
	 crySSMExFile << "b";
       else if (output_format == vectorial)
	 write_vector (crySSMExFile,  Vector<double>(*current_item - *previous_item));
       crySSMExFile << std::endl;
     }
   }
   else 
   {
     
     _state=3;
     current_item = normal_distribution(0,1,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, *current_item);
       if (output_format == symbolic)
	 crySSMExFile << "c";
       else if (output_format == vectorial)
	 write_vector (crySSMExFile,  Vector<double>(*current_item - *previous_item));
       crySSMExFile << std::endl;
     }
   }
   return current_item;
 }
 else if( _state==1)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, *previous_item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "b  ";
   
   _state=0;
   current_item = normal_distribution(0,0,_sigma,_sigma);
   if (crySSMExFile.is_open())
   {
     write_vector (crySSMExFile, *current_item);
     if (output_format == symbolic)
       crySSMExFile << "a";
     else if (output_format == vectorial)
       write_vector (crySSMExFile,  Vector<double>(*current_item - *previous_item));
     crySSMExFile << std::endl;
   }   
   return current_item;
 }
 else if (_state==2)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, *previous_item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "a  ";
   
   if (prob < _transProb)
   {
     _state=3;
     current_item = normal_distribution(0,1,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, *current_item);
       if (output_format == symbolic)
	 crySSMExFile << "c";
       else if (output_format == vectorial)
	 write_vector (crySSMExFile,  Vector<double>(*current_item - *previous_item));
       crySSMExFile << std::endl;
     }
   }
   else 
   {
     
     _state=1;
     current_item = normal_distribution(1,0,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, *current_item);
       if (output_format == symbolic)
	 crySSMExFile << "b";
       else if (output_format == vectorial)
	 write_vector (crySSMExFile,  Vector<double>(*current_item - *previous_item));
       crySSMExFile << std::endl;
     } 
   } 
   return current_item;
 }
 else if(_state==3)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, *previous_item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "c  ";
   _state=2;
   current_item = normal_distribution(0,0,_sigma,_sigma);
   if (crySSMExFile.is_open())
   {
     write_vector (crySSMExFile, *current_item);
     if (output_format == symbolic)
       crySSMExFile << "a";
     else if (output_format == vectorial)
       write_vector (crySSMExFile,  Vector<double>(*current_item - *previous_item));
     crySSMExFile << std::endl;
   }   
   return current_item;
 }
 return 0;
}

} // namespace neuralgas

#endif
