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
#include <fstream>

namespace neuralgas {

enum noisyautomata_encoding { vectorial, symbolic };

class NoisyAutomata : public DataGenerator<float>
{
      public:
             NoisyAutomata();
             void generate(const int& number);
             virtual void reset();
             void setSigma(const float&);
             void setTransProb(const float&);
             void setState(const int&);
             void openCrySSMExFile (const char*);
             void setInputFormat (unsigned int);
             void setOutputFormat (unsigned int);
      private:
             virtual  Vector<float>* generate();
             virtual  Vector<float>* next(){return generate();}
             // Vector<float>* normal_distribution(const float&, const float&);
             // float    approx_gauss(const float&);
             // void     box_muller(float&,float&);
             short    _state;
             float    _sigma;
             float    _transProb;
             std::ofstream crySSMExFile;
             unsigned int input_format;
             unsigned int output_format;
      
};

NoisyAutomata::NoisyAutomata() : DataGenerator<float>(2)
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
      crySSMExFile << "# nr input symbols" << std::endl << 0.0 << std::endl;
    if (output_format == symbolic)
    {
      crySSMExFile << "# nr output symbols" << std::endl;
      crySSMExFile << 3 << std::endl;
      crySSMExFile << "# examples" << std::endl << "a b c" << std::endl;
    }
    
    
  }
  DataGenerator<float>::generate(number);
  if (crySSMExFile.is_open())
    crySSMExFile.close();
}

void write_vector (std::ofstream& file, const Vector<float>* v)
{
  for (unsigned int i=0; i<v->size(); i++)
    file << v->at(i) << "  ";
}

Vector<float>* NoisyAutomata::generate()
{
 float prob  =float (rand() % 1000) /1000;
 // states are coded as follows
 // 0 : ba
 // 1 : ab
 // 2 : ca
 // 3 : ac

 Vector<float>* item;

 if (_data->size() == 0)
 {
   item = normal_distribution(0,0,_sigma,_sigma);
   return item; 
 }
 else
   item = _data->back();

 if ( _state==0)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "a  ";
       
   if (prob <= _transProb)
   {
     
     _state=1;
     item = normal_distribution(1,0,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, item);
       if (output_format == symbolic)
	 crySSMExFile << "b" << std::endl;
     }
     return item; 
   }
   else 
   {
     
     _state=3;
     item = normal_distribution(0,1,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, item);
       if (output_format == symbolic)
	 crySSMExFile << "c" << std::endl;
     }
     return item;
   }
 }
 else if( _state==1)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "b  ";

   _state=0;
   item = normal_distribution(0,0,_sigma,_sigma);
   if (crySSMExFile.is_open())
   {
     write_vector (crySSMExFile, item);
     if (output_format == symbolic)
       crySSMExFile << "a" << std::endl;
   }   
   return item;
 }
 else if (_state==2)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "a  ";

   if (prob <= _transProb)
   {
     _state=3;
     item = normal_distribution(0,1,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, item);
       if (output_format == symbolic)
	 crySSMExFile << "c" << std::endl;
     }
     return item;
   }
   else 
   {
     
     _state=1;
     item = normal_distribution(1,0,_sigma,_sigma);
     if (crySSMExFile.is_open())
     {
       write_vector (crySSMExFile, item);
       if (output_format == symbolic)
	 crySSMExFile << "b" << std::endl;
     } 
     return item;
   } 
 }
 else if(_state==3)
 {
   if (crySSMExFile.is_open() && input_format == vectorial)
     write_vector (crySSMExFile, item);
   else if (crySSMExFile.is_open() && input_format == symbolic)
     crySSMExFile << "c  ";
   _state=2;
   item = normal_distribution(0,0,_sigma,_sigma);
   if (crySSMExFile.is_open())
   {
     write_vector (crySSMExFile, item);
     if (output_format == symbolic)
       crySSMExFile << "a" << std::endl;
   }   
   return item;
 }
 return 0;
}

} // namespace neuralgas

#endif
