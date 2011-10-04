/** 
* \file BinaryAutomata.h
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef BINARYAUTOMATA_H
#define BINARYAUTOMATA_H

#include <math.h>
#include "DataGenerator.h"

namespace neuralgas {

//! \class BinaryAutomata
/*! \brief Binary Automata experiment proposed by Voegtlin
 */
class BinaryAutomata : public DataGenerator<double>
{
      public:
             BinaryAutomata();
             void generate(const int& number){DataGenerator<double>::generate(number);}
             virtual void reset();
      private:
             virtual  Vector<double>* generate();
             virtual  Vector<double>* next(){return generate();}
             short    _last_bit;

};

BinaryAutomata::BinaryAutomata() : DataGenerator<double>(1)
{
 srand( (unsigned)time( NULL ) );
 if ( (rand() % 7) > 3 )
    _last_bit = 0;
 else _last_bit = 1;
                            
 Vector<double> * v = new Vector<double>(1);
 (*v)[0]= _last_bit;
 _data->push_back(v);
}

void BinaryAutomata::reset()
{
 DataGenerator<double>::reset();
 if ( (rand() % 7) > 3 )
  _last_bit = 0;
 else
  _last_bit = 1;
          
 Vector<double> * v = new Vector<double>(1);
 (*v)[0]= _last_bit;
 _data->push_back(v);

}

Vector<double>* BinaryAutomata::generate()
{
 if ( _last_bit )
 {
      if ( (rand() % 10) <= 4 )
         _last_bit = 0;
 }
 else
 {
      if ( (rand() % 10) <= 3 )
         _last_bit = 1;   
 }
  Vector<double> * v = new Vector<double>(1);
 (*v)[0]= _last_bit;
 return v;

}

} // namespace neuralgas 

#endif
