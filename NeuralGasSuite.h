/** 
* \class NeuralGasSuite
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef NEURALGASSUITE_H
#define NEURALGASSUITE_H

#include <vector>
#include "DataGenerator/DataGenerator.h"
#include "GrowingNeuralGas/GNGModul.h"
#include "GrowingNeuralGas/Testing/ErrorTesting.h"
#include "GrowingNeuralGas/MergeGrowingNeuralGas/MGNGAlgorithm.h"

template <typename T,typename S> class NeuralGasSuite
{
 public:
         NeuralGasSuite();
         ~NeuralGasSuite();
         void                 setRefVectors(const int&);
         void                 add(GNGModul<T,S>*);
         void                 rm(GNGModul<T,S>*);
         void                 rm(const int&);
         void                 setDataGenerator(DataGenerator<T>*);
         void                 run();
         std::vector<T>       getErrors(const int&,const int&, const bool&);
         std::vector<T>       getErrors(GNGModul<T,S>*,const int&, const bool& );
         GNGModul<T,S>*        operator[](const int& index){return _algos[index];}
         const GNGModul<T,S>*  operator[](const int& index) const{return _algos[index];}
         const int            size() const;
 private:
         DataGenerator<T>*  _dg;
         ErrorTesting<T,S>*   _et; 
         std::vector<GNGModul<T,S> *> _algos;
};

template <typename T,typename S> NeuralGasSuite<T,S>::NeuralGasSuite()
{
         _et       =          new ErrorTesting<T,S>;
}

template <typename T,typename S> NeuralGasSuite<T,S>::~NeuralGasSuite()
{
         delete _et;
}


template <typename T,typename S> const int NeuralGasSuite<T,S>::size() const
{
         return _algos.size();
}

template <typename T,typename S> void NeuralGasSuite<T,S>::add(GNGModul<T,S>* algo)
{_algos.push_back(algo);}

template <typename T,typename S> void NeuralGasSuite<T,S>::rm(GNGModul<T,S>* algo)
{
         for (int i=0; i <_algos.size();i++)
             if ( *(_algos[i]) == *algo )
             {
                  _algos.erase(_algos.begin() + i);
                  return;
             }
}

template <typename T,typename S> void NeuralGasSuite<T,S>::rm(const int& index)
{_algos.erase(_algos.begin() + index);}

template <typename T,typename S> void NeuralGasSuite<T,S>::setDataGenerator(DataGenerator<T>* dg)
{_dg=dg;}

template <typename T,typename S> void NeuralGasSuite<T,S>::run()
{
         for (int i=0; i <_algos.size();i++)
         {
             std::cout << "Starting run "<<i <<std::endl;
             _algos[i]->setData(_dg->getData());
             _algos[i]->run();
             std::cout << "Run "<<i <<" done"<<std::endl;
         }
}

template <typename T,typename S> std::vector<T> NeuralGasSuite<T,S>::getErrors(GNGModul<T,S>* algo,const int& pastTimeSteps, const bool& random=false )
{
         
         for (int i=0; i <_algos.size();i++)
             if ( *(_algos[i]) == *algo )
             {
                _et->setGNGObject(_algos[i]);   
                break;
             }
         return _et->getErrors(pastTimeSteps,random);      
}

template <typename T,typename S> std::vector<T> NeuralGasSuite<T,S>::getErrors(const int& num_of_algo,const int& pastTimeSteps, const bool& random=false )
{
         _et->setGNGObject(_algos[num_of_algo]);   
         return _et->getErrors(pastTimeSteps,random);      
}

template <typename T,typename S> void NeuralGasSuite<T,S>::setRefVectors(const int& number)
{
         int max_value = _dg->maxRandomValue();
        
         for (int i=0; i <_algos.size();i++)
         {
              _algos[i]->setRefVectors(number,max_value);
         }
}
         
#endif
