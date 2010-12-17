/** 
* \class GNGModul
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef GNGMODUL_H
#define GNGMODUL_H


#include <vector>
#include <limits>
#include <NeuralGas.h>
#include "GNGModulGraph.h"

template<typename T,typename S> class ErrorTesting;

/** \brief This class is a modul offering common functions that are used within nearly
*   every GrowingNeuralGas Algorithms.
*
*   This modul is an abstraction of GrowingNeuralGas Algorithm that is not runnable but
*   comprises moduls that were used within nearly every GrowingNeuralGas Algorithm like 
*   the calculation of the two closest reference vectors, the arc removing etc.
*   It is intended to expand this class whenever components are found to be extractable.
*
*   \param _graphModulptr ptr to the GNGModulGraph
*/

template<typename T,typename S> class GNGModul : public NeuralGas<T,S>
{

 public:
         //cto with size of dimension as input 
                  GNGModul(const int& dim);
         //std dto
                  ~GNGModul();
           virtual void    showGraph()=0;
 protected:
           // class dependent distance function that is used within the winner function
           virtual           T      getDistance(const int&,const int&)=0;
           // func determines for the current time step / data item the two most similar nodes       
           virtual           void   getWinner( int&, int&, const int&);
           // removes all edges that have an age greater than the given value
           virtual           void   rmOldEdges(const float&);
           // removes all nodes from the graph that are not connected
           virtual           void   rmNotConnectedNodes();
           // ptr to the underlying graph structure
           GNGModulGraph<T,S>*      _graphModulptr;
 
 private:         
           // ErrorTesting is defined as friend in order to not having duplicate anything
           friend class ErrorTesting<T,S>;

};

/** \brief cto with size of dimension as input 
*/

template<typename T,typename S> GNGModul<T,S>::GNGModul(const int& dim) : NeuralGas<T,S>(dim)
{}

/** \brief std dto
*/
template<typename T,typename S> GNGModul<T,S>::~GNGModul(){}


/** \brief func determines for the current time step / data item the two most similar nodes
*
*   The metric is used (either the user defined or the preset one to determine the distance
*   from the current data to the current nodes in the network.
*   The two closest nodes are stored in the parameter given by reference.
*  
*   NOTE: It is intented to "outsource" this function and to let the user define the way
*         of determining the winner and refer to it via a function ptr.
*  
*   \param first_winner before func call an arbitrary value, after func call the closest node
*   \param second_winner before func call an arbitrary value, after func call the second closest node
*   \param time is the current time step reflecting the current data to be processed
*/

template<typename T,typename S> void GNGModul<T,S>::getWinner(int& first_winner,int& second_winner, const int& time) 
{
   // init with zero
   T distance      = this->_zero;
   T best_distance = this->_zero; 

   // best_distance set to "infinity"
   best_distance = std::numeric_limits<T>::max();
 
   for (int j = 0; j < _graphModulptr->size(); j++)
   {
    distance = getDistance(time,j);

    if (distance < best_distance)
    {
     second_winner              =       first_winner;
     first_winner               =       j;
     best_distance              =       distance;
    }       
   }
}

/** \brief Removes all edges that have an age greater than the value given by max_age
* 
* \param max_age is the maximal age that is permitted for an edge to have
*/

template<typename T,typename S> void GNGModul<T,S>::rmOldEdges(const float& max_age)
{
    for (int j = 0; j < _graphModulptr->size(); j++)
    {
     for (int k = 0; k < j; k++)
        if ( _graphModulptr->getAge(k,j) > max_age ) // age is greater than max_age
           _graphModulptr->rmEdge(k,j);              // and has to be removedd
    }  
}
/** \brief Removes all nodes from the graph that are not connected
*
*/
template<typename T,typename S> void GNGModul<T,S>::rmNotConnectedNodes()
{
    for (int j = _graphModulptr->size() - 1 ; j >= 0 ; j--)
       if (!(_graphModulptr->isConnected(j)) ) 
           _graphModulptr->rmNode(j);   
}

#endif
