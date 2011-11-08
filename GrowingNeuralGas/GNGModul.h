/** 
* \file GNGModul.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef GNGMODUL_H
#define GNGMODUL_H


#include <vector>
#include <NeuralGas.h>
#include "GNGModulGraph.h"

namespace neuralgas {

template<typename T,typename S> class ErrorTesting;

/** \class GNGModul
 *  \brief This class is a modul offering common functions that are used within nearly
 *  every GrowingNeuralGas Algorithms.
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
        GNGModul(const unsigned int& dim);
        //std dto
        ~GNGModul();
        virtual void    showGraph()=0;
        //set stopping value (maximal number of epochs)
        void setMaxEpochs (unsigned int);
protected:
        // removes all edges that have an age greater than the given value
        virtual           void   rmOldEdges(const unsigned int&);
        // removes all nodes from the graph that are not connected
        virtual           bool   rmNotConnectedNodes();
        // ptr to the underlying graph structure
        GNGModulGraph<T,S>*      _graphModulptr;
 
private:         
        // ErrorTesting is defined as friend in order to not having duplicate anything
        friend class ErrorTesting<T,S>;
protected:
        //maximal number of epochs (used as stopping criterion)
        unsigned int max_epochs;

};

/** \brief cto with size of dimension as input 
*/

template<typename T,typename S> GNGModul<T,S>::GNGModul(const unsigned int& dim) : NeuralGas<T,S>(dim)
{
  max_epochs = 1;
}

/** \brief std dto
*/
template<typename T,typename S> GNGModul<T,S>::~GNGModul(){}

/** \brief set stopping value (maximal number of epochs)
 *  \param value for maximum epochs
 */
template<typename T,typename S> void GNGModul<T,S>::setMaxEpochs (unsigned int value)
{
  max_epochs = value;
}

/** \brief Removes all edges that have an age greater than the value given by max_age
* 
* \param max_age is the maximal age that is permitted for an edge to have
*/

template<typename T,typename S> void GNGModul<T,S>::rmOldEdges(const unsigned int& max_age)
{
    for (unsigned int j = 0; j < _graphModulptr->size(); j++)
    {
     for (unsigned int k = 0; k < j; k++)
        if ( _graphModulptr->getAge(k,j) > max_age ) // age is greater than max_age
           _graphModulptr->rmEdge(k,j);              // and has to be removedd
    }  
}
/** \brief Removes all nodes from the graph that are not connected
*
*/
template<typename T,typename S> bool GNGModul<T,S>::rmNotConnectedNodes()
{
    assert (_graphModulptr->size() > 0);
    unsigned int original_graph_size = _graphModulptr->size();
    for (unsigned int j = 0; j < _graphModulptr->size(); j++)
       if (!(_graphModulptr->isConnected(j)) && _graphModulptr->size() > 2 )
       {
           _graphModulptr->rmNode(j);
	   j--;
       }

    return original_graph_size > _graphModulptr->size();
}

} // namespace neuralgas

#endif
