/** 
* \class GNGAlgorithm
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/


#ifndef GNGALGORITHM_H
#define GNGALGORITHM_H

#include "GNGGraph.h"
#include <NeuralGas.h>
#include "GNGModul.h"
#include <limits>
#include <fstream>
#include <iostream>

namespace neuralgas {

/** \brief Class implements the algorithm proposed in 
*   "Growing Neural Gas for Temporal Clustering" by Isaac J. Sledge and James M. Keller.
*
* The class is a flexible (not yet as flexible as it is intented to be) implemention of
* the algorithm proposed by the above mentionend authors.
* Since the class is derived from the abstract class NeuralGas it is necessary to define
* the precompiler macro NUM_PARAM reflecting the number of needed parameter params within the
* algorithm.
* Recap : the values are set via the super class func setFuncArray taking a func ptr of the
*  float (*)(const int&) and an index of the parameter to set. The setting is
* params[0] alpha, params[1] beta, params[2] gamma, params[3] delta, params[4] epsilon_b, 
* params[5] epsilon_n, params[6] theta params[7] eta params[8] lambda
* In order to change the number of parameters the precompile definition NUM_PARAM in file 
* NeuralGas.h has to be changed accordingly.
* The int parameter in float (*)(const int&) is the given time step allowing to define
* time / step dependend parameters.
*
* IMPORTANT : USAGE
*
* setFuncArray(index) has to be called for each parameter(#NUM_PARAM)
*
* setData 
*
* maximal random value has to be determined
*
* setRefVectors(num,max_random) inits graph with num ref vectors, their random init weight
*                               and context values are within max_random.
*                               Is the data updated by assigning new data by the setData func
*                               the call of setRefVectors erases the learned graph structure.
*                               Therefore the call is equal to a new init of the algorithm.
*
* run()                         Starts the learning algorithm.
*                               Is some new data assigned by the setData func without calling
*                               setRefVectors then another call
*                               of run() refines, adjusts respectively the graph to this newly
*                               added data.
*
* OPTIONAL:
* setMetric(T (*)(const Vector<T>& a,const Vector<T>& b)) Sets a new metric for distance measurement.
*                                                         Default is L2 euclidean metric.
*                                                         Called with NULL sets the default metric.
* 
*
*   \param _graphptr is a Base_Graph casted pointer to the thereof derived class GNGGraph
*/

template<typename T,typename S> class GNGAlgorithm : public GNGModul<T,S>
{
 public:          // cto initializing the class 
                         GNGAlgorithm(const int& dim);
                  // std dto
                         ~GNGAlgorithm();
                         
                 // runs the proposed algorithm 
                 void    run();
           
                 // sets the number of inital reference vectors
         virtual void    setRefVectors(const int&,const int&);
                
                 void    showGraph(){_graphptr->showGraph();}
                 // algorithmic dependent distance function
                 T getDistance(const Vector<T>&,const int&) const;
 private:        
                 //is a Base_Graph casted pointer to thereof derived class MGNGGraph
                 GNGGraph<T,S>*           _graphptr;
                 //defines the update rule for the in the second index given neighbor by using a given datum
                 void updateNeighbor(const int&,const int&);
                 //defines the update rule for the winner
                 void updateWinner(const int&,const int&);
                    
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S> GNGAlgorithm<T,S>::GNGAlgorithm(const int& dim): GNGModul<T,S>(dim)
{_graphptr=NULL;}

/** \brief std dto
*/
template<typename T,typename S> GNGAlgorithm<T,S>::~GNGAlgorithm()
{
//if (_graphptr!=NULL) 
  //    delete _graphptr;
   _graphptr=NULL;                          
}

/** \brief Sets the number of initial reference vectors.
*
* The function sets the initial number of reference vectors and initializes those
* with random values that are within the range of the given parameter max_value.
* \param num_of_ref_vec is the number of initial reference vectors
* \param max_value is the max value that shall be used for the random init value generation
*/
template<typename T,typename S> void GNGAlgorithm<T,S>::setRefVectors(const int& num_of_ref_vec,const int& max_value)
{
  if (_graphptr!=NULL)
      delete _graphptr;
  if (this->graphptr!=NULL)
     delete this->graphptr;
  
  _graphptr           = new GNGGraph<T,S>(this->getDimension());
  // DANGER DownCast is performed via dynamic_cast
  //_graphptr       = dynamic_cast< MGNGGraph<T,S> * >(this->_graphModulptr);
  this->graphptr      = _graphptr;
  this->_graphModulptr = _graphptr; 
  _graphptr->setMaxRandomValue(max_value);   // sets the max random value for the init of the context vector
  _graphptr->initRandomGraph(num_of_ref_vec,max_value); // creates a Graph object with given size of the 
                                                // vectors and number of ref vectors initilized with 
                                                // random values
}

/** \brief Algorithmic dependent distance function
*
*   This function returns the distance of the given datum and the given node. 
*   The distance is a algorithmic dependent function that is either
*   just the setted metric or a combination thereof.
*   Currently dist  = metric(x_t,w_j) where x_t is the data vector and w_j the node vector,
*
*   \param item datum
*   \param node_index is the node where to the distance shall be determined
*/
template<typename T,typename S> T GNGAlgorithm<T,S>::getDistance(const Vector<T>& item, const int& node_index) const
{
    // dist  = metric(x_t,w_j) instead of metric(x_t,w_j)^2 as proposed in the paper
    // since this accelerates the calculation but does not change the result 
    return metric( item, (*_graphptr)[node_index].weight);
}

/** \brief Defines the update rule for the neighbor given by the second index 
*
*   The update rule depends on the actual datum. With that datum a given
*   topological neighbor is updated by an algorithmic dependent rule.
*   w_j(new) = w_j(old) + epsilon_n * ( x_t - w_j(old) )
*
*   \param time is the data vector that is used for updating 
*   \param node_index is the index of topological neighbor that shall be updated
*/
template<typename T,typename S> void GNGAlgorithm<T,S>::updateNeighbor(const int& time,const int& node_index)
{
 (*_graphptr)[node_index].weight  += this->params[5] * ( (*this)[time]-(*_graphptr)[node_index].weight);
}

/** \brief defines the update rule for the winner
*
*   The update rule depends on the actual datum. With that datum the given
*   winner is updated by an algorithmic dependent rule.
*    w(new) = w(old) + epsilon_b * (x_t - w(old) )
*
*   \param time is the data vector that is used for updating 
*   \param winner is the index of the winner that shall be updated
*/
template<typename T,typename S> void GNGAlgorithm<T,S>::updateWinner(const int& time,const int& winner)
{
 (*_graphptr)[winner].weight  += this->params[4] * ( (*this)[time]-(*_graphptr)[winner].weight);
}


/** \brief Runs the algorithm as proposed in "Incremental Unsupervised Time Series 
*   Analysis using Merge Growing Neural Gas" by Andreakis,Hoyningen-Huene and Beets
*
*   Runs the algorithm where a greater flexibility is implemented by outsourcing the counter
*   update rule and by assigning possibly time depending functions to the parameters.
*
*   NOTE: It is intended to abstract it further by outsourcing further elements, like
*   the learning rule.
*
*/
template<typename T,typename S> void GNGAlgorithm<T,S>::run()
{
 
 if (this->getDimension()>0)
 {
   // number of total steps for the algorithm
   int tsize                       =     this->size();
 
   for(int t = 0; t < tsize; t++)
   {
    // init
    int first_winner               =     1;
    int second_winner              =     0;
 
   //params are defined by user set functions that may depend on the time step
   //params[0] alpha factor by which the two nodes with greatest error are changed after adding a new nod
   //params[1] beta factor by which the error of all nodes is changed
   //params[2] gamma maximal permitted distortion before adding a new reference
   //params[3] alpha_max maximal age for an edge 
   //params[4] epsilon_b weight in winner node update
   //params[5] epsilon_n weight in neighbor node update
   //params[6] lambda number of iterations after a new node is added
   //params[7] w_max maximal number of reference vectors / maximal size of cookbook
  
    for( int j = 0; j < NUM_PARAM; j++)  
    {
      this->params[j] =((*this)._funcArray[j])(t);
    }
    
    this->getWinner(first_winner,second_winner,(*this)[t]);

    T distance = pow(getDistance((*this)[t],first_winner),2);

    if ( _graphptr->size() < this->params[7] &&  distance >= this->params[2] )
    {
     _graphptr->addNode();
 
     int gsize = _graphptr->size()-1;
     //set the position of the new node to be the position of the current data vector
     (*_graphptr)[gsize].weight = (*this)[t];
     _graphptr->setAge(first_winner,gsize,0.0);
     //new node's error is set to the error of the winner node 
     _graphptr->setError(gsize,_graphptr->getError(first_winner) );  
    } 
    else
    {
     int gsize = _graphptr->size()-1;

     std::vector<int> first_winner_neighbors = _graphptr->getNeighbors(first_winner);
   
     for(int j=0; j < first_winner_neighbors.size();j++)
     {
      _graphptr->incAge(first_winner, first_winner_neighbors[j] ); 
      updateNeighbor(t,first_winner_neighbors[j]);
     }

     updateWinner(t,first_winner);      
     _graphptr->setError(gsize,_graphptr->getError(first_winner) + distance  );         
   
     if (first_winner != second_winner)
        _graphptr->setAge(first_winner,second_winner,0.0);

    }

    this->rmOldEdges(this->params[3]);  
    
    this->rmNotConnectedNodes(); 
    
    if( 0<int(this->params[6]) && (t % int(this->params[6])== 0) && _graphptr->size() < this->params[7])
    {
       T max_error                   = this->_zero;
       int max_error_index           = 0;
 
       for (int i=0; i < _graphptr->size(); i++)
       {
           if (_graphptr->getError(i) > max_error )
           {
            max_error               = _graphptr->getError(i);
            max_error_index         = i;
           }
       }
       
       std::vector<int> neighbors   = _graphptr->getNeighbors(max_error_index);
       T max_error_n                = this->_zero;
       int max_error_index_n        = 0;
       
       for (int i=0; i < neighbors.size(); i++)
       {

           if (_graphptr->getError(neighbors[i]) > max_error_n )
           {
            max_error_n             = _graphptr->getError(neighbors[i]);
            max_error_index_n       = i;
           }
       }
       
       _graphptr->addNode();

       int gsize                       = _graphptr->size() - 1;
       (*_graphptr)[gsize].weight      = (*_graphptr)[max_error_index_n].weight*0.5
                                         +(*_graphptr)[max_error_index].weight*0.5;
    
       _graphptr->addEdge(max_error_index,gsize);
       _graphptr->addEdge(max_error_index_n,gsize);

       _graphptr->setError(max_error_index,this->params[0]*_graphptr->getError(max_error_index));
       _graphptr->setError(max_error_index_n,this->params[0]*_graphptr->getError(max_error_index_n));       
       _graphptr->setError(gsize,_graphptr->getError(max_error_index));
   }

   for(int i = 0; i < _graphptr->size(); i++)
   {
    _graphptr->setError(i, this->params[1] * _graphptr->getError(i));
   }
  }  
 }
}

} // namespace neuralgas

#endif
