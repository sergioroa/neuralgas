/*
 *   This file is part of NeuralGas.
 *
 *   NeuralGas is free software: you can redistribute it and/or modify it
 *   under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   NeuralGas is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with NeuralGas.  If not, see <http://www.gnu.org/licenses/>.
 */

/** 
* \file MGNGAlgorithm.h
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

// TODO   params replaced by telling variable names within the run func, getDistance

#ifndef MGNGALGORITHM_H
#define MGNGALGORITHM_H

#include "MGNGGraph.h"
#include <NeuralGas.h>
#include <GrowingNeuralGas/GNGModul.h>
#include <limits>
#include <fstream>
#include <iostream>

namespace neuralgas {

/**  \brief func is the default updating rule for a node's counter
*
*    \param n is the node to be updated
*    \param value is a float value that can be included, default it is params[7] = eta
*/
template<typename T,typename S> void updateCounter(Base_Node<T,S>* n,const float& value)
{
  MGNGNode<T,S>* temp = dynamic_cast< MGNGNode<T,S>* > (n);
  temp->setCounter(temp->getCounter() * value);
}


/** \class MGNGAlgorithm
 *  \brief Class implements the algorithm proposed in "Incremental Unsupervised Time Series 
 *   Analysis using Merge Growing Neural Gas" by Andreakis,Hoyningen-Huene and Beets.
 *
 * The class is a flexible (not yet as flexible as it is intented to be) implemention of
 * the algorithm proposed by the above mentionend authors.
 * Since the class is derived from the abstract class NeuralGas it is necessary to define
 * the precompiler macro NUM_PARAM reflecting the number of needed parameter params within the
 * algorithm.
 * Recap : the values are set via the super class func setFuncArray taking a func ptr of the
 *  float (*)(const unsigned int&) and an index of the parameter to set. The setting is
 * params[0] alpha, params[1] beta, params[2] gamma, params[3] delta, params[4] epsilon_b, 
 * params[5] epsilon_n, params[6] theta params[7] eta params[8] lambda
 * In order to change the number of parameters the precompile definition NUM_PARAM in file 
 * NeuralGas.h has to be changed accordingly.
 * The int parameter in float (*)(const unsigned int&) is the given time step allowing to define
 * time / step dependend parameters.
 * A further flexiblization is func is the setFuncUpdateCounter taking a func ptr of the form
 *  void (*)(Base_Node<T,S>* n, const float&)
 * defining the update rule for each node, the second float parameter is param[7].
 *
 * IMPORTANT : USAGE
 *
 * setFuncArray(index) has to be called for each parameters(#NUM_PARAM)
 *
 * setData 
 *
 * maximal random value has to be determined
 *
 * setRefVectors(num,max_random) inits graph with num ref vectors, their random init weight
 *                               and context values is within max_random.
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
 * setFuncUpdateCounter(void (*)(Base_Node<T,S>* n,const float&))  Defining the update rule for 
 *                                                       each node, the second float parameter 
 *                                                       is param[7].
 *                                                       Called with NULL sets default rule counter*param[7]
 *
 *
 *   \param globalContextV is the vector reflecting the global time series context
 *   \param _graphptr is a Base_Graph casted pointer to thereof derived class MGNGGraph
 *   \param funcUpdateCounter is func pointer to a func with the time as paramter defining the way of updating the counter of each node
 */
template<typename T,typename S> class MGNGAlgorithm : public GNGModul<T,S>
{
public:
	// cto initializing the class 
	MGNGAlgorithm(const unsigned int& dim);
	// std dto
	~MGNGAlgorithm();
                         
	// runs the proposed algorithm 
	void    run();
           
	// sets the number of inital reference vectors
	virtual void    setRefVectors(const unsigned int&,const Vector<T>&, const Vector<T>&);
	// sets the rule for updating the node counter
	void    setFuncUpdateCounter(void (*)(Base_Node<T,S>* n,const float&));
	void    showGraph(){_graphptr->showGraph();}
	// stores the graph in myfile , just for internal use
	void    storeGraph(const unsigned int& );
	friend class MGNGGraph<T,S>;
protected:
	virtual void updateNeighbor(const unsigned int&,const unsigned int&);
	virtual void updateWinner(const unsigned int&,const unsigned int&);
        
	//vector reflecting the global time series context
	Vector<T> globalContextV;
	//is a Base_Graph casted pointer to thereof derived class MGNGGraph
	MGNGGraph<T,S>*           _graphptr;
	//is func pointer to a func with the time as paramter defining 
	// the way of updating the counter of each node
	void (*_funcUpdateCounter)(Base_Node<T,S>* n,const float&);
   
                 
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S> MGNGAlgorithm<T,S>::MGNGAlgorithm(const unsigned int& dim): GNGModul<T,S>(dim)
{_graphptr=NULL;_funcUpdateCounter=updateCounter;}

/** \brief std dto
*/
template<typename T,typename S> MGNGAlgorithm<T,S>::~MGNGAlgorithm()
{
//if (_graphptr!=NULL) 
    delete _graphptr;
   // _graphptr=NULL;                          
}

/** \brief sets the rule for updating the node counter
*   
*   After a user defined update rule was set, a call of this function with a NULL as parameter
*   resets the default update rule.
* 
*   \param *funcUpdateCounter) a function ptr pointing to the func defining the update rule for the counter
*/

template<typename T,typename S> void MGNGAlgorithm<T,S>::setFuncUpdateCounter(void (*funcUpdateCounter)(Base_Node<T,S>* n,const float&)=NULL)
{
 if (funcUpdateCounter == NULL)
    _funcUpdateCounter = updateCounter;
 else
    _funcUpdateCounter = funcUpdateCounter;
}

/** \brief Sets the number of initial reference vectors.
*
* The function sets the initial number of reference vectors and initializes those
* with random values that are within the range of the given parameter max_value.
* \param num_of_ref_vec is the number of initial reference vectors
* \param max_value is the max value that shall be used for the random init value generation
*/
template<typename T,typename S> void MGNGAlgorithm<T,S>::setRefVectors(const unsigned int& num_of_ref_vec,const Vector<T>& low_limits, const Vector<T>& high_limits)
{
  if (_graphptr!=NULL)
      delete _graphptr;
  if (this->graphptr!=NULL)
     delete this->graphptr;
  
  _graphptr           = new MGNGGraph<T,S>(this->getDimension());
  // DANGER DownCast is performed via dynamic_cast
  //_graphptr       = dynamic_cast< MGNGGraph<T,S> * >(this->_graphModulptr);
  this->graphptr      = _graphptr;
  this->_graphModulptr = _graphptr;
  _graphptr->setAlgorithm (this);

  // sets the min values for the init of the context vector
  _graphptr->setLowLimits(low_limits);
  // sets the max values for the init of the context vector
  _graphptr->setHighLimits(high_limits);

  _graphptr->initRandomGraph(num_of_ref_vec,low_limits, high_limits); // creates a Graph object with given size of the 
                                                // vectors and number of ref vectors initilized with 
                                                // random values
  for (unsigned int i=0;i < num_of_ref_vec; i++)
      (*_graphptr).setBirthday(i,0);
}


/** \brief Defines the update rule for the in the second index given neighbor 
*
*   The update rule depends on the actual datum. With that datum a given
*   topological neighbor is updated by an algorithmic dependent rule.
*   w_j(new) = w_j(old) + epsilon_n * ( x_t - w_j(old) )
*   c_j(new) = c_j(old) + epsilon_n * (C - c_j(old) )
*
*   \param time is the data vector that is used for updating 
*   \param node_index is the index of topological neighbor that shall be updated
*/
template<typename T,typename S> void MGNGAlgorithm<T,S>::updateNeighbor(const unsigned int& time,const unsigned int& index)
{
 (*_graphptr)[index].weight  += this->params[5] * ( (*this)[time] - (*_graphptr)[index].weight);
 (*_graphptr).context(index) += this->params[5] * (globalContextV - (*_graphptr).context(index));
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
template<typename T,typename S> void MGNGAlgorithm<T,S>::updateWinner(const unsigned int& time,const unsigned int& winner)
{
 (*_graphptr)[winner].weight  += this->params[4] * ( (*this)[time]-(*_graphptr)[winner].weight);
 (*_graphptr).context(winner) += this->params[4] * ( globalContextV - (*_graphptr).context(winner));
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
template<typename T,typename S> void MGNGAlgorithm<T,S>::run()
{
 
 if (this->getDimension()>0)
 {
  // the line numbers refer to the lines given in the above mentioned paper
  // line 1 - 3
  // number of total steps for the algorithm
   unsigned int tsize                       =     this->size();
  // line 4
   globalContextV.resize(this->getDimension(),this->_zero);  
  //line 5
   for(unsigned int t = 0; t < tsize; t++)
   {
    // init
    unsigned int first_winner;
    unsigned int second_winner;
    

   //params are defined by user set functions that may depend on the time step
   //params[0] alpha used as weight for the distance metric in getDistance
   //params[1] beta used as weight for the global context computation
   //params[2] gamma maximal age for an edge
   //params[3] delta factor by which the counter of a new node is changed
   //params[4] epsilon_b weight in winner node update
   //params[5] epsilon_n weight in neighbor node update
   //params[6] theta maximal number of reference nodes / size of cookbook
   //params[7] eta factor by which the counter of all nodes is changed
   //params[8] lambda number of iterations after a new node is added
    for( int j = 0; j < NUM_PARAM; j++)  
    {
      this->params[j] =((*this)._funcArray[j])(t);
    }

    // line 6
    _graphptr->getWinner(first_winner,second_winner,(*this)[t]);
    // line 7
    Vector<T> vec = (*_graphptr).context(first_winner);
    // C_t = (1 - beta) * w_r +beta * c_r
    globalContextV=( 1 - this->params[1]) * (*_graphptr)[first_winner].weight;
    globalContextV+= this->params[1] * vec;
    // line 10
    _graphptr->incCounter(first_winner);
    // line 11
    updateWinner(t,first_winner);
        
    std::vector<unsigned int> first_winner_neighbors = _graphptr->getNeighbors(first_winner);
   
    for(unsigned int j=0; j < first_winner_neighbors.size();j++)
    {
      updateNeighbor(t,first_winner_neighbors[j]);
      //line 12
      _graphptr->incAge(first_winner, first_winner_neighbors[j] ); 
    //  (*_graphptr).setBirthday(first_winner_neighbors[j],t);
      
    }
    //(*_graphptr).setBirthday(first_winner,t);
    
   // line 9
   // line 9 is realised after 10-12 for efficiency reasons, in line 12 we had otherwise
   // to check in each for-loop with an if statement whether the neighbor is the second_winner
   // in order to not increment the edge age
   // therefore we possibly (if the second_winner was already a neighbor) increment the age first
   // and set it then to 0 by creating an edge ( if there was already an edge, nothing is changed 
   // by the command )
   
    // line 8 - 9
    if (first_winner != second_winner)
       _graphptr->setAge(first_winner,second_winner,0.0);
    // line 13
    this->rmOldEdges(this->params[2]);  
    // line 14
    this->rmNotConnectedNodes();    
    // line 15
   
    if(  ((t % (unsigned int)(this->params[8]))== 0) && _graphptr->size() < this->params[6])
    {
       float max_counter                   = 0.0;
       int max_counter_index               = 0;
       // line 15a
       for (unsigned int i=0; i < _graphptr->size(); i++)
       {
           if (_graphptr->getCounter(i) > max_counter )
           {
            max_counter             = _graphptr->getCounter(i);
            max_counter_index       = i;
           }
       }
 
       std::vector<unsigned int> neighbors   = _graphptr->getNeighbors(max_counter_index);
       float max_counter_n                 = 0.0;
       int max_counter_index_n             = 0;
       // line 15b
       for (unsigned int i=0; i < neighbors.size(); i++)
       {

           if (_graphptr->getCounter(neighbors[i]) > max_counter_n )
           {
            max_counter_n             = _graphptr->getCounter(neighbors[i]);
            max_counter_index_n       = i;
           }
       }
       // line 15c
       _graphptr->addNode();

       unsigned int gsize                       = _graphptr->size() - 1;
       (*_graphptr).setBirthday(gsize,t);
       (*_graphptr)[gsize].weight      = (*_graphptr)[max_counter_index_n].weight*0.5
                                         +(*_graphptr)[max_counter_index].weight*0.5;
        _graphptr->context(gsize)       = _graphptr->context(max_counter_index_n)*0.5
                                          +_graphptr->context(max_counter_index)*0.5;
       // line 15d
       _graphptr->rmEdge(max_counter_index, max_counter_index_n);  
       _graphptr->addEdge(max_counter_index,gsize);
       _graphptr->addEdge(max_counter_index_n,gsize);
       // line 15e
       _graphptr->context(max_counter_index)   = (1- this->params[3]) * _graphptr->context(max_counter_index);                        
       _graphptr->context(max_counter_index_n) = (1- this->params[3]) * _graphptr->context(max_counter_index_n);
   }
   // line 16
     applyFunc2AllNodes(_funcUpdateCounter,this->params[7]);
  }
 }
}

template<typename T,typename S> void MGNGAlgorithm<T,S>::storeGraph(const unsigned int& t)
{
  std::ofstream myfile("graph.txt",std::ios::out | std::ios::app);
  
 for(unsigned int i=0; i < _graphptr->size(); i++)
 {
	 myfile << "connex " <<(*_graphptr)[i].num_connections << " ";
	 for(unsigned int j=0; j < (*_graphptr)[i].edges.size(); j++)
		 myfile <<(*_graphptr)[i].edges[j] << " ";    
	 myfile << "counter "<<(*_graphptr).getCounter(i)<<" ";
	 myfile << "weight ";
	 for(unsigned int j=0; j < (*_graphptr)[i].weight.size(); j++)
		 myfile << (*_graphptr)[i].weight[j]<<" ";
	 myfile << "context ";
	 for(unsigned int j=0; j < _graphptr->context(i).size(); j++)
		 myfile << _graphptr->context(i)[j]<<" ";
	 
	 myfile << std::endl;
 }
 myfile << "data ";
 for (unsigned int i=0; i < (*this)[t].size(); i++)
     myfile << (*this)[t][i] <<" ";
 myfile<<std::endl;
 myfile << std::endl;
  myfile.close();
}

} // namespace neuralgas

#endif
