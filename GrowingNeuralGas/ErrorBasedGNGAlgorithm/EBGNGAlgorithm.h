/** 
* \file EBGNGAlgorithm.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef EBGNGALGORITHM_H
#define EBGNGALGORITHM_H

#include <GrowingNeuralGas/GNGModul.h>
#include <GrowingNeuralGas/GNGGraph.h>
#include <queue>

namespace neuralgas {

/** \class EBGNGAlgorithm
 *  \brief A simple implementation of GNG that modifies learning rate with the current
           average error. Stops according to error threshold and maximum network size
 */
template<typename T,typename S> class EBGNGAlgorithm : public GNGModul<T,S>
{
public:
	// cto initializing the class 
	EBGNGAlgorithm(const unsigned int& dim);
	// std dto
	~EBGNGAlgorithm();
        
	// runs the proposed algorithm 
	void    run();
        
	// sets the number of inital reference vectors
	virtual void    setRefVectors(const unsigned int&,const Vector<T>&, const Vector<T>&);
	// check graph stability as a stopping criterion
	bool isGraphStable ();
	/// show graph for visualization
	void    showGraph(){_graphptr->showGraph();}
	// set error threshold constant
	void setErrorThreshold (T);
private:        
	//is a Base_Graph casted pointer to thereof derived class MGNGGraph
	GNGGraph<T,S>*           _graphptr;
	//defines the update rule for the in the second index given neighbor by using a given datum
	void updateNeighbor(const unsigned int&,const unsigned int&);
	//defines the update rule for the winner
	void updateWinner(const unsigned int&,const unsigned int&);
	//a learning cycle for instance
        void learning_loop ( unsigned int );

	// the average error
	T average_error;
	// rate parameter for adjusting neighbour weights
	T rate;
	// error_threshold
	T error_threshold;
        // errors queue for calculating \p average_error
        std::queue<T> errors; 
        
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S> EBGNGAlgorithm<T,S>::EBGNGAlgorithm(const unsigned int& dim): GNGModul<T,S>(dim)
{
 _graphptr=NULL;
 unsigned int num_of_elem=50;
 for(unsigned int i=0; i < num_of_elem; i++)
   errors.push(this->_zero);
 
}

/** \brief std dto
*/
template<typename T,typename S> EBGNGAlgorithm<T,S>::~EBGNGAlgorithm()
{
//if (_graphptr!=NULL) 
	delete _graphptr;
   // _graphptr=NULL;                          
}

/** \brief Sets the number of initial reference vectors. This is the first function
 * to be called.
 *
 * The function sets the initial number of reference vectors and initializes those
 * with random values that are within the range of the given parameter max_value.
 * \param num_of_ref_vec is the number of initial reference vectors
 * \param low_limit is the min value that shall be used for the random init value generat * \param high_limit is the max value that shall be used for the random init value generation
*/
template<typename T,typename S>
void EBGNGAlgorithm<T,S>::setRefVectors(const unsigned int& num_of_ref_vec,const Vector<T>& low_limits, const Vector<T>& high_limits)
{
	if (_graphptr!=NULL)
		delete _graphptr;
	if (this->graphptr!=NULL)
		delete this->graphptr;
     
	average_error       = 0.0;
	rate = 0.5;

	_graphptr           = new GNGGraph<T,S>(this->getDimension());
	this->graphptr      = _graphptr;
	this->_graphModulptr = _graphptr;
	// sets the min value for the init of the context vector
	_graphptr->setLowLimits(low_limits);
	// sets the max value for the init of the context vector
	_graphptr->setHighLimits(high_limits);
	// creates a Graph object with given size of the 
	// vectors and number of ref vectors initilized with 
	// random values
	_graphptr->initRandomGraph(num_of_ref_vec, low_limits, high_limits);
}

/** \brief Defines the update rule for the neighbor given by the second index 
*
*   The update rule depends on the actual datum. With that datum a given
*   topological neighbor is updated by an algorithmic dependent rule.
*   w_j(new) = w_j(old) + epsilon_n * ( x_t - w_j(old) )
*
*   \param item_index is the data vector index that is used for updating 
*   \param node_index is the index of topological neighbor that shall be updated
*/
template<typename T,typename S> void EBGNGAlgorithm<T,S>::updateNeighbor(const unsigned int& item_index,const unsigned int& node_index)
{
 //(*_graphptr)[node_index].weight  += this->params[5] * ( (*this)[time]-(*_graphptr)[node_index].weight);
 (*_graphptr)[node_index].weight  += rate * ( (*this)[item_index]-(*_graphptr)[node_index].weight);

}

/** \brief defines the update rule for the winner
*
*   The update rule depends on the actual datum. With that datum the given
*   winner is updated by an algorithmic dependent rule.
*    w(new) = w(old) + epsilon_b * (x_t - w(old) )
*
*   \param item_index is the data vector index that is used for updating 
*   \param winner is the index of the winner that shall be updated
*/
template<typename T,typename S> void EBGNGAlgorithm<T,S>::updateWinner(const unsigned int& item_index,const unsigned int& winner)
{
 (*_graphptr)[winner].weight  += this->params[4] * ( (*this)[item_index]-(*_graphptr)[winner].weight);
}


/** \brief Runs the algorithm */
template<typename T,typename S> void EBGNGAlgorithm<T,S>::run()
{
 
  if (this->getDimension()>0)
  {
    // number of total steps for the algorithm
    unsigned int tsize                       =     this->size();
    /*   for( int j = 0; j < NUM_PARAM; j++)  
         {
         this->params[j] =((*this)._funcArray[j])(0);
     
         }*/
    this->params[3] =((*this)._funcArray[3])(0);
    this->params[4] =((*this)._funcArray[4])(0);
    this->params[7] =((*this)._funcArray[7])(0);

    if (this->sampling_mode == sequential)
    {
	    if (this->stopping_criterion == epochs)
		    for (unsigned int e=0; e<this->max_epochs; e++)
			    for(unsigned int t = 0; t < tsize; t++)
				    learning_loop (t);
	    else if (this->stopping_criterion == stability)
		    do
		    {
			    for (unsigned int t = 0; t < tsize; t++)
				    learning_loop (t);
		    } while (!isGraphStable());
    }
    else if (this->sampling_mode == randomly)
    {
	    ::srand( (unsigned)time( NULL ) );
	    if (this->stopping_criterion == epochs)
		    for (unsigned int e=0; e<this->max_epochs; e++)
			    for (unsigned int t=0; t<tsize; t++)
				    learning_loop (::rand() % tsize);
	    else if (this->stopping_criterion == stability)
		    do {
			    for(unsigned int t = 0; t < tsize; t++)
				    learning_loop (::rand() % tsize);
		    }
		    while (!isGraphStable());
      
    }
    std::cout << "Graph size: " << _graphptr->size()<< " " << std::endl; 
  }
  std::cout << "Average error: " << average_error << " " << std::endl;
}

template<typename T,typename S> void EBGNGAlgorithm<T,S>::learning_loop ( unsigned int t)
{
  // std::cout << "t = " << t << std::endl;
  // init
    
  unsigned int first_winner               =     1;
  unsigned int second_winner              =     0;
 
  //params are defined by user set functions that may depend on the time step
  //params[0] alpha factor by which the two nodes with greatest error are changed after adding a new nod
  //params[1] beta factor by which the error of all nodes is changed
  //params[2] gamma maximal permitted distortion before adding a new reference
  //params[3] alpha_max maximal age for an edge 
  //params[4] epsilon_b weight in winner node update
  //params[5] epsilon_n weight in neighbor node update
  //params[6] lambda number of iterations after a new node is added
  //params[7] w_max maximal number of reference vectors / maximal size of cookbook
      
  _graphptr->getWinner(first_winner,second_winner,(*this)[t]);

  //   T distance = pow(_graphptr->getDistance((*this)[t],first_winner),2);
  T distance = _graphptr->getDistance((*this)[t],first_winner);
  

  if ( _graphptr->size() < this->params[7] &&  distance >= average_error)
  
  {
    _graphptr->addNode();
    int gsize = _graphptr->size()-1;
    //set the position of the new node to be the position of the current data vector
    (*_graphptr)[gsize].weight = (*this)[t];
    _graphptr->setAge(first_winner,gsize,0.0);
    _graphptr->setAge(second_winner,gsize,0.0);
     
    //new node's error is set to the error of the winner node 
    // _graphptr->setError(gsize,_graphptr->getError(first_winner) );  
  } 
  else
  {
    //int gsize = _graphptr->size()-1;

    std::vector<unsigned int> first_winner_neighbors = _graphptr->getNeighbors(first_winner);
   
    for(unsigned int j=0; j < first_winner_neighbors.size();j++)
    {
	    _graphptr->incAge(first_winner, first_winner_neighbors[j] ); 
	    updateNeighbor(t,first_winner_neighbors[j]);
    }
     
    updateWinner(t,first_winner);      
    errors.push(distance);
    float elem=errors.front();
    float av_error_old = average_error;
    average_error = average_error + ( (distance - elem) / errors.size());  
    
    errors.pop();
    if (average_error !=0.0 )
      rate =((1.0  - av_error_old  / average_error)>0) ? 1.0  - av_error_old  / average_error : -(1.0  - av_error_old  / average_error) ;
    if (first_winner != second_winner)
	    _graphptr->setAge(first_winner,second_winner,0.0);
  }

  this->rmOldEdges(this->params[3]);  
    
  this->rmNotConnectedNodes(); 


}

/*  \brief check if graph is stable by using a low error criterion and a expected maximal nr of nodes
 */
template<typename T, typename S>
bool EBGNGAlgorithm<T,S>::isGraphStable ()
{
    std::cout << "Graph size: " << _graphptr->size()<< " " << std::endl; 
    std::cout << "Average error: " << average_error << " " << std::endl;
	if (_graphptr->size() == this->params[7] && error_threshold > average_error)
		return true;
	return false;
}

/** \brief set \p error_threshold
    \param threshold given error threshold */
template<typename T, typename S>
void EBGNGAlgorithm<T,S>::setErrorThreshold (T threshold)
{
	error_threshold = threshold;
}

	
} // namespace neuralgas

#endif
