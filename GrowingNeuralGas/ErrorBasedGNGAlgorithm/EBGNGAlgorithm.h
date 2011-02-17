/** 
* \class EBGNGAlgorithm
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

template<typename T,typename S> class EBGNGAlgorithm : public GNGModul<T,S>
{
public:
	// cto initializing the class 
	EBGNGAlgorithm(const int& dim);
	// std dto
	~EBGNGAlgorithm();
        
	// runs the proposed algorithm 
	void    run();
        
	// sets the number of inital reference vectors
	virtual void    setRefVectors(const int&,const int&);
	
	void    showGraph(){_graphptr->showGraph();}
	// algorithmic dependent distance function
	T       getDistance(const Vector<T>&,const int&) const;
private:        
	//is a Base_Graph casted pointer to thereof derived class MGNGGraph
	GNGGraph<T,S>*           _graphptr;
	//defines the update rule for the in the second index given neighbor by using a given datum
	void updateNeighbor(const int&,const int&);
	//defines the update rule for the winner
	void updateWinner(const int&,const int&);
	//a learning cycle for instance
        void learning_loop ( unsigned int );

	// the average error
	float average_error;
	// rate parameter for adjusting neighbour weights
	float rate;
        // errors queue for calculating \p average_error
        std::queue<T> errors; 
        
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S> EBGNGAlgorithm<T,S>::EBGNGAlgorithm(const int& dim): GNGModul<T,S>(dim)
{
 _graphptr=NULL;
 int num_of_elem=50;
 for(int i=0; i < num_of_elem; i++)
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


/** \brief Sets the number of initial reference vectors.
*
* The function sets the initial number of reference vectors and initializes those
* with random values that are within the range of the given parameter max_value.
* \param num_of_ref_vec is the number of initial reference vectors
* \param max_value is the max value that shall be used for the random init value generation
*/
template<typename T,typename S> void EBGNGAlgorithm<T,S>::setRefVectors(const int& num_of_ref_vec,const int& max_value)
{
  if (_graphptr!=NULL)
      delete _graphptr;
  if (this->graphptr!=NULL)
     delete this->graphptr;
     
  average_error       = 0.0;
  rate = 0.5;
  
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


template<typename T,typename S> T EBGNGAlgorithm<T,S>::getDistance(const Vector<T>& item, const int& node_index) const
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
template<typename T,typename S> void EBGNGAlgorithm<T,S>::updateNeighbor(const int& time,const int& node_index)
{
 //(*_graphptr)[node_index].weight  += this->params[5] * ( (*this)[time]-(*_graphptr)[node_index].weight);
 (*_graphptr)[node_index].weight  += rate * ( (*this)[time]-(*_graphptr)[node_index].weight);

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
template<typename T,typename S> void EBGNGAlgorithm<T,S>::updateWinner(const int& time,const int& winner)
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
    }
    else if (this->sampling_mode == randomly)
    {
      ::srand( (unsigned)time( NULL ) );
      if (this->stopping_criterion == epochs)
        for (unsigned int e=0; e<this->max_epochs; e++)
          for (unsigned int t=0; t<tsize; t++)
            learning_loop (::rand() % tsize);

    }
    std::cout << "Graph size: " << _graphptr->size()<< " " << std::endl; 
  }
  std::cout << "Average error: " << average_error << " " << std::endl;
}

template<typename T,typename S> void EBGNGAlgorithm<T,S>::learning_loop ( unsigned int t)
{
  std::cout << "t = " << t << std::endl;
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
      
  this->getWinner(first_winner,second_winner,(*this)[t]);

  //   T distance = pow(getDistance((*this)[t],first_winner),2);
  T distance = getDistance((*this)[t],first_winner);

  

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

    std::vector<int> first_winner_neighbors = _graphptr->getNeighbors(first_winner);
   
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


} // namespace neuralgas

#endif
