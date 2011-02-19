/** 
* \class LLBGNGAlgorithm
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef LLBGNGALGORITHM_H
#define LLBGNGALGORITHM_H

#include <GrowingNeuralGas/GNGModul.h>
#include <GrowingNeuralGas/GNGGraph.h>

namespace neuralgas {

/** \brief Class implements some techniques proposed in the algorithm Life-long learning
 *         cell structures -- continuosly learning without catastrophic interference
 *         by Fred H. Hamker.
 *  Parameters like insertion and deletion criteria are added, whose purpose is a trade-off
 *  for dealing with bias-variance issues. Also Activation function?
*/
class LLBGNGAlgorithm : public GNGModul<T,S>
{
public:

	// cto class initialization
	LLBGNGAlgorithm (const int& dim);
	// std dto
	~LLBGNGAlgorithm ();

	// run the algorithmv
	void run ();
	//sets the number of inital reference vectors
	virtual void setRefVectors(const int&,const int&);
	// algorithmic dependent distance function
	T       getDistance(const Vector<T>&,const int&) const;
private:
	/// Base_Graph casted pointer to thereof derived class GNGGraph
	GNGGraph<T,S>*           _graphptr;
	//defines the update rule for the in the second index given neighbor by using a given datum
	void updateNeighbor(const int&,const int&);
	//defines the update rule for the winner
	void updateWinner(const int&,const int&);
	//a learning cycle for instance
        void learning_loop ( unsigned int );
	
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S>
LLBGNGAlgorithm<T,S>::LLBGNGAlgorithm(const int& dim): GNGModul<T,S>(dim)
{
	_graphptr=NULL;
}

/** \brief std dto
*/
template<typename T,typename S>
LLBGNGAlgorithm<T,S>::~LLBGNGAlgorithm()
{
	delete _graphptr;
}

/** \brief Sets the number of initial reference vectors.
*
* The function sets the initial number of reference vectors and initializes those
* with random values that are within the range of the given parameter max_value.
* \param num_of_ref_vec is the number of initial reference vectors
* \param max_value is the max value that shall be used for the random init value generation
*/
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::setRefVectors(const int& num_of_ref_vec,const int& max_value)
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
template<typename T,typename S>
T LLBGNGAlgorithm<T,S>::getDistance(const Vector<T>& item, const int& node_index) const
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
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::updateNeighbor(const int& time,const int& node_index)
{
 //(*_graphptr)[node_index].weight  += this->params[5] * ( (*this)[time]-(*_graphptr)[node_index].weight);
	(*_graphptr)[node_index].weight  += (*_graphptr)[node_index].learning_rate * ( (*this)[time]-(*_graphptr)[node_index].weight);

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
template<typename T,typename S> void LLBGNGAlgorithm<T,S>::updateWinner(const int& time,const int& winner)
{
	(*_graphptr)[winner].weight  += (*_graphptr)[winner].learning_rate * ( (*this)[time]-(*_graphptr)[winner].weight);
}

/** \brief Runs an algorithm based on Life-Long Learning Cell Structures algorithm proposed by Hamker
*/
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::run()
{
	assert (this->getDimension() > 0);
	int tsize = this->size();
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
				for(unsigned int t = 0; t < tsize; t++)
					learning_loop (::rand() % tsize);
	}
}

/** \brief a learning cycle for instance t
 *  \param t index in the dataset
 */
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::learning_loop ( unsigned int t )
{
	//initialization
	int first_winner               =     1;
	int second_winner              =     0;
	this->getWinner(first_winner,second_winner,(*this)[t]);

	

	

}

} // namespace neuralgas

#endif
