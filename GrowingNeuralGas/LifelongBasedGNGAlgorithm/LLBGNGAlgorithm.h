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
#include "LLBGNGGraph.h"

namespace neuralgas {

/** \brief Class implements some techniques proposed in the algorithm Life-long learning
 *         cell structures -- continuosly learning without catastrophic interference
 *         by Fred H. Hamker.
 *  Parameters like insertion and deletion criteria are added, whose purpose is a trade-off
 *  for dealing with bias-variance issues. Also Activation function?
*/
template<typename T, typename S>
class LLBGNGAlgorithm : public GNGModul<T,S>
{
public:

	// cto class initialization
	LLBGNGAlgorithm (const unsigned int& dim);
	// std dto
	~LLBGNGAlgorithm ();

	// run the algorithm
	void run ();
	//sets the number of inital reference vectors
	virtual void setRefVectors(const unsigned int&,const T&, const T&);
	// algorithmic dependent distance function
	T getDistance(const Vector<T>&,const unsigned int&) const;
	// set time window constants
	void setTimeWindows (unsigned int, unsigned int, unsigned int);
	// set input adaptation threshold
	void setAdaptationThreshold (T);
	// set initial learning rate constants
	void setLearningRates (T, T, T);
	// set insertion rate
	void setInsertionRate (T);
	// set insertion tolerance
	void setInsertionTolerance (T);
	// set deletion threshold
	void setDeletionThreshold (T);
	// set minimal node age constant
	void setMinimalNodeAge (T);
	// set maximal edge age constant
	void setMaximalEdgeAge (unsigned int);
	// set stabilization constant
	void setStabilization (T);
	//find node with maximal value of insertion criterion
	int maxInsertionCriterionNode ();
	//find node with maximal value of insertion quality among some nodes
	int maxInsertionQualityNode (const std::vector<unsigned int>&);
	//find node with minimal value of deletion criterion
	int minDeletionCriterionNode ();
	// for visualization
	void showGraph(){_graphptr->showGraph();}
protected:
	// calculate and return local similarity of weights for some node
	T calculateWeightsLocalSimilarity (const int);
	// return average similarity of weights
	T getWeightsAvgSimilarity ();
	//defines the update rule for the given node by using a given data vector index
	void updateWeight(const unsigned int&,const unsigned int&);
private:
	/// Base_Graph casted pointer to thereof derived class GNGGraph
	LLBGNGGraph<T,S>*           _graphptr;
	//a learning cycle for instance
        void learning_loop ( unsigned int, unsigned int );
	/// insertion rate constant
	T insertion_rate;
	
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S>
LLBGNGAlgorithm<T,S>::LLBGNGAlgorithm(const unsigned int& dim): GNGModul<T,S>(dim)
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

/** \brief Sets the number of initial reference vectors. This is the first function
 * to be called.
 *
 * The function sets the initial number of reference vectors and initializes those
 * with random values that are within the range of the given parameter max_value.
 * \param num_of_ref_vec is the number of initial reference vectors
 * \param low_limit is the min value that shall be used for the random init value generat * \param high_limit is the max value that shall be used for the random init value generation
*/
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::setRefVectors(const unsigned int& num_of_ref_vec,const T& low_limit, const T& high_limit)
{
	if (_graphptr!=NULL)
		delete _graphptr;
	if (this->graphptr!=NULL)
		delete this->graphptr;
     
	_graphptr           = new LLBGNGGraph<T,S>(this->getDimension());
	this->graphptr      = _graphptr;
	this->_graphModulptr = _graphptr;
	// sets the min value for the init of the context vector
	_graphptr->setLowLimit(low_limit);
	// sets the max value for the init of the context vector
	_graphptr->setHighLimit(high_limit);
	// creates a Graph object with given size of the 
	// vectors and number of ref vectors initilized with 
	// random values
	_graphptr->initRandomGraph(num_of_ref_vec, low_limit, high_limit);
}

/** \brief set window constants. This function should only
    be called once when creating the graph
    \param shortterm_window short term time window constant
    \param longterm_window long term time window constant
    \param age_time_window age time window constant */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setTimeWindows (unsigned int shortterm_window, unsigned int longterm_window, unsigned int age_time_window)
{
	assert (_graphptr != NULL);
	_graphptr->setTimeWindows (shortterm_window, longterm_window, age_time_window);
}

/** \brief set initial learning rate constants
 *  \param winner_learning_rate initial winner learning rate constant
 *  \param neighbors_learning_rate initial winner-neighbors learning rate constant
 *  \param insertion_learning_rate insertion threshold learning rate constant
 */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setLearningRates (T winner_learning_rate, T neighbors_learning_rate, T insertion_learning_rate)
{
	assert (_graphptr != NULL);
	_graphptr->setLearningRates (winner_learning_rate, neighbors_learning_rate, insertion_learning_rate);
}

/** \brief set adaptation threshold
    \param adaptation_threshold given adaptation threshold */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setAdaptationThreshold (T adaptation_threshold)
{
	assert (_graphptr != NULL);
	_graphptr->setAdaptationThreshold (adaptation_threshold );	
}

/** \brief set insertion tolerance
    \param insertion_tolerance insertion tolerances */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setInsertionTolerance (T insertion_tolerance)
{
	assert (_graphptr != NULL);
	_graphptr->setInsertionTolerance (insertion_tolerance);	
}

/** \brief set insertion rate
    \param rate insertion rate (used as factor of nr. of nodes) */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setInsertionRate (T rate)
{
	assert (_graphptr != NULL);
	insertion_rate = rate;	
}

/** \brief set \p deletion_threshold
    \param deletion_threshold given deletion threshold */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setDeletionThreshold (T deletion_threshold)
{
	assert (_graphptr != NULL);
	_graphptr->setDeletionThreshold (deletion_threshold);	
}

/* \brief set \p minimal_node_age constant
 * \param minimal_node_age minimal node age constant
 */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setMinimalNodeAge (T minimal_node_age)
{
	assert (_graphptr != NULL);
	_graphptr->setMinimalNodeAge (minimal_node_age);
}

/* \brief set maximal edge age constant
 * \param maximal_edge_age maximal edge age constant
 */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setMaximalEdgeAge (unsigned int maximal_edge_age)
{
	assert (_graphptr != NULL);
	_graphptr->setMaximalEdgeAge (maximal_edge_age) ;
}

/* set \p stabilization constant
 * \param stabilization learning stabilization constant
 */
template<typename T, typename S>
void LLBGNGAlgorithm<T,S>::setStabilization (T stabilization)
{
	assert (_graphptr != NULL);
	_graphptr->setStabilization( stabilization );
}

/** \brief Algorithmic dependent distance function
*
*   This function returns the distance of the given item and the given node. 
*   The distance is a algorithmic dependent function that is either
*   just the setted metric or a combination thereof.
*   Currently dist  = metric(x_t,w_j) where x_t is the data vector and w_j the node vector,
*
*   \param item data vector
*   \param node_index is the node where to the distance shall be determined
*/
template<typename T,typename S>
T LLBGNGAlgorithm<T,S>::getDistance(const Vector<T>& item, const unsigned int& node_index) const
{
    // dist  = metric(x_t,w_j) instead of metric(x_t,w_j)^2 as proposed in the paper
    // since this accelerates the calculation but does not change the result
	return metric( item, (*_graphptr)[node_index].weight);
}

/** \brief Defines the update rule for a node given by the second index 
*
*   The update rule depends on the current item. With that item a given
*   topological neighbor is updated by an algorithmic dependent rule.
*   w_j(new) = w_j(old) + epsilon_n * ( x_t - w_j(old) )
*
*   \param item_index is the data vector index that is used for updating 
*   \param node_index is the index of topological neighbor that shall be updated
*/
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::updateWeight(const unsigned int& item_index,const unsigned int& node_index)
{
	LLBGNGNode<T,S>* node = static_cast<LLBGNGNode<T,S>* >(&(*_graphptr)[node_index]);
	
	node->weight  += node->learning_rate * ( (*this)[item_index]-node->weight);

}

/** \brief Runs an algorithm based on Life-Long Learning Cell Structures algorithm proposed by Hamker
*/
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::run()
{
	assert (this->getDimension() > 0);
	//setRefVectors() must be called before
	assert (_graphptr);
	assert (insertion_rate);
	
	unsigned int tsize = this->size();
	if (this->sampling_mode == sequential)
	{
		if (this->stopping_criterion == epochs)
			for (unsigned int e=0; e<this->max_epochs; e++)
				for(unsigned int t = 0; t < tsize; t++)
					learning_loop (t, t);
	}
	else if (this->sampling_mode == randomly)
	{
		::srand( (unsigned)time( NULL ) );
		if (this->stopping_criterion == epochs)
			for (unsigned int e=0; e<this->max_epochs; e++)
				for(unsigned int t = 0; t < tsize; t++)
					learning_loop (::rand() % tsize, t);
	}
}

/** \brief a learning cycle for instance t
 *  \param t index in the dataset
 *  \param i iteration counter
 */
template<typename T,typename S>
void LLBGNGAlgorithm<T,S>::learning_loop ( unsigned int t, unsigned int i )
{
	//initialization of first winner
	unsigned int b = 1;
	//second winner
	unsigned int s = 0;
	this->getWinner(b,s,(*this)[t]);

	//learning rule for weight adaptation
	//after calculating learning quality for best matching node b and neighbors
	//based on input learning rate, input adaptation threshold and nodes age
	_graphptr->calculateLearningQuality(b);
	_graphptr->updateWinnerLearningRate (b);

	updateWeight(t,b);
	std::vector<unsigned int> b_neighbors = _graphptr->getNeighbors(b);
	
	for(unsigned int j=0; j < b_neighbors.size();j++)
	{
		unsigned int b_neighbor = b_neighbors[j];
		_graphptr->calculateLearningQuality(b_neighbor);
		_graphptr->updateNeighborLearningRate (b_neighbor);
		
		updateWeight(t,b_neighbor);
	}


	//calculate insertion quality
	if (i % (unsigned int)(insertion_rate * _graphptr->size()))
	{
		for (unsigned int j=0; j<_graphptr->size(); j++)
		{
			_graphptr->calculateInsertionQuality (j);
			_graphptr->calculateInsertionCriterion (j);
		}
		
		//find node with maximal value of insertion criterion
		int q = maxInsertionCriterionNode ();
		
		// find node among neighbours of q with maximal value of
		//quality measure of insertion
		if (q != -1)
		{
			std::vector<unsigned int> q_neighbors = _graphptr->getNeighbors(q);
			int f = maxInsertionQualityNode (q_neighbors);
			if (f != -1)
			{
				
				_graphptr->rmEdge (q, f);
				_graphptr->addNode ();
				int node_index = _graphptr->size()-1;
				_graphptr->calculateInheritedParams (node_index, q, f);
				_graphptr->setAge(q,node_index,0.0);
				_graphptr->setAge(f,node_index,0.0);
			}
		}
		//find node with minimal value of deletion criterion
		int d = minDeletionCriterionNode ();
		
		if (d != -1)
			_graphptr->rmNode (d);
			
			
	}
	
	//update long-term and short-term error variables
	_graphptr->updateAvgError (b, getDistance ((*this)[t], b));

	//decrease age of best-matching node
	_graphptr->decreaseNodeAge (b);

	//the original algorithm modifies the insertion threshold if the distribution of the error changes. Here we are assuming a stationary distribution

	//adapt edges
	bool s_neighborof_b = false;
	for(unsigned int j=0; j < b_neighbors.size();j++)
	{
		unsigned int b_neighbor = b_neighbors[j];
		if (b_neighbor == s)
		{
			s_neighborof_b = true;
			_graphptr->setAge (b, s, 0.0);
		}
		else
			_graphptr->incAge (b, b_neighbor);
	}
	if (!s_neighborof_b)
		_graphptr->setAge (b, s, 0.0);

	//remove all edges older than the maximal value for age
	this->rmOldEdges (_graphptr->getMaximalEdgeAge());

	//remove nodes without any edge
	this->rmNotConnectedNodes(); 

}

/** \brief find node with maximal value of insertion criterion
 */
template<typename T, typename S>
int LLBGNGAlgorithm<T,S>::maxInsertionCriterionNode ()
{
	assert (_graphptr->size());
	T max_value = 0;
	unsigned int q = -1;
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		T value = (static_cast<LLBGNGNode<T,S>* > (&(*_graphptr)[i]))->insertion_criterion;
		if (max_value < value)
		{
			max_value = value;
			q = i;
		}
	}
	return q;
}

/** \brief find node with maximal value of insertion quality among some nodes
 *  \param nodes vector of nodes to query
 */
template<typename T, typename S>
int LLBGNGAlgorithm<T,S>::maxInsertionQualityNode (const std::vector<unsigned int>& nodes)
{
	if (nodes.size() == 0)
		return -1;
	T max_value = (static_cast<LLBGNGNode<T,S>* > (&(*_graphptr)[nodes[0]]))->insertion_quality;
	unsigned int f = nodes[0];
	for (unsigned int i=1; i<nodes.size(); i++)
	{
		T value = (static_cast<LLBGNGNode<T,S>* > (&(*_graphptr)[nodes[i]]))->insertion_quality;
		if (max_value < value)
		{
			max_value = value;
			f = nodes[i];
		}
	}
	return f;
}

/* \brief calculate and return local similarity of weights for some node
 * \param index node index
 */
template<typename T, typename S>
T LLBGNGAlgorithm<T,S>::calculateWeightsLocalSimilarity (const int index)
{
	LLBGNGNode<T,S>* node = static_cast<LLBGNGNode<T,S>* >(&(*_graphptr)[index]);
	std::vector<unsigned int> neighbors = _graphptr->getNeighbors(index);
	node->local_similarity = this->_zero;
	for(unsigned int i=0; i < neighbors.size();i++)
		node->local_similarity += metric (node->weight, (static_cast<LLBGNGNode<T,S>* > (&(*_graphptr)[neighbors[i]]))->weight);

	node->local_similarity /= (T) neighbors.size();
	return node->local_similarity;
	
}

/* \brief return average similarity of weights. Uses \p calculateWeightsLocalSimilarity
   to store the local weight similarities for all nodes, which will be used by
   \p minDeletionCriterionNode
 */
template<typename T, typename S>
T LLBGNGAlgorithm<T,S>::getWeightsAvgSimilarity ()
{
	T similarity = this->_zero;
	for(unsigned int i=0; i < _graphptr->size();i++)
		similarity += calculateWeightsLocalSimilarity (i);

	similarity /= (T) _graphptr->size();
	return similarity;

}

/* \brief find node with minimal value of deletion criterion
 */
template<typename T, typename S>
int LLBGNGAlgorithm<T,S>::minDeletionCriterionNode ()
{
	assert (_graphptr->size());
	T avg_similarity = getWeightsAvgSimilarity ();
	T min_value = _graphptr->getDeletionThreshold ();
	int d = -1;
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		LLBGNGNode<T,S>* node = static_cast<LLBGNGNode<T,S>* > (&(*_graphptr)[i]);
		if (_graphptr->getNeighborsSize(i) >= 2 && node->age < _graphptr->getMinimalNodeAge() && node->learning_quality < _graphptr->getStabilization())
		{
			T deletion_criterion = node->local_similarity / avg_similarity;
			if (min_value > deletion_criterion)
			{
				min_value = deletion_criterion;
				d = i;
			}
		}
	}
	return d;
		
}



} // namespace neuralgas

#endif
