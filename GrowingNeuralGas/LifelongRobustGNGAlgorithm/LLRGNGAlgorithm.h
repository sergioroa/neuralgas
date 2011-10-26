/** 
* \file LLRGNGAlgorithm.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef LLRGNGALGORITHM_H
#define LLRGNGALGORITHM_H

#include <GrowingNeuralGas/GNGModul.h>
#include <GrowingNeuralGas/GNGGraph.h>
#include "LLRGNGGraph.h"
#include <algorithm>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <VoronoiDiagramm/VoronoiWidget.h>

namespace neuralgas {

//! \class LLRGNGThread
/*! \brief This class is used to provide threading facilities
 */
class LLRGNGThread : public QThread {
	Q_OBJECT
	
public:
	LLRGNGThread ()
	{
		debugging = false;
	}
	~LLRGNGThread ()
	{
		mutex.lock();
		condition.wakeOne();
		mutex.unlock();
	}
	virtual void run () = 0;
	void setDebugging (bool debug)
	{
		debugging = debug;
	}
	void begin ()
	{
		if (!isRunning ()) {
			start(LowPriority);
		}
		else {
			condition.wakeOne();
		}
		
	}
	QMutex* getMutex ()
	{
		return &mutex;
	}
	QWaitCondition* getWaitCondition ()
	{
		return &condition;
	}
signals:
	void updateData (SeqNeurons* neurons);
	void initializeData (SeqData* data, SeqNeurons* neurons, unsigned int sidesize = 1000);
protected:
	QMutex mutex;
	QWaitCondition condition;
	bool debugging;

};

/** \class LLRGNGAlgorithm
 *  \brief Class implements some techniques based on the algorithms explained in 
 *  Robust growing neural gas algorithm with application in cluster analysis
 *  by A.K. Qin and P.N. Suganthan and
 *  Life-long learning cell structures -- continuosly learning without catastrophic interference
 *         by Fred H. Hamker.
 *  In this new algorithm, there is no output layer. The learning process is unsupervised. There are stopping criteria based on Minimum Description Length and Learning progress.
*/
template<typename T, typename S>
class LLRGNGAlgorithm : public GNGModul<T,S>, public LLRGNGThread
{
public:

	// cto class initialization
	LLRGNGAlgorithm (const unsigned int& dim, const unsigned int&);
	// std dto
	~LLRGNGAlgorithm ();

	// run the algorithm
	void run ();

	//sets the number of initial reference vectors
	virtual void setRefVectors(const unsigned int&,const Vector<T>&, const Vector<T>&);
	// set time window constants
	void setTimeWindows (unsigned int, unsigned int, unsigned int);
	// set input adaptation threshold
	void setAdaptationThreshold (T);
	// set initial learning rate constants
	void setLearningRates (T, T);
	// set insertion rate
	void setInsertionRate (T);
	// set maximal edge age constant
	void setMaximalEdgeAge (unsigned int);
	// set data accuracy (quantization) constant
	void setDataAccuracy (T);
	// set model efficiency contribution constant
	void setModelEfficiencyConst (T);
	// find node with maximal value of insertion criterion
	int maxInsertionCriterionNode ();
	// find two nodes with maximal value of insertion criterion
	std::vector<int> maxInsertionCriterionNodes ();
	// find node with maximal value of insertion quality among some nodes
	int maxInsertionQualityNode (const std::vector<unsigned int>&);
	/// show graph for visualization
	void showGraph(){_graphptr->showGraph();}
	// sets a maximal partition
	void setMaxNodes (unsigned int nr) { max_nodes = nr;}
	// gets maximal partition
	unsigned int getMaxNodes () const;
	// calculate the Minimum Description Length of the current graph and dataset
	T calculateMinimumDescriptionLength ();
	// set maximum nr of epochs after avg error reduction is expected
	void setMaxEpochsErrorReduction (unsigned int);
	// set maximum nr of epochs after mdl reduction is expected
	void setMaxEpochsMDLReduction (unsigned int);
protected:
	//check if minimal error for all nodes has not changed for more than \p max_epochs_improvement
	bool minimalAvgErrors ();
	// update the graph with minimal description length
	void updateMinimalGraphMDL ();
	// check if minimal MDL has not changed for more than \p max_epochs_mdl_reduction
	bool minimalMDL ();
	// find graph with a node that is not properly located according to MDL principle
	LLRGNGGraph<T,S>* findDislocatedNodeGraph (unsigned int&, unsigned int&, unsigned int&);
	// calculate model efficiency (error) of a graph
	T calculateModelEfficiency (LLRGNGGraph<T,S>* graph);
	// mark the current graph as stable in order to stop the algorithm
	void markAsStableGraph ();
	// calculate value range of data
	void calculateValueRange ();
	// calculate initial restrictive distances for every node
	void calculateInitialRestrictingDistances ();
	//defines the update rule for the winner node by using a given data vector index
	void updateWinnerWeight(const unsigned int&, const unsigned int&, T& distance);
	//defines the update rule for a neighboring node by using a given data vector index
	void updateNeighborWeight(const unsigned int&, const unsigned int&, const unsigned int&, std::map<unsigned int, T>&);
	/// Base_Graph casted pointer to thereof derived class GNGGraph
	LLRGNGGraph<T,S>* _graphptr;
	//a learning cycle for instance
        void learning_loop ( unsigned int, unsigned int );
	// update params like age, avg errors, etc.
	void updateParams (unsigned int&, unsigned int&, unsigned int&);
	/// insertion rate constant
	T insertion_rate;
	/// maximal number of nodes if it is set (not obligatory)
	unsigned int max_nodes;
	/// accuracy constant
	T data_accuracy;
	/// current training epoch
	unsigned int epoch;
	/// current minimum description length value
	T mdl;
	/// current minimal mdl
	T min_mdl;
	/// current minimal graph (mdl criterion)
	LLRGNGGraph<T,S>* min_mdl_graphptr;
	/// maximal nr of epochs for checking for mean errors non-reduction
	unsigned int max_epochs_error_reduction;
	/// maximal nr of epochs for checking \p mdl reduction
	unsigned int max_epochs_mdl_reduction;
	/// last epoch when mdl reduction found
	unsigned int last_epoch_mdl_reduction;
	/// used as stopping criterion (graph stability)
	bool stable_graph;
	/// value range of data
	T value_range;
	/// constant to balance contribution of model efficiency in minimum description length calculation
	T model_efficiency_const;
	/// number of bits needed to encode a vector (calculated in \p calculateValueRange)
	T bits_vector;
	unsigned int nr_of_insertions;
	/// flag for tracking if a dislocated node was found in previous insertion epoch
	bool insert_this_iter;
};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S>
LLRGNGAlgorithm<T,S>::LLRGNGAlgorithm(const unsigned int& dim, const unsigned int& window = 200):
	GNGModul<T,S>(dim),
	max_nodes (0),
	data_accuracy (0.001),
	mdl (1e10),
	min_mdl (mdl),
	min_mdl_graphptr (NULL),
	max_epochs_error_reduction (5),
	max_epochs_mdl_reduction (80),
	last_epoch_mdl_reduction (0),
	model_efficiency_const (1.0)
{
	_graphptr = new LLRGNGGraph<T,S>(dim, window);
	this->graphptr = _graphptr;
	this->_graphModulptr = _graphptr;
}

/** \brief std dto
*/
template<typename T,typename S>
LLRGNGAlgorithm<T,S>::~LLRGNGAlgorithm()
{
	delete _graphptr;
	if (this->stopping_criterion == stability && min_mdl_graphptr != NULL)
		delete min_mdl_graphptr;
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
void LLRGNGAlgorithm<T,S>::setRefVectors(const unsigned int& num_of_ref_vec,const Vector<T>& low_limits, const Vector<T>& high_limits)
{
	assert (num_of_ref_vec >= 2);
	assert (_graphptr->size() == 0);
	// sets the min value for the init of the context vector
	_graphptr->setLowLimits(low_limits);
	// sets the max value for the init of the context vector
	_graphptr->setHighLimits(high_limits);
	// creates a Graph object with given size of the 
	// vectors and number of ref vectors initilized with 
	// random values
	_graphptr->initRandomGraph(num_of_ref_vec, low_limits, high_limits);

	calculateInitialRestrictingDistances ();

	if (debugging)
		emit initializeData (this->_data, _graphptr->getNodes());

}

/** \brief set window constants. This function should only
    be called once when creating the graph
    \param smoothing smoothing time window constant
    \param longterm error time window constant
    \param age age time window constant */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setTimeWindows (unsigned int smoothing, unsigned int error, unsigned int age)
{
	_graphptr->setTimeWindows (smoothing, error, age);
}


/** \brief set initial learning rate constants
 *  \param winner_learning_rate initial winner learning rate constant
 *  \param neighbors_learning_rate initial winner-neighbors learning rate constant
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setLearningRates (T winner_learning_rate, T neighbors_learning_rate)
{
	_graphptr->setLearningRates (winner_learning_rate, neighbors_learning_rate);
}

/** \brief set adaptation threshold
    \param adaptation_threshold given adaptation threshold */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setAdaptationThreshold (T adaptation_threshold)
{
	_graphptr->setAdaptationThreshold (adaptation_threshold );	
}

/** \brief set insertion rate
    \param rate insertion rate (used as factor of nr. of nodes) */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setInsertionRate (T rate)
{
	insertion_rate = rate;	
}

/** \brief set maximal edge age constant
 *  \param maximal_edge_age maximal edge age constant
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setMaximalEdgeAge (unsigned int maximal_edge_age)
{
	_graphptr->setMaximalEdgeAge (maximal_edge_age) ;
}

/** set \p data_accuracy constant
 *  \param accuracy maximal data precision or minimal quantization constant
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setDataAccuracy (T accuracy)
{
	data_accuracy = accuracy;
}

/** set \p model_efficiency_const constant to balance the contribution of model efficiency for MDL calculation
 * \param constant set model complexity contribution
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setModelEfficiencyConst (T constant)
{
	model_efficiency_const = constant;
}

/** \brief Defines the update rule for a node given by the second index 
*   See Qin and Suganthan (2004) update rule.
*
*   \param item_index is the data vector index that is used for updating 
*   \param node_index is the index of topological neighbor that shall be updated
*   \param distances_winner distances of winner neighbors (used only for neighbors updating as a repulsive force
*/
template<typename T,typename S>
void LLRGNGAlgorithm<T,S>::updateWinnerWeight(const unsigned int& item_index,const unsigned int& node_index, T& distance)
{
	LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* >(&(*_graphptr)[node_index]);
	// T distance = getDistance ((*this)[item_index], node_index);
	T amplitude;
	if (distance >= node->prev_restricting_distance)
		amplitude = node->restricting_distance;
	else
		amplitude = distance;
	
	node->weight += node->learning_rate * amplitude * ( (*this)[item_index]-node->weight) / distance;

}

/** \brief Defines the update rule for a node given by the second index 
*
*   See Qin and Suganthan (2004) update rule, taking into account
*   repulsive forces.
*
*   \param item_index is the data vector index that is used for updating 
*   \param node_index is the index of topological neighbor that shall be updated
*   \param winner_index is the index of the winner node
*   \param distances_winner distances of winner neighbors (used for repulsive force)
*/
template<typename T,typename S>
void LLRGNGAlgorithm<T,S>::updateNeighborWeight(const unsigned int& item_index,const unsigned int& node_index, const unsigned int& winner_index, std::map<unsigned int, T>& distances_winner)
{
	LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* >(&(*_graphptr)[node_index]);
	// LLRGNGNode<T,S>* winner = static_cast<LLRGNGNode<T,S>* >(&(*_graphptr)[winner_index]);

	T distance = _graphptr->getDistance ((*this)[item_index], node_index);
	T amplitude;
	if (distance >= node->prev_restricting_distance)
		amplitude = node->restricting_distance;
	else
		amplitude = distance;
	
	node->weight += node->learning_rate * amplitude * ( (*this)[item_index]-node->weight) / distance;

	T dist_avg;
	typename std::map<unsigned int, T>::iterator it;
	for (it = distances_winner.begin(); it != distances_winner.end(); it++)
		dist_avg += it->second;
	
	dist_avg /= distances_winner.size();

	node->weight += exp (-distances_winner[node_index] / node->repulsion) * 2 * dist_avg * (node->weight - (*_graphptr)[winner_index].weight) / distances_winner[node_index];

}

/** \brief Runs an algorithm based on different implementations of GNG and own ideas
*/
template<typename T,typename S>
void LLRGNGAlgorithm<T,S>::run()
{
	assert (this->getDimension() > 0);
	//setRefVectors() must be called before
	assert (_graphptr);
	assert (insertion_rate);

	unsigned int tsize = this->size();
	if (this->sampling_mode == sequential)
	{
		if (this->stopping_criterion == epochs)
			for (epoch = 0; epoch < this->max_epochs; epoch++)
				for(unsigned int t = 0; t < tsize; t++)
					learning_loop (t, t);
		else if (this->stopping_criterion == stability)
		{
			epoch = 0;
			stable_graph = false;
			nr_of_insertions = 1;
			insert_this_iter = true;
			calculateValueRange ();
			do
			{
				for(unsigned int t = 0; t < tsize; t++)
				{
					learning_loop (t, t);
					if (stable_graph)
						break;
				}
				if (stable_graph)
					break;
				epoch++;
			}
			while (true);
			std::cout << "Finished!" << std::endl;
		}
	}
	else if (this->sampling_mode == randomly)
	{
		::srand( (unsigned)time( NULL ) );
		if (this->stopping_criterion == epochs)
			for (epoch = 0; epoch < this->max_epochs; epoch++)
				for(unsigned int t = 0; t < tsize; t++)
					learning_loop (::rand() % tsize, t);
		else if (this->stopping_criterion == stability)
		{
			epoch = 0;
			stable_graph = false;
			nr_of_insertions = 1;
			insert_this_iter = true;
			calculateValueRange ();
			do {
				for(unsigned int t = 0; t < tsize; t++)
				{
					learning_loop (::rand() % tsize, t);
					if (stable_graph)
						break;
				}
				if (stable_graph)
					break;
				epoch++;
				
			}
			while (true);
			std::cout << "Finished!" << std::endl;
		}

	}
}

/** \brief a learning cycle for instance t
 *  \param t index in the dataset
 *  \param i iteration counter
 */
template<typename T,typename S>
void LLRGNGAlgorithm<T,S>::learning_loop ( unsigned int t, unsigned int i )
{
	//initialization of first winner
	unsigned int b;
	//second winner
	unsigned int s;
	T distance = _graphptr->getWinner(b,s,(*this)[t]);
	// _graphptr->increaseItemsCounter (b);
	
	//learning rule for weight adaptation
	//after calculating learning quality for best matching node b and neighbors
	//based on input learning rate, input adaptation threshold and nodes age
	_graphptr->calculateLearningQuality(b);
	_graphptr->updateWinnerLearningRate (b);

	updateWinnerWeight(t,b,distance);
	std::vector<unsigned int> b_neighbors = _graphptr->getNeighbors(b);
	std::map<unsigned int, T> b_neighbors_distances = _graphptr->get1toMDistances (b, b_neighbors);
	
	
	for(unsigned int j=0; j < b_neighbors.size();j++)
	{
		unsigned int b_neighbor = b_neighbors[j];
		_graphptr->calculateLearningQuality(b_neighbor);
		_graphptr->updateNeighborLearningRate (b_neighbor);
		
		updateNeighborWeight(t,b_neighbor,b, b_neighbors_distances);
	}

	//calculate insertion quality
	if (i % (unsigned int)(insertion_rate /** _graphptr->size()*/) == 0 )
	{
		if (debugging)
		{
			mutex.lock();
			emit updateData(_graphptr->getNodes());
			condition.wait (&mutex);
			mutex.unlock();
		}
		if (min_mdl_graphptr != NULL)
			updateMinimalGraphMDL ();
		else
		{
			std::cout << "mdl: " << calculateMinimumDescriptionLength () << std::endl;
			min_mdl = mdl;
			min_mdl_graphptr = new LLRGNGGraph<T,S>(*_graphptr);
				
		}
		if (minimalMDL ())
		{
			markAsStableGraph ();
			return;
		}
		//find node with maximal value of insertion criterion when the graph is not improving more
		if (minimalAvgErrors ())
		{
			unsigned int dislocated_node;
			if (_graphptr->size() > 2)
			{
				if (_graphptr->deleteInactiveNodes(b, s))
				{
					updateMinimalGraphMDL();

					if (debugging)
					{
						mutex.lock();
						emit updateData(_graphptr->getNodes());
						condition.wait (&mutex);
						mutex.unlock();
					}

					if (minimalMDL ())
					{
						markAsStableGraph ();
						return;
					}

				}
				if (nr_of_insertions == 1)
				{
					LLRGNGGraph<T,S>* dislocated_node_graphptr = findDislocatedNodeGraph (b, s, dislocated_node);
					if (dislocated_node_graphptr != NULL)
					{
						std::cout << "deleting nodes..." << std::endl;
						delete _graphptr;
						_graphptr = dislocated_node_graphptr;
						this->graphptr = _graphptr;
						this->_graphModulptr = _graphptr;
						insert_this_iter = false;
						updateMinimalGraphMDL();
						if (debugging)
						{
							mutex.lock();
							emit updateData(_graphptr->getNodes());
							condition.wait (&mutex);
							mutex.unlock();
						}

						if (minimalMDL ())
						{
							markAsStableGraph ();
							return;
						}

					}
					else delete dislocated_node_graphptr;
						
				}
			}
			// std::vector<unsigned int> dislocated_nodes = findDislocatedNodes ();
			// std::cout << "nr. of dislocated nodes: " << dislocated_nodes.size() << std::endl;

			// unsigned int nr_of_insertions = 1;
			// if (dislocated_node_graphptr != NULL)
			// 	nr_of_insertions = 2;
			if (insert_this_iter == false)
			{
				insert_this_iter = true;
				updateParams (t,b,s);
				nr_of_insertions = 2;
				return;
			}
			else
			{
				for (unsigned int j=0; j<_graphptr->size(); j++)
				{
					_graphptr->calculateInsertionQuality (j);
					_graphptr->calculateInsertionCriterion (j);
				}
				
				std::vector<int> q_nodes = maxInsertionCriterionNodes();
				//for (unsigned int d=0; d<dislocated_nodes.size(); d++)
				//{
				for (unsigned int n=0; n<nr_of_insertions; n++)
				{
					// int q = maxInsertionCriterionNode ();
					int q = q_nodes[n];
		
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
							std::cout << "adding node " << node_index << "..." << std::endl;
							_graphptr->calculateInheritedParams (node_index, q, f);
							_graphptr->setAge(q,node_index,0.0);
							_graphptr->setAge(f,node_index,0.0);

							calculateInitialRestrictingDistances ();
							if (debugging)
							{
								mutex.lock();
								emit updateData(_graphptr->getNodes());
								condition.wait (&mutex);
								mutex.unlock();
							}

						}
					}
				}

				nr_of_insertions = 1;
			}
		}
		
	}

	if (!stable_graph)
	{
		updateParams (t,b,s);
	}


	//remove all edges older than the maximal value for age
	this->rmOldEdges (_graphptr->getMaximalEdgeAge());

	//remove nodes without any edge
	this->rmNotConnectedNodes();


}

//! \brief update params like age, avg errors
/*! 
  \param t current time step
  \param b current winner
  \param s current second winner
*/template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::updateParams (unsigned int& t, unsigned int& b, unsigned int& s)
{
	if (b < _graphptr->size())
	{
		T min_error = _graphptr->getNodeMinLastAvgError (b);
		T distance = _graphptr->getDistance ((*this)[t], b);
		_graphptr->updateAvgError (b, distance);
		_graphptr->updateRestrictingDistance (b, distance);
		std::vector<unsigned int> b_neighbors = _graphptr->getNeighbors(b);

		//for (unsigned int j=0; j < b_neighbors.size(); j++)
		//_graphptr->updateRestrictingDistance (j, getDistance((*this)[t], j));
		if (min_error > _graphptr->getNodeMinLastAvgError (b))
			_graphptr->setLastEpochImprovement (b, epoch);
		
		//decrease age of best-matching node
		_graphptr->decreaseNodeAge (b);

		//the LLG (Hamker) algorithm modifies the insertion threshold if the distribution of the error changes. Here we are assuming a stationary distribution

		//adapt edges
		if (s < _graphptr->size() && s != b)
		{
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
		}
	}

}

/** \brief find node with maximal value of insertion criterion
 */
template<typename T, typename S>
int LLRGNGAlgorithm<T,S>::maxInsertionCriterionNode ()
{
	assert (_graphptr->size());
	T max_value = 0;
	int q = -1;
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		T value = (static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]))->insertion_criterion;
		if (max_value < value)
		{
			max_value = value;
			q = i;
		}
	}
	return q;
}

/** \brief find two nodes with maximal values of insertion criterion
 */
template<typename T, typename S>
std::vector<int> LLRGNGAlgorithm<T,S>::maxInsertionCriterionNodes ()
{
	assert (_graphptr->size());
	T max_value = (static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[0]))->insertion_criterion;
	std::vector<int> q (2);
	q[0] = 0;
	for (unsigned int i=1; i < _graphptr->size(); i++)
	{
		T value = (static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]))->insertion_criterion;
		if (max_value < value)
		{
			max_value = value;
			q[1] = q[0];
			q[0] = i;
		}
	}
	return q;
}

 
/** \brief find node with maximal value of insertion quality among some nodes
 *  \param nodes vector of nodes to query
 */
template<typename T, typename S>
int LLRGNGAlgorithm<T,S>::maxInsertionQualityNode (const std::vector<unsigned int>& nodes)
{
	if (nodes.size() == 0)
		return -1;
	T max_value = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[nodes[0]])->insertion_quality;
	unsigned int f = nodes[0];
	for (unsigned int i=1; i<nodes.size(); i++)
	{
		T value = (static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[nodes[i]]))->insertion_quality;
		if (max_value < value)
		{
			max_value = value;
			f = nodes[i];
		}
	}
	return f;
}

//! \brief calculate model efficiency of a graph.
/*! \param graph 
  \return model efficiency
*/
template<typename T, typename S>
T LLRGNGAlgorithm<T,S>::calculateModelEfficiency (LLRGNGGraph<T,S>* graph)
{
	graph->resetMDLCounters ();
	for(unsigned int t = 0; t < this->size(); t++)
	{
		//winners
		unsigned int b, s;
		/*T distance = */graph->getWinner(b,s,(*this)[t]);
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*graph)[b]);
		node->items_counter++;
		// node->utifactor += distance;
		
		// winner_per_item[t] = b;
		T node_item_efficiency = 0;
		for (unsigned int i=0; i<this->getDimension(); i++)
		{
			node_item_efficiency += std::max(log2((fabs((*this)[t][i] -  (*graph)[b].weight[i])) / data_accuracy), 1.0);
			node->efficiency += node_item_efficiency;
		}
		graph->model_efficiency += node_item_efficiency;
	}
	/*for (unsigned int i=0; i < graph->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*graph)[i]);
		if (node->items_counter > 0)
			node->utifactor /= (T) node->items_counter;
		graph->utifactor += node->utifactor;
	}*/
	return graph->model_efficiency;


}

//! \brief Calculate MDL for the current \p _graphptr
template<typename T, typename S>
T LLRGNGAlgorithm<T,S>::calculateMinimumDescriptionLength ()
{
	// items_per_winner = std::vector<unsigned int> (_graphptr->size());
	//std::vector<int> winner_per_item (this->size(), -1);
	// std::vector<T> mdl_per_item (this->size());
	// mdl = 0;
	calculateModelEfficiency (_graphptr);
	mdl = model_efficiency_const * _graphptr->model_efficiency;
	T L_inliers = this->size() * log2 (_graphptr->size());
	
	mdl += _graphptr->size() * bits_vector + L_inliers;
	std::cout << "K: " << bits_vector << std::endl;
	std::cout << "L_inliers: " << L_inliers << std::endl;
	std::cout << "efficiency: " << mdl - (_graphptr->size() * bits_vector + L_inliers ) << std::endl;
	return mdl;

}

//! \brief find nodes that are not properly located according to MDL principle

/*! \return graph_neg_mdl_change graph with more negative mdl change containing a node that is dislocated
*/
template<typename T, typename S>
LLRGNGGraph<T,S>* LLRGNGAlgorithm<T,S>::findDislocatedNodeGraph (unsigned int& winner, unsigned int& snd_winner, unsigned int& dislocated_node)
{
	// test dislocated nodes creating pruned graphs
	std::vector<LLRGNGGraph<T,S>* > pruned_graphs;
	if (_graphptr->size() <= 2)
		return NULL;
	for (unsigned int i=0; i<_graphptr->size(); i++)
	{
		LLRGNGGraph<T,S>* pruned_graph = new LLRGNGGraph<T,S>(*_graphptr);
		pruned_graphs.push_back (pruned_graph);
		pruned_graph->rmNode (i);
		calculateModelEfficiency (pruned_graph);
	}
	dislocated_node = _graphptr->size();
	T min_change = 0;
	T model_complexity_change = -bits_vector + this->size() * (log2 (_graphptr->size() - 1) - log2 (_graphptr->size()));
	std::cout << "mod complex. change: " << model_complexity_change << std::endl;
	
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		T change = model_complexity_change + model_efficiency_const * (pruned_graphs[i]->model_efficiency - _graphptr->model_efficiency);
		if (change < min_change)
		{
			
			std::cout << "eff pruned: " << pruned_graphs[i]->model_efficiency << std::endl;
			std::cout << "change: " << change << std::endl;
			dislocated_node = i;
		}
	}
	LLRGNGGraph<T,S>* graph_neg_mdl_change;
	if (dislocated_node != _graphptr->size())
	{
		graph_neg_mdl_change = pruned_graphs[dislocated_node];
		if (winner > dislocated_node)
			winner--;
		else if (winner == dislocated_node)
			winner = _graphptr->size();
		if (snd_winner > dislocated_node)
			snd_winner--;
		else if (snd_winner == dislocated_node)
			snd_winner = _graphptr->size();
	}
	else
		graph_neg_mdl_change = NULL;
	for (unsigned int i=0; i < _graphptr->size(); i++)
		if (i != dislocated_node)
		{
			LLRGNGGraph<T,S>* graph_to_del = pruned_graphs[i];
			delete graph_to_del;
		}
	pruned_graphs.clear ();
	return graph_neg_mdl_change;
}



/** \brief get maximal number of partitions
 */
template<typename T, typename S>
unsigned int LLRGNGAlgorithm<T,S>::getMaxNodes () const
{
	if (max_nodes == 0)
		return _graphptr->size ();
	else
		return max_nodes;
}

/** \brief set maximum nr of epochs after avg error reduction is expected
    \param epochs maximum nr of epochs
*/
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setMaxEpochsErrorReduction (unsigned int epochs)
{
	max_epochs_error_reduction = epochs;
}

/** \brief set maximum nr of epochs after mdl reduction is expected
    \param epochs maximum nr of epochs
*/
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setMaxEpochsMDLReduction (unsigned int epochs)
{
	max_epochs_mdl_reduction = epochs;
}


/** \brief check if minimal error for all nodes has not changed for more than \p max_epochs_error_reduction
 */
template<typename T, typename S>
bool LLRGNGAlgorithm<T,S>::minimalAvgErrors ()
{
 	for (unsigned int i=1; i < _graphptr->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]);
		if (epoch - node->last_epoch_improvement < max_epochs_error_reduction)
			return false;
		
	}
	return true;
}


/** \brief update the graph with minimal description length
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::updateMinimalGraphMDL ()
{

	std::cout << "Nr. of nodes now: " <<  _graphptr->size() << std::endl;
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]);
		// std::cout << "node " << i << " age: " << node->age << std::endl;
		std::cout << "node " << i << " learning q.: " << node->learning_quality << std::endl;
		std::cout << "node " << i << " last avg e.: " << node->last_avgerror << std::endl;
		std::cout << "node " << i << " prev avg e.: " << node->prev_avgerror << std::endl;
		std::cout << "node " << i << " repulsion: " << node->repulsion << std::endl;
		std::cout << "node " << i << " restricting d.: " << node->restricting_distance << std::endl;
		// std::cout << "node " << i << " insertion c.: " << node->insertion_criterion << std::endl;
		std::cout << "node " << i << " learning rate: " << node->learning_rate << std::endl;
		// std::cout << "node " << i << " insertion q.: " << node->insertion_quality << std::endl;
	}

	std::cout << "mdl: " << calculateMinimumDescriptionLength ();

	if (min_mdl > mdl )
	{
		std::cout << "\t\t\t\t <-- is minimal!" << std::endl;
		min_mdl = mdl;
		delete min_mdl_graphptr;
		min_mdl_graphptr = new LLRGNGGraph<T,S>(*_graphptr);
		last_epoch_mdl_reduction = epoch;
	}
	else
		std::cout << std::endl;
	

}

/** \brief check if minimal MDL has not changed for more than \p max_epochs_mdl_reduction
 */
template<typename T, typename S>
bool LLRGNGAlgorithm<T,S>::minimalMDL ()
{
	return (epoch - last_epoch_mdl_reduction >= max_epochs_mdl_reduction);

}


/** \brief mark the current graph as stable in order to stop the algorithm
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::markAsStableGraph ()
{
	delete _graphptr;
	_graphptr = new LLRGNGGraph<T,S>(*min_mdl_graphptr);
	this->graphptr = _graphptr;
	this->_graphModulptr = _graphptr;
	stable_graph = true;

	if (debugging)
	{
		mutex.lock();
		emit updateData(_graphptr->getNodes());
		condition.wait (&mutex);
		mutex.unlock();
	}

}

/** \brief calculate value range of data
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::calculateValueRange ()
{
	T avg_values;
	for(unsigned int t = 0; t < this->size(); t++)
		for (unsigned int i=0; i<this->getDimension(); i++)
			avg_values += (*this)[t][i];
	avg_values /= T (this->size());
	value_range = avg_values - _graphptr->getLowLimit ();

	std::cout << "value range: " << value_range << std::endl;

	bits_vector = ceil (log2(value_range / data_accuracy));

}

/** \brief calculate initial restrictive distances for every node
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::calculateInitialRestrictingDistances ()
{
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]);
		node->restricting_distance = 0;

		for (unsigned int t=0; t<this->size(); t++)
		{
			node->restricting_distance += 1.0 / metric (node->weight, (*this)[t]);
		}
		node->restricting_distance /= this->size();
		node->restricting_distance = 1.0 / node->restricting_distance;
		node->prev_restricting_distance = node->restricting_distance;
	}
	
}


} // namespace neuralgas

#endif
