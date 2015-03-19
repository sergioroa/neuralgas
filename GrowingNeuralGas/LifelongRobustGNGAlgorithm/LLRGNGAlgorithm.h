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
		visualizing = false;
	}
	~LLRGNGThread ()
	{
		mutex.lock();
		condition.wakeOne();
		mutex.unlock();
	}
	virtual void run () = 0;
	void setVisualizing (bool visualize)
	{
		visualizing = visualize;
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
	bool visualizing;

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
	friend class boost::serialization::access;
public:

	// default cto.
	LLRGNGAlgorithm ();
	// cto class initialization
	LLRGNGAlgorithm (const unsigned int&, const unsigned int&);
	// copy cto
	LLRGNGAlgorithm (const LLRGNGAlgorithm&);
	// std dto
	~LLRGNGAlgorithm ();

	// run the algorithm
	void run ();

	// sets the number of initial reference vectors
	virtual void setRefVectors(const unsigned int&);
	// set time window constants
	void setTimeWindows (unsigned int, unsigned int, unsigned int);
	// reset max errors size
	void setMaxErrorsSize (unsigned int);
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
	// find node with maximal value of insertion quality among some nodes
	int maxInsertionQualityNode (const std::vector<unsigned int>&);
	/// show graph for visualization
	void showGraph(){_graphptr->showGraph();}
	// sets a maximal partition
	void setMaxNodes (unsigned int nr) { max_nodes = nr;}
	// gets maximal partition
	unsigned int getMaxNodes () const;
	// calculate the Minimum Description Length of the current graph and dataset
	T calculateMinimumDescriptionLength (bool calculate_model_efficiency = true);
	// set maximum nr of epochs after avg error reduction is expected
	void setMaxEpochsErrorReduction (unsigned int);
	// set maximum nr of epochs after mdl reduction is expected
	void setMaxEpochsMDLReduction (unsigned int);
	// set mean distance calculation mode
	void setMeanDistanceMode (unsigned int);
	// save MDL history in text files
	void saveMDLHistory (std::string, bool append = false);
	// resets some learning related variables when restarting a experiment
	void resetLearning ();
	// recursively find nodes that are not properly located according to MDL principle
	void findDislocatedNodesStableGraph ();
	/// get average error among all graph nodes
	T getAvgError ();
	// get model efficiency (error) of current graph
	T getModelEfficiency ();
	// calculate value range of data
	void calculateValueRange ();
protected:
	//check if minimal error for all nodes has not changed for more than \p max_epochs_improvement
	bool minimalAvgErrors ();
	// update the graph with minimal description length
	void updateMinimalGraphMDL (bool calculate_model_efficiency = true);
	// check if minimal MDL has not changed for more than \p max_epochs_mdl_reduction
	bool minimalMDL ();
	// calculate model efficiency (error) of a graph
	T calculateModelEfficiency (LLRGNGGraph<T,S>* graph);
	T calculateModelEfficiency (LLRGNGGraph<T,S>* graph, unsigned int&);
	// mark the current graph as stable in order to stop the algorithm
	void markAsStableGraph ();
	// calculate initial restricting distances for every node
	void calculateInitialRestrictingDistances ();
	// check if restricting distances are overflowed
	void checkOverflowedRestrictingDistances ();
	// save MDL history in text files
	void saveMDLHistory ();
	// close MDL history file
	void closeMDLHistory ();
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
	// UGraph<T,S>* min_mdl_graphptr;
	/// temporal graph to store a graph without a dislocated node of \p _graphptr that minimizes \p mdl
	LLRGNGGraph<T,S>* dislocated_node_graphptr;
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
	/// mode for calculating mean distances
	unsigned int mean_distance_mode;
	/// number of bits needed to encode a vector (calculated in \p calculateValueRange)
	T bits_vector;
	/// file for saving MDL history
	std::ofstream* mdl_history;
private:
	template<class Archive>
	void serialize(Archive & ar, const unsigned int);

};

/** \brief cto initializing the class 
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S>
LLRGNGAlgorithm<T,S>::LLRGNGAlgorithm(const unsigned int& dim, const unsigned int& window = 81):
	GNGModul<T,S>(dim),
	max_nodes (0),
	data_accuracy (0.001),
	mdl (1e10),
	min_mdl (mdl),
	min_mdl_graphptr (NULL),
	max_epochs_error_reduction (5),
	max_epochs_mdl_reduction (80),
	last_epoch_mdl_reduction (0),
	model_efficiency_const (1.0),
	mean_distance_mode (harmonic),
	mdl_history (NULL)
{
	_graphptr = new LLRGNGGraph<T,S>(dim, window);
	this->graphptr = _graphptr;
	this->_graphModulptr = _graphptr;
}

/** \brief Default constructor (Should only be used for serialization)
*
* \param dim is the dimension of the node weights
*/
template<typename T,typename S>
LLRGNGAlgorithm<T,S>::LLRGNGAlgorithm():
	GNGModul<T,S>(0),
	max_nodes (0),
	data_accuracy (0),
	mdl (0),
	min_mdl (0),
	min_mdl_graphptr (0),
	max_epochs_error_reduction (0),
	max_epochs_mdl_reduction (0),
	last_epoch_mdl_reduction (0),
	model_efficiency_const (0),
	mean_distance_mode (harmonic),
	mdl_history (0)
{
}

/** \brief Copy constructor
 */
template<typename T,typename S>
LLRGNGAlgorithm<T,S>::LLRGNGAlgorithm(const LLRGNGAlgorithm& l) :
	GNGModul<T,S>(l),
	max_nodes (l.max_nodes),
	data_accuracy (l.data_accuracy),
	mdl (l.mdl),
	min_mdl (l.min_mdl),
	min_mdl_graphptr (0),
	max_epochs_error_reduction (l.max_epochs_error_reduction),
	max_epochs_mdl_reduction (l.max_epochs_mdl_reduction),
	last_epoch_mdl_reduction (l.last_epoch_mdl_reduction),
	model_efficiency_const (l.model_efficiency_const),
	mean_distance_mode (l.mean_distance_mode),
	mdl_history (0)
{
	_graphptr = new LLRGNGGraph<T,S>(*l._graphptr);
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
void LLRGNGAlgorithm<T,S>::setRefVectors(const unsigned int& num_of_ref_vec)
{
	assert (num_of_ref_vec >= 2);
	assert (_graphptr->size() == 0);
	assert (this->size());
	// sets the min value for the init of the context vector
	_graphptr->setLowLimits(this->minValues());
	// sets the max value for the init of the context vector
	_graphptr->setHighLimits(this->maxValues());
	// creates a Graph object with given size of the 
	// vectors and number of ref vectors initilized with 
	// random values
	_graphptr->initRandomGraph(num_of_ref_vec);

	if (mean_distance_mode == harmonic)
		calculateInitialRestrictingDistances ();

	if (visualizing)
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

/** \brief reset max errors size
    \param size
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setMaxErrorsSize (unsigned int size = 81)
{
	_graphptr->setMaxErrorsSize (size);
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

/** \brief set \p data_accuracy constant
 *  \param accuracy maximal data precision or minimal quantization constant
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setDataAccuracy (T accuracy)
{
	data_accuracy = accuracy;
}

/** \brief set \p model_efficiency_const constant to balance the contribution of model efficiency for MDL calculation
 *  \param constant set model complexity contribution
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setModelEfficiencyConst (T constant)
{
	model_efficiency_const = constant;
}

/** \brief set \p mdl_history file
 *  \param filename file name for saving the MDL history
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::saveMDLHistory (std::string filename, bool append)
{
	if (append)
		mdl_history = new std::ofstream (filename.c_str(), std::ios::out | std::ios::app );
	else
		mdl_history = new std::ofstream (filename.c_str());

	if (!mdl_history->is_open())
	{
		closeMDLHistory ();
	}
}

/** \brief store new value in \p mdl_history file
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::saveMDLHistory ()
{
	*mdl_history << mdl;
	if (epoch == last_epoch_mdl_reduction)
		*mdl_history << " " << min_mdl;
	*mdl_history << std::endl;
}

/** \brief close \p mdl_history file
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::closeMDLHistory ()
{
	mdl_history->close();
	delete mdl_history;
	mdl_history = NULL;
}

/** \brief reset some learning related variables when restarting a experiment
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::resetLearning ()
{
	// data have to be added before
	// sets the min value for the init of the context vector
	_graphptr->setLowLimits(this->minValues());
	// sets the max value for the init of the context vector
	_graphptr->setHighLimits(this->maxValues());
	if (mean_distance_mode == harmonic)
		calculateInitialRestrictingDistances ();

	this->epoch = 0;
	this->stable_graph = false;
	if (min_mdl_graphptr != NULL)
		delete min_mdl_graphptr;
	min_mdl_graphptr = NULL;
	this->last_epoch_mdl_reduction = 0;
	for (unsigned int i=0; i<this->_graphptr->size(); i++)
	{
		this->_graphptr->setLastEpochImprovement (i, 0);
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*this->_graphptr)[i]);
		node->min_last_avgerror = metric (_graphptr->getHighLimits (), _graphptr->getLowLimits());
		node->last_avgerror = node->prev_avgerror = -1;
		node->errors.clear();
		node->dim_errors.clear();
		for (unsigned int j=0; j<node->edges.size(); j++)
			if ( node->edges[j] != NULL )
				_graphptr->setAge (i, j, 0.0);
	}

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
	if (distance >= node->prev_restricting_distance && mean_distance_mode == harmonic)
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
	if (distance >= node->prev_restricting_distance && mean_distance_mode == harmonic)
		amplitude = node->restricting_distance;
	else
		amplitude = distance;
	
	node->weight += node->learning_rate * amplitude * ( (*this)[item_index]-node->weight) / distance;

	// T dist_avg;
	// typename std::map<unsigned int, T>::iterator it;
	// for (it = distances_winner.begin(); it != distances_winner.end(); it++)
	// 	dist_avg += it->second;
		
	// dist_avg /= distances_winner.size();
		
	// node->weight += exp (-distances_winner[node_index] / node->repulsion) * 2 * dist_avg * (node->weight - (*_graphptr)[winner_index].weight) / distances_winner[node_index];

	// node->weight += (node->weight - exp (-distances_winner[node_index] / winner->repulsion) * winner->weight);
	// node->weight += exp (-distances_winner[node_index] / winner->repulsion) * 2 * dist_avg * (node->weight - winner->weight) / distances_winner[node_index];


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
	calculateValueRange ();
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
			do
			{
				for(unsigned int t = 0; t < tsize; t++)
				{
					learning_loop (t, t);
					if (stable_graph)
						break;
				}
				if (mdl_history != NULL)
					saveMDLHistory ();
				if (stable_graph)
					break;
				epoch++;
			}
			while (true);
			*(this->out) << "Finished!" << std::endl;
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
			do {
				for(unsigned int t = 0; t < tsize; t++)
				{
					learning_loop (::rand() % tsize, t);
					if (stable_graph)
						break;
				}
				if (mdl_history != NULL)
					saveMDLHistory ();
				if (stable_graph)
					break;
				epoch++;
				
			}
			while (true);
			*(this->out) << "Finished!" << std::endl;
		}

	}
	if (mdl_history != NULL)
		closeMDLHistory ();
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
		// calculateInitialRestrictingDistances ();
		if (visualizing)
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
			*(this->out) << "mdl: " << calculateMinimumDescriptionLength () << std::endl;
			min_mdl = mdl;
			min_mdl_graphptr = new LLRGNGGraph<T,S>(*_graphptr);
			// min_mdl_graphptr = new UGraph<T,S>(dynamic_cast<UGraph<T,S>& >(*_graphptr));
				
		}
		if (minimalMDL ())
		{
			markAsStableGraph ();
			return;
		}
		//find node with maximal value of insertion criterion when the graph is not improving more
		if (minimalAvgErrors ())
		{
			if (_graphptr->size() > 2)
			{
				if (_graphptr->deleteInactiveNodes(b, s))
				{
					if (mean_distance_mode == harmonic)
						calculateInitialRestrictingDistances ();
					updateMinimalGraphMDL();

					if (visualizing)
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
			}

			for (unsigned int j=0; j<_graphptr->size(); j++)
			{
				_graphptr->calculateInsertionQuality (j);
				_graphptr->calculateInsertionCriterion (j);
			}
				
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
					*(this->out) << "adding node " << node_index << "..." << std::endl;
					_graphptr->calculateInheritedParams (node_index, q, f);
					_graphptr->setLastEpochImprovement (node_index, epoch);
					_graphptr->setAge(q,node_index,0.0);
					_graphptr->setAge(f,node_index,0.0);

							
					if (mean_distance_mode == harmonic)
						calculateInitialRestrictingDistances ();
					if (visualizing)
					{
						mutex.lock();
						emit updateData(_graphptr->getNodes());
						condition.wait (&mutex);
						mutex.unlock();
					}
				}
			}
		}		
	}

	if (!stable_graph)
	{
		updateParams (t,b,s);
	}


	//remove all edges older than the maximal value for age
	this->rmOldEdges (_graphptr->getMaximalEdgeAge());

	// //remove nodes without any edge
	// if (this->rmNotConnectedNodes() && mean_distance_mode == harmonic)
	// 	calculateInitialRestrictingDistances ();

	// checkOverflowedRestrictingDistances ();

}

//! \brief update params like age, avg errors
/*! 
  \param t current time step
  \param b current winner
  \param s current second winner
*/
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::updateParams (unsigned int& t, unsigned int& b, unsigned int& s)
{
	if (b < _graphptr->size())
	{
		Vector<T> dim_distances (this->getDimension());
		for (unsigned int i=0; i<this->getDimension(); i++)
			dim_distances[i] = T(fabs ((*this)[t][i] - (*_graphptr)[b].weight[i]));
		T min_error = _graphptr->getNodeMinLastAvgError (b);
		T distance = _graphptr->getDistance ((*this)[t], b);
		_graphptr->updateAvgError (b, distance, dim_distances);
		if (mean_distance_mode == harmonic)
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
		node->data.push_back (t);
		T node_item_efficiency = 0;
		for (unsigned int i=0; i<this->getDimension(); i++)
		{
			node_item_efficiency += std::max(log2((fabs((*this)[t][i] -  (*graph)[b].weight[i])) / data_accuracy), 1.0);
		}
		node->efficiency += node_item_efficiency;
		graph->model_efficiency += node_item_efficiency;
	}
	return graph->model_efficiency;


}

//! \brief recalculate model efficiency of a graph when looking for dislocated nodes
/*! \param graph
    \param rmnode_index index of the removed node
  \return model efficiency
*/
template<typename T, typename S>
T LLRGNGAlgorithm<T,S>::calculateModelEfficiency (LLRGNGGraph<T,S>* graph, unsigned int& rmnode_index)
{
	LLRGNGNode<T,S>* rmnode = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[rmnode_index]);
	for(unsigned int t = 0; t < rmnode->data.size(); t++)
	{
		unsigned int b, s;
		graph->getWinner(b,s,(*this)[rmnode->data[t]]);
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*graph)[b]);
		node->items_counter++;
		node->data.push_back (rmnode->data[t]);
		T node_item_efficiency = 0;
		for (unsigned int i=0; i<this->getDimension(); i++)
		{
			node_item_efficiency += std::max(log2((fabs((*this)[rmnode->data[t]][i] -  (*graph)[b].weight[i])) / data_accuracy), 1.0);
		}
		node->efficiency += node_item_efficiency;	
	}
	graph->model_efficiency = 0;
	for (unsigned int i=0; i < graph->size(); i++)
		graph->model_efficiency += static_cast<LLRGNGNode<T,S>& > ((*graph)[i]).efficiency;
	return graph->model_efficiency;
}

//! \brief Get model efficiency of the current \p _graphptr
template<typename T, typename S>
T LLRGNGAlgorithm<T,S>::getModelEfficiency ()
{
	return _graphptr->model_efficiency;
}

//! \brief Calculate MDL for the current \p _graphptr
template<typename T, typename S>
T LLRGNGAlgorithm<T,S>::calculateMinimumDescriptionLength (bool calculate_model_efficiency)
{
	// items_per_winner = std::vector<unsigned int> (_graphptr->size());
	//std::vector<int> winner_per_item (this->size(), -1);
	// std::vector<T> mdl_per_item (this->size());
	// mdl = 0;
	if (calculate_model_efficiency)
		calculateModelEfficiency (_graphptr);
	mdl = model_efficiency_const * _graphptr->model_efficiency;
	T L_inliers = this->size() * log2 (_graphptr->size());
	
	mdl += _graphptr->size() * bits_vector + L_inliers;
	*(this->out) << "K: " << bits_vector << std::endl;
	*(this->out) << "L_inliers: " << L_inliers << std::endl;
	*(this->out) << "efficiency: " << mdl - (_graphptr->size() * bits_vector + L_inliers ) << std::endl;
	return mdl;

}

/** \brief recursively find nodes that are not properly located according to MDL principle
    at the end of learning phase.
    Sets \p dislocated_node_graphptr and _graphptr to
    a graph without a dislocated node which minimizes MDL.
*/
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::findDislocatedNodesStableGraph ()
{
	// test dislocated nodes creating pruned graphs
	if (_graphptr->size() <= 2)
	{
		this->graphptr = _graphptr;
		this->_graphModulptr = _graphptr;
		return;
	}

	T min_change = 0;
	T model_complexity_change = -bits_vector + this->size() * (log2 (_graphptr->size() - 1) - log2 (_graphptr->size()));
	*(this->out) << "mod complex. change: " << model_complexity_change << std::endl;
	unsigned int dislocated_node = _graphptr->size();
	dislocated_node_graphptr = NULL;

	LLRGNGGraph<T,S>* pruned_graph; 
	for (unsigned int i=0; i<_graphptr->size(); i++)
	{
		pruned_graph = new LLRGNGGraph<T,S>(*_graphptr);
		pruned_graph->rmNode (i);
		calculateModelEfficiency (pruned_graph, i);

		T change = model_complexity_change + model_efficiency_const * (pruned_graph->model_efficiency - _graphptr->model_efficiency);
		if (change < min_change)
		{
			
			*(this->out) << "eff pruned: " << pruned_graph->model_efficiency << std::endl;
			*(this->out) << "change: " << change << std::endl;
			dislocated_node = i;
			min_change = change;
			if (dislocated_node_graphptr != NULL)
				delete dislocated_node_graphptr;
			dislocated_node_graphptr = new LLRGNGGraph<T,S> (*pruned_graph);
		}
		delete pruned_graph;
	}
	if (dislocated_node != _graphptr->size())
	{
		delete _graphptr;
		_graphptr = dislocated_node_graphptr;
		calculateMinimumDescriptionLength (false);
		findDislocatedNodesStableGraph ();
	}
	this->graphptr = _graphptr;
	this->_graphModulptr = _graphptr;

	return;
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

/** \brief set mean distance calculation mode
    \param mode mean distance calculation mode
*/
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::setMeanDistanceMode (unsigned int mode)
{
	mean_distance_mode = mode;
	_graphptr->setMeanDistanceMode (mode);
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
void LLRGNGAlgorithm<T,S>::updateMinimalGraphMDL (bool calculate_model_efficiency)
{

	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]);
		// *(this->out) << "node " << i << " age: " << node->age << std::endl;
		*(this->out) << "node " << i << " learning q.: " << node->learning_quality << std::endl;
		*(this->out) << "node " << i << " last avg e.: " << node->last_avgerror << std::endl;
		*(this->out) << "node " << i << " prev avg e.: " << node->prev_avgerror << std::endl;
		// *(this->out) << "node " << i << " repulsion: " << node->repulsion << std::endl;
		if (mean_distance_mode == harmonic)
			*(this->out) << "node " << i << " restricting d.: " << node->restricting_distance << std::endl;
		// *(this->out) << "node " << i << " insertion c.: " << node->insertion_criterion << std::endl;
		*(this->out) << "node " << i << " learning rate: " << node->learning_rate << std::endl;
		// *(this->out) << "node " << i << " insertion q.: " << node->insertion_quality << std::endl;
	}
	*(this->out) << "Nr. of nodes now: " <<  _graphptr->size() << std::endl;

	*(this->out) << "mdl: " << calculateMinimumDescriptionLength (calculate_model_efficiency);

	if (min_mdl > mdl )
	{
		*(this->out) << "\t\t\t\t <-- is minimal!" << std::endl;
		min_mdl = mdl;
		delete min_mdl_graphptr;
		min_mdl_graphptr = new LLRGNGGraph<T,S>(*_graphptr);
		// min_mdl_graphptr = new UGraph<T,S>(dynamic_cast<UGraph<T,S>& >(*_graphptr));
		last_epoch_mdl_reduction = epoch;
	}
	else
		*(this->out) << std::endl;
	

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
	// Vector<T> low_limits = _graphptr->getLowLimits();
	// Vector<T> high_limits = _graphptr->getHighLimits();
	delete _graphptr;
	// _graphptr = new LLRGNGGraph<T,S>(this->getDimension());
	// _graphptr->setLowLimits (low_limits);
	// _graphptr->setHighLimits (high_limits);
	// _graphptr->setNodes (min_mdl_graphptr->getNodes());
	_graphptr = new LLRGNGGraph<T,S>(*min_mdl_graphptr);
	findDislocatedNodesStableGraph ();
	stable_graph = true;

	if (visualizing)
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
	T avg_values = 0;
	for(unsigned int t = 0; t < this->size(); t++)
		for (unsigned int i=0; i<this->getDimension(); i++)
			avg_values += (*this)[t][i];
	avg_values /= T (this->size() * this->getDimension());
	value_range = avg_values - this->minValue ();

	*(this->out) << "value range: " << value_range << std::endl;

	bits_vector = ceil (log2(value_range / data_accuracy));

}

/** \brief calculate initial restricting distances for every node
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
			node->restricting_distance += 1.0 / (this->metric (node->weight, (*this)[t]) + 0.01);
		}
		node->restricting_distance /= this->size();
		node->restricting_distance = 1.0 / node->restricting_distance;
		node->prev_restricting_distance = node->restricting_distance;
	}
	
}


/** \brief Get average error among all graph nodes
 */
template<typename T, typename S>
T LLRGNGAlgorithm<T,S>::getAvgError ()
{
	T avgerror = 0;
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]);
		avgerror += node->last_avgerror;
	}
	avgerror /= _graphptr->size ();
	return avgerror;
}

/** \brief check if restricting distances are overflowed and if so recalculate them
 */
template<typename T, typename S>
void LLRGNGAlgorithm<T,S>::checkOverflowedRestrictingDistances ()
{
	for (unsigned int i=0; i < _graphptr->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*_graphptr)[i]);
		if (isnan(node->restricting_distance))
		{
			calculateInitialRestrictingDistances ();
			return;
		}
	}

}

template<typename T, typename S>
template<class Archive>
void 
LLRGNGAlgorithm<T,S>::serialize(Archive & ar, const unsigned int /* file_version */) 
{
	ar & boost::serialization::base_object<NeuralGas<T, S> >(*this);
	ar & BOOST_SERIALIZATION_NVP(_graphptr);
	ar & BOOST_SERIALIZATION_NVP(data_accuracy);
	ar & BOOST_SERIALIZATION_NVP(mdl);
	ar & BOOST_SERIALIZATION_NVP(min_mdl);
	ar & BOOST_SERIALIZATION_NVP(max_epochs_error_reduction);
	ar & BOOST_SERIALIZATION_NVP(max_epochs_mdl_reduction);
	ar & BOOST_SERIALIZATION_NVP(model_efficiency_const);
	ar & BOOST_SERIALIZATION_NVP(mean_distance_mode);
}

} // namespace neuralgas

#endif
