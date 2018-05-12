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
* \file LLRGNGGraph.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef LLRGNGGRAPH_H
#define LLRGNGGRAPH_H

#include <GrowingNeuralGas/GNGModulGraph.h>
#include <tools/metrics.h>

namespace neuralgas {

/** \class LLRGNGNode
 *  \brief The derived node for a LLRGNGGraph
 *
 *  The node permits having a quality measure for learning when inserting and deleting nodes,
 *  calculated from a last mean error and a previous mean error in a similar way as in the paper
 *  "Life-long learning cell structures -- continuosly learning without catastrophic interference"
 *  by Fred H. Hamker,
 *  but using the average error variables proposed in
 * "Intrinsic Motivation Systems for Autonomous Mental Development"
 * by Oudeyer, Kaplan, Hafner.
 *
 */
template<typename T,typename S>
struct LLRGNGNode : Base_Node<T,S>
{
	friend class boost::serialization::access;
	// default \p LLRGNGNode cto
	LLRGNGNode ();
	// default \p LLRGNGNode dto
	virtual ~LLRGNGNode ();
	/// new operator overloading
	static inline void* operator new( std::size_t sz )
	{ return pool.allocate() ; }
	/// delete operator overloading
	static inline void operator delete( void* p )
	{ pool.deallocate( static_cast<LLRGNGNode<T,S>* >(p) ) ; }
	/// errors vector for calculating different parameters. Its size is
	/// set by the maximum allowable error window
	std::vector<T> errors;
	/// errors vector for each dimension
	std::vector<Vector<T> > dim_errors;
	// calculate \p learning_quality measure
	void calculateLearningQuality ();
	/// quality measure for learning
	T learning_quality;
	// calculate \p insertion_quality measure
	void calculateInsertionQuality ();
	/// quality measure for insertion
	T insertion_quality;
	// calculate \p insertion_criterion
	void calculateInsertionCriterion ();
	/// insertion criterion
	T insertion_criterion;
	// calculate \p prev_avgerror and \p last_avgerror
	void updateAvgError (T&, Vector<T>&, const unsigned int&, const unsigned int&, const unsigned int&);
	// update \p restricting_distance
	void updateRestrictingDistance (T);
	/// previous mean error counter
	T prev_avgerror;
	/// last mean error counter
	T last_avgerror;
	/// restricting distance parameter to avoid outlier influence
	T restricting_distance;
	/// previous restricting distance parameter
	T prev_restricting_distance;
	/// last minimal mean error counter found
	T min_last_avgerror;
	/// last epoch where avg error was reduced
	unsigned int last_epoch_improvement;
	/// last mean error for each dimension
	Vector<T> dim_last_avgerror;
	// decrease node age
	void decreaseAge (unsigned int);
	/// age of the node
	T age;
	// update \p learning_rate
	void updateLearningRate (T&, T&);
	/// learning rate of the node
	T learning_rate;
	/// counter for the nr. of items in the receptive field of this node
	unsigned int items_counter;
	/// model efficiency contribution to total MDL
	T efficiency;
	/// repulsion constant for updating weights
	// T repulsion;
	/// mode for calculating mean distances
	unsigned int mean_distance_mode;
	/// data set indices that the node covers (used for calculating mdl values and for active learning)
        std::vector< unsigned int, boost::pool_allocator<unsigned int> > data;
	// memory pool for node objects
        static boost::fast_pool_allocator<LLRGNGNode<T,S> > pool;
private:
	template<class Archive>
	void serialize(Archive & ar, const unsigned int) {
		ar & boost::serialization::base_object<Base_Node<T, S> >(*this);
		ar & BOOST_SERIALIZATION_NVP(mean_distance_mode);
	}
};

/// \brief default \p LLRGNGNode cto
template<typename T, typename S>
LLRGNGNode<T,S>::LLRGNGNode () :
	last_epoch_improvement (0),
	age (1),
	learning_rate (0),
	items_counter (0),
	efficiency (0)
{
}

/// \brief default \p LLRGNGNode dto
template<typename T, typename S>
LLRGNGNode<T,S>::~LLRGNGNode ()
{
	this->edges.clear();
}


/// \brief calculate \p learning_quality measure
template<typename T, typename S>
void LLRGNGNode<T,S>::calculateLearningQuality ()
{
	// learning_quality = (last_avgerror /*+ 1*/) / (prev_avgerror /*+ 1*/);
	// learning_quality = (last_avgerror - prev_avgerror);
	// learning_quality = (last_avgerror / prev_avgerror - 1) * (last_avgerror / prev_avgerror - 1);
	// learning_quality = exp (last_avgerror / prev_avgerror - 1);
	learning_quality = - (log (last_avgerror) - log( prev_avgerror));

}

/// \brief calculate \p insertion_quality measure
template<typename T, typename S>
void LLRGNGNode<T,S>::calculateInsertionQuality ()
{
	insertion_quality = last_avgerror;

}

/// \brief calculate \p insertion_criterion
template<typename T, typename S>
void LLRGNGNode<T,S>::calculateInsertionCriterion ()
{
	insertion_criterion = insertion_quality /*- age*/;
}



/** \brief update \p learning_quality measure
 *  \param adaptation_threshold cut-off value of learning
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::updateLearningRate (T& adaptation_threshold, T& default_rate)
{

	T rate_decision_boundary;

	rate_decision_boundary = learning_quality / (1 + adaptation_threshold) /*+ age*/;
	// rate_decision_boundary = learning_quality;

	// if (rate_decision_boundary < 0)
	// 	learning_rate = 0;
	// else if (rate_decision_boundary <= 1)
	// 	learning_rate = rate_decision_boundary * default_rate;
	// else
	//         learning_rate = default_rate;
	if (rate_decision_boundary > 1)
		learning_rate = default_rate;
	else if (rate_decision_boundary > 0.1)
		learning_rate = rate_decision_boundary * default_rate;	
	else if (rate_decision_boundary <= 0.1)
		learning_rate = 0.1 * default_rate;
	// else
	// 	learning_rate = rate_decision_boundary * default_rate;
	
}


/** \brief calculate \p last_avgerror and \p prev_avgerror,
 *         using the strategy in Oudeyer et al. Harmonic
 *         distances are used.
 *  \param last_error last calculated distance to some data item
 *  \param smoothing smoothing window constant
 *  \param timewindow error time window constant
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::updateAvgError (T& last_error, Vector<T>& dim_last_error, const unsigned int& smoothing, const unsigned int& timewindow, const unsigned int& max_errors_size)
{
	if (errors.size() == max_errors_size)
	{
		errors.erase (errors.begin());
		dim_errors.erase (dim_errors.begin());
	}

	errors.push_back (last_error);
	dim_errors.push_back (dim_last_error);
	
	unsigned int errors_size = errors.size();
	
	T timewindow_ratio = timewindow / T(smoothing + timewindow);
	T smoothing_ratio = smoothing / T(smoothing + timewindow);

	unsigned int smoothing_last = smoothing + 1;
	unsigned int smoothing_prev = smoothing + 1;

	unsigned int windowbegin_prev_avgerror;
	unsigned int windowlast_prev_avgerror;
	unsigned int windowbegin_last_avgerror;
	if (errors_size <= smoothing + timewindow) {
		windowbegin_prev_avgerror = 0;
	}
	else
		windowbegin_prev_avgerror = errors_size -1 - (timewindow + smoothing);

	if (errors_size <= timewindow)
		windowlast_prev_avgerror = (unsigned int)ceil((errors_size - 1) * smoothing_ratio);
	else
		windowlast_prev_avgerror = errors_size - 1 - timewindow;

	
	if (errors_size <= smoothing) {
		windowbegin_last_avgerror = (unsigned int)ceil((errors_size - 1) * timewindow_ratio);
		smoothing_last = smoothing_prev = errors_size - windowbegin_last_avgerror;
		// smoothing_last = errors_size - windowbegin_last_avgerror;
	}
	else
		windowbegin_last_avgerror = errors_size - 1 - smoothing;

	assert (windowbegin_prev_avgerror >= 0 && windowlast_prev_avgerror >= 0 && windowbegin_last_avgerror >= 0);
	prev_avgerror = 0.0;
	last_avgerror = 0.0;

	if (mean_distance_mode == harmonic)
	{
		for (unsigned int i=windowbegin_last_avgerror; i < errors_size; i++) 
			last_avgerror += 1.0 / errors[i];
		for (unsigned int i=windowbegin_prev_avgerror; i<=windowlast_prev_avgerror; i++)
			prev_avgerror += 1.0 / errors[i];	
		
		prev_avgerror = smoothing_prev / prev_avgerror;
		last_avgerror = smoothing_last / last_avgerror;
		
		for (unsigned int i=0; i<dim_last_avgerror.size(); i++)
		{
			dim_last_avgerror[i] = 0.0;
			for (unsigned int j=windowbegin_last_avgerror; j < errors_size; j++)
				dim_last_avgerror[i] += 1.0 / dim_errors[j][i];
			dim_last_avgerror[i] = smoothing_last / dim_last_avgerror[i];
		}
	}
	else if (mean_distance_mode == arithmetic)
	{
		for (unsigned int i=windowbegin_last_avgerror; i < errors_size; i++) 
			last_avgerror += errors[i];
		for (unsigned int i=windowbegin_prev_avgerror; i<=windowlast_prev_avgerror; i++)
			prev_avgerror += errors[i];	
		
		prev_avgerror = prev_avgerror / smoothing_prev;
		last_avgerror = last_avgerror / smoothing_last;

		for (unsigned int i=0; i<dim_last_avgerror.size(); i++)
		{
			dim_last_avgerror[i] = 0.0;
			for (unsigned int j=windowbegin_last_avgerror; j < errors_size; j++)
				dim_last_avgerror[i] += dim_errors[j][i];
			dim_last_avgerror[i] = dim_last_avgerror[i] / smoothing_last;
		}
	}

	if (min_last_avgerror > last_avgerror)
		min_last_avgerror = last_avgerror;

	// repulsion = last_avgerror * 0.3;
	
	// learningProgressHistory.push_back (-(last_avgerror - prev_avgerror));

	
	// cout << "\tLearning progress: " << endl;
	// cout << "\t" << learningProgressHistory.back() << endl;

}

/** \brief update \p restricting_distance
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::updateRestrictingDistance (T last_error)
{
	//update restricting distance
	prev_restricting_distance = restricting_distance;
	if (last_error >= restricting_distance)
		restricting_distance = 1 / (0.5 * (1.0 / restricting_distance + 1.0 / (last_error + 0.01)));
	else
		restricting_distance = 0.5 * (restricting_distance + last_error + 0.01);


}

/** \brief decrease node age
 *  \param age_time_window time window parameter
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::decreaseAge (unsigned int age_time_window)
{
	age = exp (-1/(T) age_time_window) * age;
}

template<typename T, typename S> class LLRGNGAlgorithm;

/** \class LLRGNGGraph
 *  \brief This class provides some additional learning strategies, some of them proposed in
 *   the algorithm Life-long learning cell structures -- continuosly learning
 *   without catastrophic interference by Fred H. Hamker.
 * and Robust growing neural gas algorithm with application in cluster analysis by Qin and Suganthan.
 *   The algorithm tries to avoid bias/variance issues in the original GNG.
*/
template<typename T, typename S>
class LLRGNGGraph : public GNGModulGraph<T,S>
{
	friend class boost::serialization::access;
public:
	// default cto.
	LLRGNGGraph ();
	/// cto Graph creation (node and edges weights share dimensionality)
	LLRGNGGraph (const unsigned int&, const unsigned int&);
	/// copy constructor
	LLRGNGGraph (const LLRGNGGraph&);
	/// std dto
	virtual ~LLRGNGGraph();
	/// new operator overloading
	static inline void* operator new( std::size_t sz )
	{ return llrgngpool.allocate() ; }
	/// delete operator overloading
	static inline void operator delete( void* p )
	{ llrgngpool.deallocate( static_cast<LLRGNGGraph<T,S>* >(p) ) ; }
	// removes the node given by the index, removes its edges and updates the number of connections of its neighbors
	virtual void rmNode(const unsigned int&); 
	// set time window constants
	void setTimeWindows (unsigned int, unsigned int, unsigned int);
	// reset max errors size
	void setMaxErrorsSize (unsigned int);
	// calculate inherited variables for a node to be inserted between two nodes
	void calculateInheritedParams (const unsigned int, const unsigned int, const unsigned int);
	// calculate long term and short term error for some node
	void updateAvgError (const unsigned int, T&, Vector<T>&);
	// update restricting distance value for some node
	void updateRestrictingDistance (const unsigned int, T);
	// calculate learning quality for some node
	void calculateLearningQuality (const unsigned int);
	// calculate insertion quality for some node
	void calculateInsertionQuality (const unsigned int);
	// calculate insertion criterion for some node
	void calculateInsertionCriterion (const unsigned int);
	// update learning rate for some winner node
	void updateWinnerLearningRate (const unsigned int);
	// update learning rate for some winner-neighboring node
	void updateNeighborLearningRate (const unsigned int);
	// set input adaptation threshold
	void setAdaptationThreshold (T);
	// set initial learning rate constants
	void setLearningRates (T, T);
	// set maximal edge age constant
	void setMaximalEdgeAge (unsigned int&);
	// get maximal edge age constant
	unsigned int getMaximalEdgeAge () const;
	// decrease age of some node
	void decreaseNodeAge (const unsigned int);
	// update items counter in the receptive field for some node
	void increaseItemsCounter (const unsigned int);
	// delete nodes that do not have items in their receptive fields
	bool deleteInactiveNodes (unsigned int&, unsigned int&);
	// find node with less items in their receptive fields
	int findLessItemsNode ();
	// reset items counter and efficiency contribution for MDL for all nodes
	void resetMDLCounters ();
	// get last stored minimal average error for some node
	T getNodeMinLastAvgError (const unsigned int);
	// set last epoch where an error reduction for a node was achieved
	void setLastEpochImprovement (const unsigned int, unsigned int);
	// set mean distance calculation mode
	void setMeanDistanceMode (unsigned int);
	friend class LLRGNGAlgorithm<T,S>;

protected:
	// returns a pointer to a edge of a type that is currently used by the graph
	virtual LLRGNGNode<T,S>* newNode(void);
	/// adaptation threshold constant
	T adaptation_threshold;
	/// initial winner learning rate constant
	T winner_learning_rate;
	/// initial winner-neighbors learning rate constant
	T neighbors_learning_rate;
	/// maximal edge age constant
	unsigned int maximal_edge_age;
	/// smoothing window constant
	unsigned int smoothing_window;
	/// error time window constant
	unsigned int error_time_window;
	/// age time window constant
	unsigned int age_time_window;
	/// maximal size of error vector
	unsigned int max_errors_size;
	/// current model efficiency value
	T model_efficiency;
	/// mode for calculating mean distances
	unsigned int mean_distance_mode;
	// memory pool for graph objects
        static boost::fast_pool_allocator<LLRGNGGraph<T,S> > llrgngpool;
private:
	template<class Archive>
	void serialize(Archive & ar, const unsigned int);
		
};

/** \brief cto Graph creation (node and edges weights share dimensionality)
    \param dim dimensionality
 */
template<typename T, typename S>
LLRGNGGraph<T,S>::LLRGNGGraph (const unsigned int &dim, const unsigned int& max_error_window) :
	Base_Graph<T,S>(dim),
	UGraph<T,S>(dim),
	TGraph<T,S>(dim),
	GNGModulGraph<T,S>(dim),
	adaptation_threshold (0),
	winner_learning_rate (0),
	neighbors_learning_rate (0),
	maximal_edge_age (0),
	smoothing_window (1),
	error_time_window (1),
	age_time_window (1),
	max_errors_size (max_error_window),
	model_efficiency (0),
	mean_distance_mode (harmonic)
{
}

/** \brief Default constructor (Should only be used for serialization)
 */
template<typename T, typename S>
LLRGNGGraph<T,S>::LLRGNGGraph () :
	Base_Graph<T,S>(0),
	UGraph<T,S>(0),
	TGraph<T,S>(0),
	GNGModulGraph<T,S>(0),
	adaptation_threshold (0),
	winner_learning_rate (0),
	neighbors_learning_rate (0),
	maximal_edge_age (0),
	smoothing_window (0),
	error_time_window (0),
	age_time_window (0),
	max_errors_size (81),
	model_efficiency (0),
	mean_distance_mode (harmonic)
{
}


/** \brief copy constructor.
 *   Calls UGraph and Base_Graph customized copy constructors (they are empty)
 *  \param g graph to be copied from. 
 */
template<typename T,typename S>
LLRGNGGraph<T,S>::LLRGNGGraph (const LLRGNGGraph& g) :
	Base_Graph<T,S>(),
	UGraph<T,S>(),
	TGraph<T,S>(),
	GNGModulGraph<T,S>(),
	adaptation_threshold (g.adaptation_threshold),
	winner_learning_rate (g.winner_learning_rate),
	neighbors_learning_rate (g.neighbors_learning_rate),
	maximal_edge_age (g.maximal_edge_age),
	smoothing_window (g.smoothing_window),
	error_time_window (g.error_time_window),
	age_time_window (g.age_time_window),
	max_errors_size (g.max_errors_size),
	model_efficiency (g.model_efficiency),
	mean_distance_mode (g.mean_distance_mode)
{

	this->_dimNode = g._dimNode;
	this->_dimEdge = g._dimEdge;
	this->_metric_to_use = g._metric_to_use;
	this->low_limit = g.low_limit;
	this->high_limit = g.high_limit;
	this->low_limits = g.low_limits;
	this->high_limits = g.high_limits;
	unsigned int gsize = g.size();
	
	for (unsigned int i=0; i < gsize; i++)
	{
		this->addNode ();
		this->_nodes[i]->weight = g._nodes[i]->weight;
		
	}
	for (unsigned int i=0; i < gsize; i++)
		for (unsigned int j=0; j<g._nodes[i]->edges.size(); j++)
			if ( g._nodes[i]->edges[j] != NULL )
			{
				this->addEdge (i, j);
				this->_nodes[i]->edges[j]->weight = g._nodes[i]->edges[j]->weight;
				static_cast<TEdge<S,T>* >(this->_nodes[i]->edges[j])->age = static_cast<TEdge<S,T>* >(g._nodes[i]->edges[j])->age;
			}
	
	for (unsigned int i=0; i < gsize; i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* >(this->_nodes[i]);
		LLRGNGNode<T,S>* copynode = static_cast<LLRGNGNode<T,S>* >(g._nodes[i]);

		node->learning_quality = copynode->learning_quality;
		node->insertion_quality = copynode->insertion_quality;
		node->insertion_criterion = copynode->insertion_criterion;
		node->prev_avgerror = copynode->prev_avgerror;
		node->last_avgerror = copynode->last_avgerror;
		node->restricting_distance = copynode->restricting_distance;
		node->prev_restricting_distance = copynode->prev_restricting_distance;
		// node->repulsion = copynode->repulsion;
		node->min_last_avgerror = copynode->min_last_avgerror;
		node->last_epoch_improvement = copynode->last_epoch_improvement;
		node->age = copynode->age;
		node->learning_rate = copynode->learning_rate;
		node->items_counter = copynode->items_counter;
		node->errors = copynode->errors;
		node->dim_errors = copynode->dim_errors;
		node->efficiency = copynode->efficiency;
		node->mean_distance_mode = copynode->mean_distance_mode;
		node->data = copynode->data;
	}
}

/** \brief dto Graph deletion
 */
template<typename T, typename S>
LLRGNGGraph<T,S>::~LLRGNGGraph ()
{			     
}


/** \brief overriden function from \p Base_Graph
    to create a node of type \p LLRGNGNode */
template<typename T,typename S>
LLRGNGNode<T,S>* LLRGNGGraph<T,S>::newNode(void)
{
	LLRGNGNode<T,S>* n = new LLRGNGNode<T,S>;
	if (this->high_limits.size() != 0 && this->low_limits.size() != 0)
		n->last_avgerror = this->metric (this->high_limits, this->low_limits);
	//default inherited errors (used when initializing reference vectors)
	n->prev_avgerror = n->last_avgerror;
	// n->repulsion = 0.001;
	n->errors.push_back (n->last_avgerror);
	n->min_last_avgerror = n->last_avgerror;
	n->dim_last_avgerror.reserve (this->_dimNode);
	n->dim_last_avgerror.resize (this->_dimNode);
	n->dim_errors.reserve (max_errors_size);
	n->dim_errors.resize (1);
	n->dim_errors[0].reserve (this->_dimNode);
	n->dim_errors[0].resize (this->_dimNode);
	n->mean_distance_mode = mean_distance_mode;
	if (this->high_limits.size() != 0 && this->low_limits.size() != 0)
		for (unsigned int i=0; i<this->_dimNode; i++)
			n->dim_errors[0][i] = (this->high_limits[i] - this->low_limits[i])*(this->high_limits[i] - this->low_limits[i]);
	return n; 
}

/** \brief Removes the node given by the index, removes its edges and updates the number 
 * of connections of its neighbors. Connects the nodes that were its neighbors.
 *
 * \param index is the node that shall be deleted
 */
template<typename T,typename S>
void LLRGNGGraph<T,S>::rmNode(const unsigned int& index)
{
	// function does as follows
	// checks whether the neighbors have an (directed) edge to index and deletes them
	// deletes node in the _nodes array  
 
	unsigned int     nsize=this->size();
	assert ( index < nsize );
	std::vector<unsigned int> neighbors = this->getNeighbors(index);
   
	for(unsigned int i=0; i < neighbors.size(); i++)
		if ( this->_nodes[ neighbors[i] ]->edges[index]!=NULL ) // edge from i to index
		{
			delete  this->_nodes[  neighbors[i] ]->edges[index];
			this->_nodes[ neighbors[i] ]->edges[index]=NULL; // (really needed?)
			this->_nodes[ neighbors[i] ]->num_connections--;
		}  
   
	for(unsigned int i=0; i < nsize; i++)
		this->_nodes[ i ]->edges.erase( this->_nodes[ i ]->edges.begin() + index ); 

	delete this->_nodes[index];                                 // delete ptrs to the nodes (really needed?)
	this->_nodes[index] = NULL;  // (really needed?)
   
	this->_nodes.erase(this->_nodes.begin() + index);                 // erase node within the array

	//reindex
	for(unsigned int i=0; i < neighbors.size(); i++)
		if (neighbors[i] > index)
			neighbors[i]--;
	
	// connect all neighbor nodes
	for(unsigned int i=0; i < neighbors.size(); i++)
		for(unsigned int j=0; j < neighbors.size(); j++)
			if (i != j)
				this->setAge (neighbors[i], neighbors[j], 0.0);
}


/** \brief set window constants. This function should only
    be called once when creating the graph
    \param smoothing smoothing time window constant
    \param longterm error time window constant
    \param age age time window constant */
template<typename T, typename S>
void LLRGNGGraph<T,S>::setTimeWindows (unsigned int smoothing, unsigned int error, unsigned int age)
{
	assert (max_errors_size > error + smoothing);

	smoothing_window = smoothing;
	error_time_window = error;
	age_time_window = age;
}

/** \brief reset max errors size
    \param size
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::setMaxErrorsSize (unsigned int size)
{
	max_errors_size = size;
}


/** \brief calculate inherited variables for a node to be inserted between two nodes
 *  \param index node index
 *  \param first_index 1st node index
 *  \param snd_index 2nd node index
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::calculateInheritedParams (const unsigned int index, const unsigned int first_index, const unsigned int snd_index)
{
	LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* >(this->_nodes[index]);
	LLRGNGNode<T,S>* first_node = static_cast<LLRGNGNode<T,S>* >(this->_nodes[first_index]);
	// LLRGNGNode<T,S>* snd_node = static_cast<LLRGNGNode<T,S>* >(this->_nodes[snd_index]);

	// node->weight = first_node->weight + (snd_node->weight - node->weight) * 0.25;
	// node->weight = 2 * first_node->weight + snd_node->weight / 3;
	// node->weight = 2 * first_node->weight + snd_node->weight / -3;
	// node->weight = (first_node->weight + snd_node->weight) / 2;
	node->weight = first_node->weight + 2 * first_node->dim_last_avgerror;
	// node->prev_avgerror = (first_node->prev_avgerror + snd_node->prev_avgerror) / 2;
	// node->last_avgerror = (first_node->last_avgerror + snd_node->last_avgerror) / 2;
	node->prev_avgerror = node->last_avgerror = 0;
	// node->repulsion = node->last_avgerror * 0.3;

	node->min_last_avgerror = node->last_avgerror;
	first_node->last_avgerror = first_node->prev_avgerror = 0;
	first_node->errors.clear();
	first_node->dim_errors.clear();
	
	assert (node->errors.size() == 1);
	// node->errors.front() = node->last_avgerror;
	node->errors.pop_back();
	node->dim_errors.pop_back();

	
}



/** \brief calculate last and previous mean error for a given node
    \param index node index */
template<typename T, typename S>
void LLRGNGGraph<T,S>::updateAvgError (const unsigned int index, T& last_error, Vector<T>& dim_last_error)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->updateAvgError(last_error, dim_last_error, smoothing_window, error_time_window, max_errors_size);

}

//! \brief update restricting distance value for some node
/*! 
  
  \param index node index
  \param last_error last error to calculate current value
*/
template<typename T, typename S>
void LLRGNGGraph<T,S>::updateRestrictingDistance (const unsigned int index, T last_error)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->updateRestrictingDistance(last_error);

}


/** \brief calculate learning quality for a node
    \param index node index */
template<typename T, typename S>
void LLRGNGGraph<T,S>::calculateLearningQuality (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->calculateLearningQuality();
	
}

/** \brief calculate insertion quality for a node
    \param index node index */
template<typename T, typename S>
void LLRGNGGraph<T,S>::calculateInsertionQuality (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->calculateInsertionQuality();
	
}

/** \brief calculate insertion criterion for a node
    \param index node index */
template<typename T, typename S>
void LLRGNGGraph<T,S>::calculateInsertionCriterion (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->calculateInsertionCriterion();
	
}

/** \brief set initial learning rate constants
 *  \param winner initial winner learning rate constant
 *  \param neighbors initial winner-neighbors learning rate constant
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::setLearningRates (T winner, T neighbors)
{
	winner_learning_rate = winner;
	neighbors_learning_rate = neighbors;
}


/** \brief update learning rate for a winner node
 *  \param index node index
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::updateWinnerLearningRate (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->updateLearningRate(adaptation_threshold, winner_learning_rate);	
}

/** \brief update learning rate for a winner-neighboring node
 *  \param index node index
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::updateNeighborLearningRate (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->updateLearningRate(adaptation_threshold, neighbors_learning_rate);	
}


/** \brief set adaptation threshold
    \param threshold given adaptation threshold */
template<typename T, typename S>
void LLRGNGGraph<T,S>::setAdaptationThreshold (T threshold)
{
	adaptation_threshold = threshold;	
}

/* \brief set maximal edge age constant
 * \param age maximal edge age constant
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::setMaximalEdgeAge (unsigned int& age)
{
	maximal_edge_age = age;
}

/* \brief get \p maximal_edge_age constant
 */
template<typename T, typename S>
unsigned int LLRGNGGraph<T,S>::getMaximalEdgeAge () const
{
	return maximal_edge_age;
}

/* \brief decrease age of a node
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::decreaseNodeAge (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->decreaseAge (age_time_window);
}

/** \brief update items counter for the receptive field of some node
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::increaseItemsCounter (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->increaseItemsCounter ();
}

/* \brief delete nodes that are do not have items in their receptive fields
   \param winner winner node index that needs to be reindexed if necessary
   \param snd_winner winner node index that needs to be reindexed if necessary
 */
template<typename T, typename S>
bool LLRGNGGraph<T,S>::deleteInactiveNodes (unsigned int& winner, unsigned int& snd_winner)
{
	bool node_deleted = false;
	if (this->size() <= 2)
		return false;

	for (unsigned int i=0; i < this->size(); i++)
	{
		if (static_cast<LLRGNGNode<T,S>* > (this->_nodes[i])->items_counter == 0)
		{
			if (this->size() > 2 && i != winner)
			{
				// std::cout << "deleting node " << i << std::endl;
				rmNode (i);
				node_deleted = true;
				if (winner > i)
					winner--;
				if (snd_winner > i)
					snd_winner--;
				else if (snd_winner == i)
					snd_winner = this->size();
				i--;
			}
		}
		
	}
	return node_deleted;
}

/* \brief find node with less items in their receptive fields
 */
template<typename T, typename S>
int LLRGNGGraph<T,S>::findLessItemsNode ()
{
	unsigned int min_value = static_cast<LLRGNGNode<T,S>* > (this->_nodes[0])->activations_counter;
	int n = 0;
	
 	for (unsigned int i=1; i < this->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (this->_nodes[i]);
		if (node->items_counter < min_value)
		{
			min_value = node->items_counter;
			n = i;
		}
	}
	return n;
}

/** \brief reset items counter for all nodes and efficiency contribution to MDL
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::resetMDLCounters ()
{
	model_efficiency = 0;
	for (unsigned int i=0; i < this->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (this->_nodes[i]);
		node->items_counter = 0;
		node->efficiency = 0;
		node->data.clear ();
	}
}

/* \brief get last stored minimal average error for some node
 */
template<typename T, typename S>
T LLRGNGGraph<T,S>::getNodeMinLastAvgError (const unsigned int index)
{
	return static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->min_last_avgerror;
}

/* \brief set last epoch where an error reduction for a node was achieved
   \param index node index
   \param epoch current training epoch
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::setLastEpochImprovement (const unsigned int index, unsigned int epoch)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->last_epoch_improvement = epoch;
}

/** \brief set mean distance calculation mode
    \param mode mean distance calculation mode
*/
template<typename T, typename S>
void LLRGNGGraph<T,S>::setMeanDistanceMode (unsigned int mode)
{
	mean_distance_mode = mode;
	for (unsigned int i=0; i < this->size(); i++)
	{
		LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (this->_nodes[i]);
		node->mean_distance_mode = mode;
	}
}

template<typename T, typename S>
template<class Archive>
void 
LLRGNGGraph<T,S>::serialize(Archive & ar, const unsigned int /* file_version */) 
{
	ar & boost::serialization::base_object<Base_Graph<T, S> >(*this);
	ar & BOOST_SERIALIZATION_NVP(adaptation_threshold);
	ar & BOOST_SERIALIZATION_NVP(winner_learning_rate);
	ar & BOOST_SERIALIZATION_NVP(neighbors_learning_rate);
	ar & BOOST_SERIALIZATION_NVP(maximal_edge_age);
	ar & BOOST_SERIALIZATION_NVP(smoothing_window);
	ar & BOOST_SERIALIZATION_NVP(error_time_window);
	ar & BOOST_SERIALIZATION_NVP(age_time_window);
	ar & BOOST_SERIALIZATION_NVP(max_errors_size);
	ar & BOOST_SERIALIZATION_NVP(model_efficiency);
	ar & BOOST_SERIALIZATION_NVP(mean_distance_mode);

}


} // namespace neuralgas

#endif
