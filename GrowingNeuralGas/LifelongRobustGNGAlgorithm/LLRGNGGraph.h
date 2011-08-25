/** 
* \class LLRGNGGraph
* \author Sergio Roa
* \author Manuel Noll
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


/** \brief The derived node for a LLRGNGGraph
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
	// default \p LLRGNGNode cto
	LLRGNGNode ();
	// default \p LLRGNGNode dto
	~LLRGNGNode ();
	/// errors vector for calculating different parameters. Its size is
	/// set by the maximum allowable error window
	std::vector<T> errors;
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
	void updateAvgError (T, const unsigned int&, const unsigned int&, const unsigned int&);
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
	// decrease node age
	void decreaseAge (unsigned int);
	/// age of the node
	T age;
	// update \p learning_rate
	void updateLearningRate (T&, T&);
	/// learning rate of the node
	T learning_rate;
	// update activations counter
	void increaseActivationsCounter ();
	/// counter for the nr. of activations in a training epoch
	unsigned int activations_counter;
	/// repulsion constant for updating weights
	T repulsion;
};

/// \brief default \p LLRGNGNode cto
template<typename T, typename S>
LLRGNGNode<T,S>::LLRGNGNode () :
	age (1),
	learning_rate (0),
	activations_counter (0)
{
}

/// \brief default \p LLRGNGNode dto
template<typename T, typename S>
LLRGNGNode<T,S>::~LLRGNGNode ()
{
	errors.clear ();
}


/// \brief calculate \p learning_quality measure
template<typename T, typename S>
void LLRGNGNode<T,S>::calculateLearningQuality ()
{
	// learning_quality = (last_avgerror /*+ 1*/) / (prev_avgerror /*+ 1*/);
	// learning_quality = (last_avgerror - prev_avgerror);
	// learning_quality = (last_avgerror / prev_avgerror - 1) * (last_avgerror / prev_avgerror - 1);
	learning_quality = exp (last_avgerror / prev_avgerror - 1);
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
	else
		learning_rate = rate_decision_boundary * default_rate;
	
}


/** \brief calculate \p last_avgerror and \p prev_avgerror,
 *         using the strategy in Oudeyer et al. Harmonic
 *         distances are used.
 *  \param last_error last calculated distance to some data item
 *  \param smoothing smoothing window constant
 *  \param timewindow error time window constant
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::updateAvgError (T last_error, const unsigned int& smoothing, const unsigned int& timewindow, const unsigned int& max_errors_size)
{
	errors.push_back (last_error);
	if (errors.size() > max_errors_size)
		errors.erase (errors.begin());


	unsigned int errors_size = errors.size();
	
	T timewindow_ratio = timewindow / T(smoothing + timewindow);
	T smoothing_ratio = smoothing / T(smoothing + timewindow);

	unsigned int smoothing_last = smoothing + 1;
	unsigned int smoothing_prev = smoothing + 1;

	unsigned int windowbegin_prev_avgerror;
	unsigned int windowlast_prev_avgerror;
	unsigned int windowbegin_last_avgerror;
	if (errors.size () <= smoothing + timewindow) {
		windowbegin_prev_avgerror = 0;
	}
	else
		windowbegin_prev_avgerror = errors_size -1 - (timewindow + smoothing);

	if (errors.size () <= timewindow)
		windowlast_prev_avgerror = (unsigned int)ceil((errors_size - 1) * smoothing_ratio);
	else
		windowlast_prev_avgerror = errors_size - 1 - timewindow;

	
	if (errors.size () <= smoothing) {
		windowbegin_last_avgerror = (unsigned int)ceil((errors_size - 1) * timewindow_ratio);
		smoothing_prev = errors_size - windowbegin_last_avgerror;
		smoothing_last = errors_size - windowbegin_last_avgerror;
	}
	else
		windowbegin_last_avgerror = errors_size - 1 - smoothing;

	assert (windowbegin_prev_avgerror >= 0 && windowlast_prev_avgerror >= 0 && windowbegin_last_avgerror >= 0);
	prev_avgerror = 0.0;
	last_avgerror = 0.0;

	
	for (unsigned int i=windowbegin_last_avgerror; i < errors_size; i++) 
		last_avgerror += 1.0 / errors[i];
	for (unsigned int i=windowbegin_prev_avgerror; i<windowlast_prev_avgerror; i++)
		prev_avgerror += 1.0 / errors[i];	
	
	prev_avgerror /= smoothing_prev;
	last_avgerror /= smoothing_last;
	prev_avgerror = 1.0 / prev_avgerror;
	last_avgerror = 1.0 / last_avgerror;

	if (min_last_avgerror > last_avgerror)
		min_last_avgerror = last_avgerror;

	repulsion = last_avgerror * 0.2;

	// learningProgressHistory.push_back (-(last_avgerror - prev_avgerror));

	
	// cout << "\tLearning progress: " << endl;
	// cout << "\t" << learningProgressHistory.back() << endl;

}

template<typename T, typename S>
void LLRGNGNode<T,S>::updateRestrictingDistance (T last_error)
{
	//update restricting distance
	prev_restricting_distance = restricting_distance;
	if (last_error >= restricting_distance)
		restricting_distance = 1 / (0.5 * (1.0 / restricting_distance + 1.0 / last_error));
	else
		restricting_distance = 0.5 * (restricting_distance + last_error);


}

/** \brief decrease node age
 *  \param age_time_window time window parameter
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::decreaseAge (unsigned int age_time_window)
{
	age = exp (-1/(T) age_time_window) * age;
}

/** \brief update \p activations_counter
 */
template<typename T, typename S>
void LLRGNGNode<T,S>::increaseActivationsCounter ()
{
	activations_counter++;
}

/** \brief LLRGNGGraph provides some additional learning strategies, some of them proposed in
 *   the algorithm Life-long learning cell structures -- continuosly learning
 *   without catastrophic interference by Fred H. Hamker.
 * and Robust growing neural gas algorithm with application in cluster analysis by Qin and Suganthan.
 *   The algorithm tries to avoid bias/variance issues in the original GNG.
*/
template<typename T, typename S>
class LLRGNGGraph : public GNGModulGraph<T,S>
{
public:
	/// cto Graph creation (node and edges weights share dimensionality)
	LLRGNGGraph (const unsigned int&, const unsigned int&);
	/// copy constructor
	LLRGNGGraph (const LLRGNGGraph&);
	/// std dto
	~LLRGNGGraph(){}
	// set time window constants
	void setTimeWindows (unsigned int, unsigned int, unsigned int);
	// calculate inherited variables for a node to be inserted between two nodes
	void calculateInheritedParams (const unsigned int, const unsigned int, const unsigned int);
	// calculate long term and short term error for some node
	void updateAvgError (const unsigned int, T);
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
	void setAdaptationThreshold (T&);
	// set initial learning rate constants
	void setLearningRates (T&, T&);
	// set maximal edge age constant
	void setMaximalEdgeAge (unsigned int&);
	// get maximal edge age constant
	unsigned int getMaximalEdgeAge () const;
	// decrease age of some node
	void decreaseNodeAge (const unsigned int);
	// update activations counter for some node
	void increaseActivationsCounter (const unsigned int);
	// reset activations counters for all nodes
	void resetActivationsCounters ();
	// get last stored minimal average error for some node
	T getNodeMinLastAvgError (const unsigned int);
	// set last epoch where an error reduction for a node was achieved
	void setLastEpochImprovement (const unsigned int, unsigned int);
private:
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
	const unsigned int max_errors_size;
};

/** \brief cto Graph creation (node and edges weights share dimensionality)
    \param dim dimensionality
 */
template<typename T, typename S>
LLRGNGGraph<T,S>::LLRGNGGraph (const unsigned int &dim, const unsigned int& max_error_window = 200) :
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
	max_errors_size (max_error_window)
{
}

/** \brief copy constructor.
 *   Calls UGraph and Base_Graph customized copy constructors (they are empty)
 *  \param g graph to be copied from. 
 */
template<typename T,typename S>
LLRGNGGraph<T,S>::LLRGNGGraph (const LLRGNGGraph& g) :
	Base_Graph<T,S>(g),
	UGraph<T,S>(g),
	TGraph<T,S>(g),
	GNGModulGraph<T,S>(g),
	adaptation_threshold (g.adaptation_threshold),
	winner_learning_rate (g.winner_learning_rate),
	neighbors_learning_rate (g.neighbors_learning_rate),
	maximal_edge_age (g.maximal_edge_age),
	smoothing_window (g.smoothing_window),
	error_time_window (g.error_time_window),
	age_time_window (g.age_time_window),
	max_errors_size (g.max_errors_size)
{

	this->_dimNode = g._dimNode;
	this->_dimEdge = g._dimEdge;
	this->_metric_to_use = g._metric_to_use;
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
		node->repulsion = copynode->repulsion;
		node->min_last_avgerror = copynode->min_last_avgerror;
		node->last_epoch_improvement = copynode->last_epoch_improvement;
		node->age = copynode->age;
		node->learning_rate = copynode->learning_rate;
		node->activations_counter = copynode->activations_counter;
		node->errors = copynode->errors;
	}
}



/** \brief overriden function from \p Base_Graph
    to create a node of type \p LLRGNGNode */
template<typename T,typename S>
LLRGNGNode<T,S>* LLRGNGGraph<T,S>::newNode(void)
{
	LLRGNGNode<T,S>* n = new LLRGNGNode<T,S>;
	if (this->high_limits.size() != 0 && this->low_limits.size() != 0)
		n->last_avgerror = metric (this->high_limits, this->low_limits);
	else
		n->last_avgerror = (this->high_limit - this->low_limit) * (this->high_limit - this->low_limit);		
	//default inherited errors (used when initializing reference vectors)
	n->prev_avgerror = n->last_avgerror;
	n->repulsion = 0.001;
	n->errors.push_back (n->last_avgerror);
	n->min_last_avgerror = n->last_avgerror;
	return n; 
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
	LLRGNGNode<T,S>* snd_node = static_cast<LLRGNGNode<T,S>* >(this->_nodes[snd_index]);

	// node->weight = first_node->weight + (snd_node->weight - node->weight) * 0.25;
	// node->weight = 2 * first_node->weight + snd_node->weight / 3;
	// node->weight = 2 * first_node->weight + snd_node->weight / -3;
	node->weight = (first_node->weight + snd_node->weight) / 2;
	node->prev_avgerror = (first_node->prev_avgerror + snd_node->prev_avgerror) / 2;
	node->last_avgerror = (first_node->last_avgerror + snd_node->last_avgerror) / 2;
	node->repulsion = node->last_avgerror * 0.2;

	node->min_last_avgerror = node->last_avgerror;
	
	assert (node->errors.size() == 1);
	node->errors.front() = node->last_avgerror;

	
}



/** \brief calculate last and previous mean error for a given node
    \param index node index */
template<typename T, typename S>
void LLRGNGGraph<T,S>::updateAvgError (const unsigned int index, T last_error)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->updateAvgError(last_error, smoothing_window, error_time_window, max_errors_size);

}

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
void LLRGNGGraph<T,S>::setLearningRates (T& winner, T& neighbors)
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
void LLRGNGGraph<T,S>::setAdaptationThreshold (T& threshold)
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

/** \brief update activations counter for some node
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::increaseActivationsCounter (const unsigned int index)
{
	static_cast<LLRGNGNode<T,S>* > (this->_nodes[index])->increaseActivationsCounter ();
}

/** \brief reset activations counter for all nodes
 */
template<typename T, typename S>
void LLRGNGGraph<T,S>::resetActivationsCounters ()
{
	for (unsigned int i=0; i < this->size(); i++)
		static_cast<LLRGNGNode<T,S>* > (this->_nodes[i])->activations_counter = 0;
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

} // namespace neuralgas

#endif