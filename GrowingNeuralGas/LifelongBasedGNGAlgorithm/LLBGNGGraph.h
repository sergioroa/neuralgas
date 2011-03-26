/** 
* \class LLBGNGGraph
* \author Sergio Roa
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef LLBGNGGRAPH_H
#define LLBGNGGRAPH_H

#include <GrowingNeuralGas/GNGModulGraph.h>
#include <queue>
#include <tools/metrics.h>

namespace neuralgas {


/** \brief The derived node for a LLBGNGGraph
 *
 *  The node permits having a quality measure for learning when inserting and deleting nodes,
 *  calculated from a short-term error and a long-term error as explained in the paper
 *  "Life-long learning cell structures -- continuosly learning without catastrophic interference
 *  by Fred H. Hamker.
 *
 */
template<typename T,typename S>
struct LLBGNGNode : Base_Node<T,S>
{
	// default \p LLBGNGNode cto
	LLBGNGNode ();
	// default \p LLBGNGNode dto
	~LLBGNGNode ();
	/// errors queue for calculating different parameters. Its size is
	/// set by the maximum allowable long term error window
	std::queue<T> errors;
	// calculate \p learning_quality measure
	void calculateLearningQuality ();
	/// quality measure for learning
	T learning_quality;
	// calculate \p insertion_quality measure
	void calculateInsertionQuality (T&);
	/// quality measure for insertion
	T insertion_quality;
	// calculate \p insertion_criterion
	void calculateInsertionCriterion ();
	/// insertion criterion
	T insertion_criterion;
	// calculate \p shortterm_avgerror and \p longterm_avgerror
	void updateAvgError (T, const T&, const T&);
	/// short-term error counter
	T shortterm_avgerror;
	/// long-term error counter
	T longterm_avgerror;
	/// inherited error calculating during node insertion
	T inherited_error;
	/// An insertion is only allowed if the long-term error exceeds this threshold
	T insertion_threshold;
	// decrease node age
	void decreaseAge (unsigned int);
	/// age of the node
	T age;
	// update \p learning_rate
	void updateLearningRate (T&, T&);
	/// learning rate of the node
	T learning_rate;
	/// local similarity of weights
	T local_similarity;
	// update activations counter
	void increaseActivationsCounter ();
	/// counter for the nr. of activations in a training epoch
	unsigned int activations_counter;
};

/// \brief default \p LLBGNGNode cto
template<typename T, typename S>
LLBGNGNode<T,S>::LLBGNGNode () :
	insertion_threshold (0),
	age (1),
	learning_rate (0),
	local_similarity (0),
	activations_counter (0)
{
}

/// \brief default \p LLBGNGNode dto
template<typename T, typename S>
LLBGNGNode<T,S>::~LLBGNGNode ()
{
	errors.clear ();
}


/// \brief calculate \p learning_quality measure
template<typename T, typename S>
void LLBGNGNode<T,S>::calculateLearningQuality ()
{
	learning_quality = (shortterm_avgerror + 1) / (longterm_avgerror + 1);
}

/// \brief calculate \p insertion_quality measure
template<typename T, typename S>
void LLBGNGNode<T,S>::calculateInsertionQuality (T& insertion_tolerance)
{
	insertion_quality = longterm_avgerror - insertion_threshold * (1 + insertion_tolerance);
}

/// \brief calculate \p insertion_criterion
template<typename T, typename S>
void LLBGNGNode<T,S>::calculateInsertionCriterion ()
{
	insertion_criterion = insertion_quality - age;
}



/** \brief update \p learning_quality measure
 *  \param adaptation_threshold cut-off value of learning
 */
template<typename T, typename S>
void LLBGNGNode<T,S>::updateLearningRate (T& adaptation_threshold, T& default_rate)
{
	learning_quality = (shortterm_avgerror + 1) / (longterm_avgerror + 1);

	T rate_decision_boundary;

	rate_decision_boundary = learning_quality / (1 + adaptation_threshold) + age -1;

	if (rate_decision_boundary < 0)
		learning_rate = 0;
	else if (rate_decision_boundary <= 1)
		learning_rate = rate_decision_boundary * default_rate;
	else
		learning_rate = default_rate;
}


/** \brief calculate \p shortterm_avgerror and \p longterm_avgerror
 *  \param last_error last calculated distance to some data item
 *  \param shortterm_window short term window constant
 *  \param longterm_window long term window constant
 */
template<typename T, typename S>
void LLBGNGNode<T,S>::updateAvgError (T last_error, const T& shortterm_window, const T& longterm_window)
{
	errors.push (last_error);
	T front_error = errors.front();
	T shortterm_window_scaled = shortterm_window;
	T longterm_window_scaled = longterm_window;
	if (errors.size() < longterm_window)
	{
		T errorwindows_ratio = longterm_window / shortterm_window;
		// std::cout << "errowindows ratio: " << errorwindows_ratio << std::endl;
		shortterm_window_scaled = errors.size() / errorwindows_ratio;
		longterm_window_scaled = errors.size();
 		// std::cout << "shorterm window: " << shortterm_window_scaled << std::endl;
		// std::cout << "longterm window: " << longterm_window_scaled << std::endl;
	}
	
	shortterm_avgerror += ( (last_error - front_error) / shortterm_window_scaled);
	// std::cout << "errors size: " << errors.size() << std::endl;
	// std::cout << "shortterm update: " << shortterm_avgerror << std::endl;
	longterm_avgerror += ( (last_error - front_error) / longterm_window_scaled);

	if (errors.size() > longterm_window)
		errors.pop ();
}

/** \brief decrease node age
 *  \param age_time_window time window parameter
 */
template<typename T, typename S>
void LLBGNGNode<T,S>::decreaseAge (unsigned int age_time_window)
{
	age = exp (-1/(T) age_time_window) * age;
}

/** \brief update \p activations_counter
 */
template<typename T, typename S>
void LLBGNGNode<T,S>::increaseActivationsCounter ()
{
	activations_counter++;
}

/** \brief LLBGNGGraph provides some additional learning strategies proposed in
 *   the algorithm Life-long learning cell structures -- continuosly learning
 *   without catastrophic interference by Fred H. Hamker.
 *   The algorithm tries to avoid bias/variance issues in the original GNG.
*/
template<typename T, typename S>
class LLBGNGGraph : public GNGModulGraph<T,S>
{
public:
	/// cto Graph creation (node and edges weights share dimensionality)
	LLBGNGGraph (const unsigned int&);
	/// std dto
	~LLBGNGGraph(){}
	// set time window constants
	void setTimeWindows (unsigned int, unsigned int, unsigned int);
	// calculate inherited variables for a node to be inserted between two nodes
	void calculateInheritedParams (const unsigned int, const unsigned int, const unsigned int);
	// calculate long term and short term error for some node
	void updateAvgError (const unsigned int, T last_error);
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
	// set insertion tolerance
	void setInsertionTolerance (T&);
	// set initial learning rate constants
	void setLearningRates (T&, T&, T&);
	// set deletion threshold
	void setDeletionThreshold (T&);
	// get deletion_threshold
	T getDeletionThreshold () const;
	// set minimal node age constant
	void setMinimalNodeAge (T&);
	// set minimal node age constant
	T getMinimalNodeAge () const;
	// set maximal edge age constant
	void setMaximalEdgeAge (unsigned int&);
	// get maximal edge age constant
	unsigned int getMaximalEdgeAge () const;
	// set stabilization constant
	void setStabilization (T&);
	// get stabilization constant
	T getStabilization () const;
	// decrease age of some node
	void decreaseNodeAge (const unsigned int);
	// update activations counter for some node
	void increaseActivationsCounter (const unsigned int);
private:
	// returns a pointer to a edge of a type that is currently used by the graph
	virtual LLBGNGNode<T,S>* newNode(void);	
	/// adaptation threshold constant
	T adaptation_threshold;
	/// insertion tolerance constant
	T insertion_tolerance;
	/// initial winner learning rate constant
	T winner_learning_rate;
	/// initial winner-neighbors learning rate constant
	T neighbors_learning_rate;
	/// initial learning rate constant of the insertion threshold
	T insertion_learning_rate;
	/// minimal node age constant
	T minimal_node_age;
	/// maximal edge age constant
	unsigned int maximal_edge_age;
	/// stabilization constant
	T stabilization;
	/// deletion threshold
	T deletion_threshold;
	/// short term window constant
	unsigned int shortterm_window;
	/// long term window constant
	unsigned int longterm_window;
	/// age time window constant
	unsigned int age_time_window;
};

/** \brief cto Graph creation (node and edges weights share dimensionality)
    \param dim dimensionality
 */
template<typename T, typename S>
LLBGNGGraph<T,S>::LLBGNGGraph (const unsigned int &dim) :
	Base_Graph<T,S>(dim),
	UGraph<T,S>(dim),
	TGraph<T,S>(dim),
	GNGModulGraph<T,S>(dim),
	adaptation_threshold (0),
	insertion_tolerance (0),
	winner_learning_rate (0),
	neighbors_learning_rate (0),
	insertion_learning_rate (0),
	minimal_node_age (0),
	maximal_edge_age (0),
	stabilization (0),
	deletion_threshold (0),
	shortterm_window (1),
	longterm_window (1),
	age_time_window (1)
{
}


/** \brief overriden function from \p Base_Graph
    to create a node of type \p LLBGNGNode */
template<typename T,typename S>
LLBGNGNode<T,S>* LLBGNGGraph<T,S>::newNode(void)
{
	LLBGNGNode<T,S>* n = new LLBGNGNode<T,S>;
	if (this->high_limits.size() != 0 && this->low_limits.size() != 0)
		n->inherited_error = euclidean<T,S> (this->getHighLimits(), this->getLowLimits());
	else
		n->inherited_error = this->high_limit - this->low_limit;		
	//default inherited errors (used at the beginning of the algorithm)
	n->shortterm_avgerror = n->inherited_error;
	n->longterm_avgerror = n->inherited_error;
	n->errors.push (n->inherited_error);
	return n; 
}

/** \brief set window constants. This function should only
    be called once when creating the graph
    \param shortterm short term time window constant
    \param longterm long term time window constant
    \param age age time window constant */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setTimeWindows (unsigned int shortterm, unsigned int longterm, unsigned int age)
{
	assert (shortterm <= longterm);
	shortterm_window = shortterm;
	longterm_window = longterm;
	age_time_window = age;
}

/** \brief calculate inherited variables for a node to be inserted between two nodes
 *  \param index node index
 *  \param first_index 1st node index
 *  \param snd_index 2nd node index
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::calculateInheritedParams (const unsigned int index, const unsigned int first_index, const unsigned int snd_index)
{
	LLBGNGNode<T,S>* node = static_cast<LLBGNGNode<T,S>* >(this->_nodes[index]);
	LLBGNGNode<T,S>* first_node = static_cast<LLBGNGNode<T,S>* >(this->_nodes[first_index]);
	LLBGNGNode<T,S>* snd_node = static_cast<LLBGNGNode<T,S>* >(this->_nodes[snd_index]);
	
	node->weight = (first_node->weight + snd_node->weight) / (T)2;
	node->longterm_avgerror = (first_node->longterm_avgerror + snd_node->longterm_avgerror) / (T)2;
	node->shortterm_avgerror = (first_node->shortterm_avgerror + snd_node->shortterm_avgerror) / (T)2;
	node->inherited_error = (first_node->inherited_error + snd_node->inherited_error) / (T)2;
	node->insertion_threshold = (first_node->insertion_threshold + snd_node->insertion_threshold) / (T)2;

	//check if insertion is successful
	if (first_node->longterm_avgerror >= first_node->inherited_error * (1 - insertion_tolerance)) {
		first_node->insertion_threshold += insertion_learning_rate * (first_node->longterm_avgerror - first_node->insertion_threshold * (1 - insertion_tolerance));
	}
	
	if (snd_node->longterm_avgerror >= snd_node->inherited_error * (1 - insertion_tolerance)) {
		snd_node->insertion_threshold += insertion_learning_rate * (snd_node->longterm_avgerror - snd_node->insertion_threshold * (1 - insertion_tolerance));
	}
	if(node->longterm_avgerror >= node->inherited_error * (1 - insertion_tolerance)) {
		node->insertion_threshold += insertion_learning_rate * (node->longterm_avgerror - node->insertion_threshold * (1 - insertion_tolerance));
		
	}
	first_node->inherited_error = first_node->longterm_avgerror;	
	snd_node->inherited_error = snd_node->longterm_avgerror;
	node->inherited_error = node->longterm_avgerror;

	assert (node->errors.size() == 1);
	node->errors.front() = node->longterm_avgerror;

	
}

/** \brief calculate long term and short term error for a given node
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph<T,S>::updateAvgError (const unsigned int index, T last_error)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->updateAvgError(last_error, shortterm_window, longterm_window);

}

/** \brief calculate learning quality for a node
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph<T,S>::calculateLearningQuality (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->calculateLearningQuality();
	
}

/** \brief calculate insertion quality for a node
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph<T,S>::calculateInsertionQuality (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->calculateInsertionQuality(insertion_tolerance);
	
}

/** \brief calculate insertion criterion for a node
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph<T,S>::calculateInsertionCriterion (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->calculateInsertionCriterion();
	
}

/** \brief set initial learning rate constants
 *  \param winner initial winner learning rate constant
 *  \param neighbors initial winner-neighbors learning rate constant
 *  \param insertion_threshold learning rate constant
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setLearningRates (T& winner, T& neighbors, T& insertion_threshold)
{
	winner_learning_rate = winner;
	neighbors_learning_rate = neighbors;
	insertion_learning_rate = insertion_threshold;
}


/** \brief update learning rate for a winner node
 *  \param index node index
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::updateWinnerLearningRate (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->updateLearningRate(adaptation_threshold, winner_learning_rate);	
}

/** \brief update learning rate for a winner-neighboring node
 *  \param index node index
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::updateNeighborLearningRate (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->updateLearningRate(adaptation_threshold, neighbors_learning_rate);	
}


/** \brief set adaptation threshold
    \param threshold given adaptation threshold */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setAdaptationThreshold (T& threshold)
{
	adaptation_threshold = threshold;	
}

/** \brief set insertion tolerance
    \param tolerance insertion tolerances */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setInsertionTolerance (T& tolerance)
{
	insertion_tolerance = tolerance;
}

/** \brief set \p deletion_threshold
    \param threshold given deletion threshold */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setDeletionThreshold (T& threshold)
{
	deletion_threshold = threshold;
}

/** \brief get \p deletion_threshold
 */
template<typename T, typename S>
T LLBGNGGraph<T,S>::getDeletionThreshold () const
{
	return deletion_threshold;
}


/* \brief set \p minimal_node_age constant
 * \param age minimal node age constant
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setMinimalNodeAge (T& age)
{
	minimal_node_age = age;
}

/* \brief get \p minimal_node_age constant
 */
template<typename T, typename S>
T LLBGNGGraph<T,S>::getMinimalNodeAge () const
{
	return minimal_node_age;
}

/* \brief set maximal edge age constant
 * \param age maximal edge age constant
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setMaximalEdgeAge (unsigned int& age)
{
	maximal_edge_age = age;
}

/* \brief get \p maximal_edge_age constant
 */
template<typename T, typename S>
unsigned int LLBGNGGraph<T,S>::getMaximalEdgeAge () const
{
	return maximal_edge_age;
}


/* set \p stabilization constant
 * \param s stabilization constant
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::setStabilization (T& s)
{
	stabilization = s;
}

/* get \p stabilization constant
 */
template<typename T, typename S>
T LLBGNGGraph<T,S>::getStabilization () const
{
	return stabilization;
}

/* \brief decrease age of a node
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::decreaseNodeAge (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->decreaseAge (age_time_window);
}

/** \brief update activations counter for some node
 */
template<typename T, typename S>
void LLBGNGGraph<T,S>::increaseActivationsCounter (const unsigned int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->increaseActivationsCounter ();
}


} // namespace neuralgas

#endif
