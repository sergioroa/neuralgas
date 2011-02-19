/** 
* \class LLBGNGGraph
* \author Sergio Roa
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef LLBGNGGRAPH_H
#define LLBGNGGRAPH_H

#include "GNGModulGraph.h"

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
	/// errors queue for calculating different parameters. Its size is
	/// set by the maximum allowable long term error window
	std::queue<T> errors;
	// calculate \p learning_quality measure
	void calculateLearningQuality ();
	/// quality measure for learning
	T learning_quality;
	// calculate \p shortterm_avgerror and \p longterm_avgerror
	void updateAvgError ();
	/// short-term error counter
	T shortterm_avgerror;
	/// long-term error counter
	T longterm_avgerror;
	/// An insertion is only allowed if the long-term error exceeds this threshold
	T insertion_threshold;
	/// short term window constant
	unsigned int shortterm_window;
	/// long term window constant
	unsigned int longterm_window;
	/// \p longterm_window / \p shortterm_window ratio
	unsigned int errorwindows_ratio;
	/// age of the node
	unsigned int age;
	// update \p learning_rate
	void updateLearningRate (T&);
	/// learning rate of the node
	T learning_rate;
};

/// \brief default \p LLBGNGNode cto
template<typename T, typename S>
void LLBGNGNode::LLBGNGNode () :
	shortterm_window (0),
	shortterm_avgerror (0),
	longterm_avgerror (0)
{
	//this value should be modified to a proper value
	erroritems_counter = 0;
	errorwindows_ratio = 1;
}

/// \brief calculate \p learning_quality measure
template<typename T, typename S>
void LLBGNGNode::calculateLearningQuality ()
{
	learning_quality = (shortterm_avgerror + 1) / (longterm_avgerror + 1);
}

/** \brief update \p learning_quality measure
 *  \param adaptation_threshold cut-off value of learning
 */
template<typename T, typename S>
void LLBGNGNode::updateLearningRate (T& adaptation_threshold)
{
	learning_quality = (shortterm_avgerror + T(1)) / (longterm_avgerror + T(1));

	T rate_decision_boundary;

	rate_change_factor = learning_quality / (T(1) + adaptation_threshold) + age -1;

	if (rate_decision_boundary < 0)
		learning_rate = T(0);
	else if (rate_decision_boundary <= 1)
		learning_rate = rate_change_factor * learning_rate;
}


/// calculate \p shortterm_avgerror and \p longterm_avgerror
template<typename T, typename S>
void LLBGNGNode::updateAvgError (T& last_error)
{
	errors.push (last_error);
	float front_error = errors.front();
	unsigned int shortterm_window_scaled = shortterm_window;
	unsigned int longterm_window_scaled = longterm_window;
	if (errors.size() < longterm_window)
	{
		shortterm_window_scaled = errors.size() / errorwindows_ratio;
		longterm_window_scaled = errors.size();
	}
	
	shortterm_avgerror += ( (last_error - front_error) / shortterm_window_scaled);
	longterm_avgerror += ( (last_error - front_error) / longterm_window_scaled);

	if (errors.size() > longterm_window)
		errors.pop ();
}


/** \brief LLBGNGGraph provides some additional learning strategies proposed in
 *   the algorithm Life-long learning cell structures -- continuosly learning
 *   without catastrophic interference by Fred H. Hamker.
 *   The algorithm tries to avoid bias/variance issues in the original GNG.
*/
template<typename T, typename S>
class LLBGNGGraph : public GNGModulGraph
{
public:
	/// cto Graph creation (node and edges weights share dimensionality)
	LLBGNGGraph (const int &dim) : Base_Graph<T,S>(dim),UGraph<T,S>(dim),TGraph<T,S>(dim),GNGModulGraph<T,S>(dim) {}
	/// std dto
	~GNGGraph(){}	
	// set long term window constant for some node
	void setErrorTimeWindows (const int, unsigned int);
	// set insertion threshold for some node
	void setInsertionThreshold (const int, T&);
	// set inherited error variables for some node
	void setInheritedParams (const int, T&);
	// calculate long term and short term error for some node
	void updateAvgError (const int);
	// calculate learning quality for some node
	void calculateLearningQuality (const int);
	// update learning rate for some node
	void updateLearningRate (const int);
	// set input adaptation threshold
	void setAdaptationThreshold (T&);

private:
	// returns a pointer to a edge of a type that is currently used by the graph
	virtual LLBGNGNode<T,S>* newNode(void);	
	/// adaptation threshold constant
	T adaptation_threshold;
	
};

/** \brief overriden function from \p Base_Graph
    to create a node of type \p LLBGNGNode */
template<typename T,typename S>
LLBGNGNode<T,S>* LLBGNGGraph<T,S>::newNode(void)
{
	LLBGNGNode<T,S>* n = new LLBGNGNode<T,S>;    
	return n; 
}

/** \brief set long term window constant for a node. This function should only
    be called for the first nodes created in the graph
    \param index node index
    \param shortterm_window short term time constant
    \param longterm_window long term time constant */
template<typename T, typename S>
void LLBGNGGraph::setErrorTimeWindows (const int index, unsigned int shortterm_window, unsigned int longterm_window)
{
	assert (shortterm_window <= longterm_window);
	LLBGNGNode *node = static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]);
	node->shortterm_window = shortterm_window;
	node->longterm_window = longterm_window;
	node->errorwindows_ratio = longterm_window / shortterm_window;
}

/** \brief set insertion threshold for a node
    \param index node index
    \param threshold insertion threshold constant */
template<typename T, typename S>
void LLBGNGGraph::setInsertionThreshold (const int index, T& threshold)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->insertion_threshold = threshold;
}

/** \brief set inherited error variables for a node
    \param index node index
    \param firstparent_index 1st parent node index
    \param sndparent_index 2nd parent node index */
template<typename T, typename S>
void LLBGNGGraph::setInheritedParams (const int index, const int firstparent_index, const int sndparent_index)
{
	LLBGNGNode* node = static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]);
	LLBGNGNode* firstparent_node = static_cast<LLBGNGNode<T,S>* > (this->_nodes[firstparent_index]);
	LLBGNGNode* sndparent_node = static_cast<LLBGNGNode<T,S>* > (this->_nodes[sndparent_index]);
	node->longterm_avgerror = (firstparent_node->longterm_avgerror + sndparent_node->longterm_avgerror) / (T)2;
	node->shortterm_avgerror = (firstparent_node->shortterm_avgerror + sndparent_node->shortterm_avgerror) / (T)2;
	node->shortterm_window = firstparent_node->shortterm_window;
	node->longterm_window = firstparent_node->longterm_window;
	node->errorwindows_ratio = longterm_window / shortterm_window;
	node->errors.push (longterm_avgerror);
}

/** \brief calculate short term error for a node given \p shortterm_error
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph::updateAvgError (const int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->updateAvgError();

}

/** \brief calculate learning quality for a node
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph::calculateLearningQuality (const int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->calculateLearningQuality();
	
}

/** \brief update learning rate for a node
    \param index node index */
template<typename T, typename S>
void LLBGNGGraph::updateLearningRate (const int index)
{
	(static_cast<LLBGNGNode<T,S>* > (this->_nodes[index]))->updateLearningRate(adaptation_threshold);
	
}


/** \brief set adaptation threshold
    \param threshold given adaptation threshold */
template<typename T, typename S>
void LLBGNGGraph::setAdaptationThreshold (T& threshold)
{
	adaptation_threshold = threshold;
	
}


} // namespace neuralgas
