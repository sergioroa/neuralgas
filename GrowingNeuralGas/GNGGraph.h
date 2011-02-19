/** 
* \class GNGGraph
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/


#ifndef GNGGRAPH_H
#define GNGGRAPH_H

#include <iostream>
#include "GNGModulGraph.h"

namespace neuralgas {

/** \brief The derived node for a GNGGraph
*
*   The node permits having an error implementing the algorithm described in 
*   "Growing Neural Gas for Temporal Clustering" by Isaac J. Sledge and James M. Keller.
*   The error represents the distortion that was measured between the node weight 
*   and the actual data vector by an user defined metric.
*   There are four functions, decCounter, incCounter, setCounter and getCounter
*   for manipulating the value. 
*   Furthermore it has a context vector as suggested in the above paper to reflect
*   the history and the time series evolution of the node.
*   \param counter the belief in the node
*   \param vector reflecting the time series evolution
*
*/


template<typename T,typename S> struct GNGNode : Base_Node<T,S>
{
	// sets the value of the error
	void setError(const T& newError){error=newError;}
	// returns the error
	T getError() const {return error;}
	//accumalted error
	T error;
};

/** \brief GNGGraph provides the graph structure for the algorithm proposed in
*   "Growing Neural Gas for Temporal Clustering" by Isaac J. Sledge and James M. Keller.
*
* The class is a undirected time graph permiting the access to the GNGnode 
* via potential dangerous downcastings.
* The following paragraph is important for derived classes. 
* It is possible to use a user defined type of node for the graph.
* In this case that newly defined node has to be derived from the struct Base_Node<T,S>
* and the virtual function newNode() has to be overloaded in the derived graph class such that 
* it returns a pointer to that newly defined but derived node data type. 
* The same holds for the edges, where the edge has to be derived from the timed edge TEdge
* in the TGraph class.
*
*/

template<typename T,typename S> class GNGGraph : public GNGModulGraph<T,S>
{
public:
	//cto creating a graph with the same dimension for node and edge weight vectors
	GNGGraph(const int& dim) :  Base_Graph<T,S>(dim),UGraph<T,S>(dim),TGraph<T,S>(dim),GNGModulGraph<T,S>(dim){}
	//cto creating a graph with the different dimension for node and edge weight vectors
	GNGGraph(const int& dimNode,const int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge),UGraph<T,S>(dimNode,dimEdge),TGraph<T,S>(dimNode,dimEdge){}
	// std dto
	~GNGGraph(){}
	// sets a new counter value for the given node
	void inline setError(const int&, const T&);
	// gets the counter value for the given node
	T inline getError(const int&) const;
	
private:
	// returns a pointer to a node of a type that is currently used by the graph
	virtual GNGNode<T,S>*             newNode(void);
};


/** \brief returns a pointer to a node of a type that is currently used by the graph
*
* This virtual function allows the use of a (sub)class specific node type.
* If subclasses want to use another type of node, then this node has to be derived from
* the struct GNGNode<T,S> and the function newNode() has to be reimplemented in the subclass
* such that the reimplemented function returns a pointer to that new defined and derived 
* node type.
* This function is called by the function addNode() and inheritance rules guarantee
* that the function of the subclass and not of the superclass is called within the function
* addNode(),resulting in the use of the used defined node type as node for the graph structure.
*/

template<typename T,typename S> GNGNode<T,S>* GNGGraph<T,S>::newNode(void)
{
  GNGNode<T,S>* n = new GNGNode<T,S>;
  n->error=0;
    
  return n; 
}

/** \brief sets a new error value for the given node
*
* \param index is the index of the node that shall get a new error value
* \param newError is the value to be set
*/
template<typename T,typename S> void inline GNGGraph<T,S>::setError(const int& index, const T& newError)
{
 (static_cast< GNGNode<T,S>* > (this->_nodes[index]))->setError(newError);
}
/** \brief gets the error value for the given node
*
* \param index is the index of the node wherefor the error is desired
*/
 
template<typename T,typename S> T inline GNGGraph<T,S>::getError(const int& index) const
{ 
  return (static_cast< GNGNode<T,S>* > (this->_nodes[index]))->getError();
}

} // namespace neuralgas

#endif
