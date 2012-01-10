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
* \file MGNGGraph.h
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/


#ifndef MGNGGRAPH_H
#define MGNGGRAPH_H

#include <iostream>
#include <GrowingNeuralGas/GNGModulGraph.h>

namespace neuralgas {

template<typename T, typename S> class MGNGAlgorithm;

/** \class MGNGNode
 *  \brief The derived node for a MGNGGraph
 *
 *   The node permits having a counter implementing the algorithm described in 
 *   "Incremental Unsupervised Time Series Analysis using Merge Growing Neural Gas"
 *   by Andreakis,Hoyningen-Huene and Beets.
 *   The counter represents the amount of winning in this algorithm. The higher the value
 *   the more reliable and important is the node.
 *   There are four functions, decCounter, incCounter, setCounter and getCounter
 *   for manipulating the value. 
 *   Furthermore it has a context vector as suggested in the above mentioned paper to reflect
 *   the history and the time series evolution of the node.
 *   \param counter the belief in the node
 *   \param context vector reflecting the time series evolution
 *
 */
template<typename T,typename S> struct MGNGNode : Base_Node<T,S>
{
 virtual ~MGNGNode () {}
	
 // sets the value of the counter
 void setCounter(const float& newCounter){counter=newCounter;}
 // returns the counter
 float getCounter() const {return counter;}
 // dec the counter by one
 void  decCounter() {--counter;}
 // inc the counter by one
 void  incCounter() {++counter;}
 // vector representing the time series evolution of that node
 Vector<T> context;
 //belief in the node
 float counter;
 
 int getBirthday() const {return birthday;}
 void setBirthday(const int& date){birthday=date;}
 int birthday;
 
 float temperature;
 void setTemp(const float& temp){temperature=temp;}
 float getTemp(){return temperature;} 
};

/** \class MGNGGraph
 * \brief provides the graph structure for the algorithm proposed in
 *   "Incremental Unsupervised Time Series Analysis using Merge Growing Neural Gas"
 *   by Andreakis,Hoyningen-Huene and Beets.
 *
 * The class is a undirected time graph permiting the access to the MGNGnode 
 * via potential dangerous downcastings.
 * The following paragraph is important for derived classes. 
 * It is possible to use a user defined type of node for the graph.
 * In this case that newly defined node has to be derived from the struct Base_Node<T,S>
 * and the virtual function newNode() has to be overloaded in the derived graph class such that 
 * it returns a pointer to that newly defined but derived node data type. 
 * The same holds for the edges,where the edge has to be derived from the timed edge TEdge
 * in the TGraph class.
 *
 */
template<typename T,typename S> class MGNGGraph : public GNGModulGraph<T,S>
{
 public:
   //cto creating a graph with the same dimension for node and edge weight vectors
   MGNGGraph(const int& dim) :  Base_Graph<T,S>(dim),UGraph<T,S>(dim),TGraph<T,S>(dim),GNGModulGraph<T,S>(dim){}
   //cto creating a graph with the different dimension for node and edge weight vectors
   MGNGGraph(const int& dimNode,const int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge),UGraph<T,S>(dimNode,dimEdge),TGraph<T,S>(dimNode,dimEdge){}
   // std dto
   ~MGNGGraph();
   // sets a new counter value for the given node
   void inline setCounter(const int&, const float&);
   // decreases the counter by one
   void inline decCounter(const int&);
   // increases the counter by one
   void inline incCounter(const int&);
   // gets the counter value for the given node
   float inline getCounter(const int&) const;
   // returns a reference to the given context vector
   Vector<T>& context(const int&) const;
   
   int getBirthday(const int& index) const 
   {return (dynamic_cast< MGNGNode<T,S>* > (this->_nodes[index]))->getBirthday();}
   void setBirthday(const int& index,const int& date)
   {(dynamic_cast< MGNGNode<T,S>* > (this->_nodes[index]))->setBirthday(date);}


   float getTemp(const int& index) const 
   {return (dynamic_cast< MGNGNode<T,S>* > (this->_nodes[index]))->getTemp();}
   void setTemp(const int& index,const float& temp)
   {(dynamic_cast< MGNGNode<T,S>* > (this->_nodes[index]))->setTemp(temp);}
   // set algorithm to use its parameters for, e.g., distance calculation
   void setAlgorithm (MGNGAlgorithm<T,S>* alg) { algorithm = alg; }
   // Algorithmic dependent distance function
   T getDistance(const Vector<T>&,const unsigned int&) const;

protected:
   // pointer to algorithm
   MGNGAlgorithm<T,S>* algorithm;

    
 private:
   // returns a pointer to a edge of a type that is currently used by the graph
   virtual MGNGNode<T,S>*             newNode(void);
};


/** \brief dto Graph deletion
 */
template<typename T, typename S>
MGNGGraph<T,S>::~MGNGGraph ()
{
	for (unsigned int i=0; i<this->size(); i++)
		for (unsigned int j=i+1; j<this->_nodes[i]->edges.size(); j++)
			if (this->_nodes[i]->edges[j] != NULL)
				delete this->_nodes[i]->edges[j];
	
	for (unsigned int i=0; i<this->size(); i++)
		this->_nodes[i]->edges.clear();
}

/** \brief returns a pointer to a node of a type that is currently used by the graph
*
* This virtual function allows the use of a (sub)class specific node type.
* If subclasses want to use another type of node, then this node has to be derived from
* the struct MGNGNode<T,S> and the function newNode() has to be reimplemented in the subclass
* such that the reimplemented function returns a pointer to that new defined and derived 
* node type.
* This function is called by the function addNode() and inheritance rules guarantee
* that the function of the subclass and not of the superclass is called within the function
* addNode(),resulting in the use of the used defined node type as node for the graph structure.
*/

template<typename T,typename S> MGNGNode<T,S>* MGNGGraph<T,S>::newNode(void)
{
  MGNGNode<T,S>* n = new MGNGNode<T,S>;
  n->temperature=0.0;
  n->context.resize(this->_dimNode);          // sets dimension of the weight vector
  n->counter=0.0;
  for(unsigned int j = 0; j < this->_dimNode; j++)          
    //n->context[j] = (T)(rand() % this->getMaxRandomValue() ); //sets the value of the weights to random values
  n->context[j] = (T) ((rand() / (static_cast<T>(RAND_MAX) + 1.0)) * (this->high_limits[j] - this->low_limits[j]) + this->low_limits[j] );
  
  return n; 
}
/** \brief sets a new counter value for the given node
*
* \param index is the index of the node that shall get a new counter value
* \param newCounter is the value to be set
*/
template<typename T,typename S> void inline MGNGGraph<T,S>::setCounter(const int& index, const float& newCounter)
{
 (static_cast< MGNGNode<T,S>* > (this->_nodes[index]))->setCounter(newCounter);
}
/** \brief gets the counter value for the given node
*
* \param index is the index of the node wherefor the counter is desired
*/
 
template<typename T,typename S> float inline MGNGGraph<T,S>::getCounter(const int& index) const
{ 
return (static_cast< MGNGNode<T,S>* > (this->_nodes[index]))->getCounter();
}
/** \brief decreases the counter by one
*
* \param index is the index of the node that shall be decreased
*/
 
template<typename T,typename S> void inline MGNGGraph<T,S>::decCounter(const int& index)
{
(static_cast< MGNGNode<T,S>* > (this->_nodes[index]))->decCounter();
}
/** \brief increases the counter by one
*
* \param index is the index of the node that shall be increased
*/
 
template<typename T,typename S> void inline MGNGGraph<T,S>::incCounter(const int& index)
{ 
(static_cast< MGNGNode<T,S>* > (this->_nodes[index]))->incCounter();
}
/** \brief returns a reference to the given context vector
*
* \param index is the index of the node wherefore the context is desired
*/
  
template<typename T,typename S> Vector<T>& MGNGGraph<T,S>::context(const int& index) const
{
return (static_cast< MGNGNode<T,S>* > (this->_nodes[index]))->context;
}

/** \brief Algorithmic dependent distance function
*
*   This function returns the distance of the given datum and the given node. 
*   The distance is a algorithmic dependent function that is either
*   just the setted metric or a combination thereof.
*   Currently dist  = (1-a)*metric(x_t,w_j)^2+a*metric(C,c_j)^2 where
*   x_t is the data vector and w_j the node vector, C the global and c_j the local
*   context vector where the latter one belongs to w_j.
*   
*   \param item datum
*   \param node_index is the node where to the distance shall be determined
*/
template<typename T,typename S> T MGNGGraph<T,S>::getDistance(const Vector<T>& item, const unsigned int& node_index) const
{
    // dist  = (1-a)*metric(x_t,w_j)^2+a*metric(C,c_j)^2
    T distance = (1 - algorithm->params[0]);
    distance*=pow(metric( item, this->_nodes[node_index]->weight),2);
    distance+=algorithm->params[0]*pow(metric(algorithm->globalContextV,context(node_index)) ,2);

    return distance;
}


} // namespace neuralgas

#endif
