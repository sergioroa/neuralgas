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
* \file GNGModulGraph.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  Copyright(c) 2011 Sergio Roa
*  \version 1.0
*  \date    2010
*/

#ifndef GNGMODULGRAPH_H
#define GNGMODULGRAPH_H

#include <Graphs/UGraph.h>
#include <Graphs/TGraph.h>
#include <limits>

namespace neuralgas {

template<typename T,typename S> class ErrorTesting;

/** \class GNGModulGraph
 *  \brief provides the basic graph structure for MergeGrowingNeuralGas
 *   algorithms and their extensions.
 *
 * The class is a undirected time graph.
 * The following paragraph is important for derived classes. 
 * It is possible to use a user defined type of node for the graph.
 * In this case that newly defined node has to be derived from the struct Base_Node<T,S>
 * and the virtual function newNode() has to be overloaded in the derived graph class such that 
 * it returns a pointer to that newly defined but derived node data type. 
 * The same holds for the edges, but in the edge has to be derived from the struct TEdge<S,T>
 * contained in the class TGraph.
 *
 * \param low_limit min value for the random initializiation of the context vector 
 * \param high_limit max value for the random initializiation of the context vector 
 */
template<typename T,typename S> class GNGModulGraph : public virtual UGraph<T,S>, public virtual TGraph<T,S>
{

public:
	//cto creating a graph with the same dimension for node and edge weight vectors
	GNGModulGraph(const unsigned int& dim) : Base_Graph<T,S>(dim),UGraph<T,S>(dim),TGraph<T,S>(dim){}
	//cto creating a graph with the different dimension for node and edge weight vectors
	GNGModulGraph(const unsigned int& dimNode,const unsigned int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge),UGraph<T,S>(dimNode,dimEdge),TGraph<T,S>(dimNode,dimEdge){}
	/// dummy cto
	GNGModulGraph () {};
	/// std dto
	virtual ~GNGModulGraph(){}
	// adds an edge between the nodes given by their indeces by calling corresponding upper class method of UGraph
	virtual void inline  addEdge(const unsigned int& x,const unsigned int& y){if (x!=y) UGraph<T,S>::addEdge(x,y); }
	// removes an edge between the nodes given by their indeces if there exists one by calling corresponding upper class method of UGraph    
	virtual void inline  rmEdge(const unsigned int& x,const unsigned int& y){if (x!=y) UGraph<T,S>::rmEdge(x,y);}

	/// sets the minimal limit values
	void setLowLimits(Vector<T> low);
	/// sets the maximal limit values
	void setHighLimits(Vector<T> high);
	/// returns the minimal limit values
	Vector<T> getLowLimits() const{return low_limits;}
	/// returns the maximal limit values
	Vector<T> getHighLimits() const{return high_limits;}
	/// returns minimal limit value
	T getLowLimit() const{return low_limit;}
	/// returns maximal limit value
	T getHighLimit() const{return high_limit;}
	// // class dependent distance function that is used within the winner function
	// //virtual T getDistance(const Vector<T>&,const unsigned int&)=0;
	virtual T getDistance(const Vector<T>&,const unsigned int&) const;
	// func determines for the current time step / data item the two most similar nodes
	virtual T getWinner( unsigned int&, unsigned int&, const Vector<T>&) const;
   

protected:

	/// low limit min value for the random initializiation of the context vector
	T low_limit;
	/// high_limit max value for the random initializiation of the context vector
	T high_limit;
	/// low limit min value for the random initializiation of the context vector for each dimension
	Vector<T> low_limits;
	/// high_limit max value for the random initializiation of the context vector for each dimension
	Vector<T> high_limits;
private:
	/// ErrorTesting is defined as friend in order to not having duplicate anything
	friend class ErrorTesting<T,S>;

};

/** \brief func determines for the current data item the two most similar nodes
*
*   The metric is used (either the user defined or the preset one to determine the distance
*   from the current data to the current nodes in the network.
*   The two closest nodes are stored in the parameter given by reference.
*  
*   NOTE: It is intented to "outsource" this function and to let the user define the way
*         of determining the winner and refer to it via a function ptr.
*  
*   \param first_winner before func call an arbitrary value, after func call the closest node
*   \param second_winner before func call an arbitrary value, after func call the second closest node
*   \param time is the current time step reflecting the current data to be processed
*/
template<typename T,typename S>
T GNGModulGraph<T,S>::getWinner(unsigned int& first_winner, unsigned int& second_winner, const Vector<T>& item) const
{
	T distance;
	// best_distance set to "infinity"
	T best_distance = std::numeric_limits<T>::max();

	assert (this->size() >= 2);

	first_winner = 1;
	second_winner = 0;
 
	for (unsigned int j = 0; j < this->size(); j++)
	{
		distance = getDistance(item,j);
		if (distance < best_distance)
		{
			second_winner              =       first_winner;
			first_winner               =       j;
			best_distance              =       distance;
		}       
	}
	return best_distance;
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
inline T GNGModulGraph<T,S>::getDistance(const Vector<T>& item, const unsigned int& node_index) const
{
	return metric( item, this->_nodes[node_index]->weight);
}


template<typename T, typename S>
void GNGModulGraph<T,S>::setHighLimits (Vector<T> high)
{
	high_limits = high;
	high_limit = high_limits[0];
	for (unsigned int i=1; i<high_limits.size(); i++)
		if (high_limits[i] > high_limit)
			high_limit = high_limits[i];
	
}

template<typename T, typename S>
void GNGModulGraph<T,S>::setLowLimits (Vector<T> low)
{
	low_limits = low;
	low_limit = low_limits[0];
	for (unsigned int i=1; i<low_limits.size(); i++)
		if (low_limits[i] < low_limit)
			low_limit = low_limits[i];
	
}


} // namespace neuralgas

#endif
