/** 
* \class GNGModulGraph
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef GNGMODULGRAPH_H
#define GNGMODULGRAPH_H

#include <Graphs/UGraph.h>
#include <Graphs/TGraph.h>


/** \brief GNGModulGraph provides the basic graph structure for MergeGrowingNeuralGas
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
* \param maxRandomValue max value for the random initializiation of the context vector 
*/

template<typename T,typename S> class GNGModulGraph : public virtual UGraph<T,S>, public virtual TGraph<T,S>
{

 public:
   //cto creating a graph with the same dimension for node and edge weight vectors
   GNGModulGraph(const int& dim) : Base_Graph<T,S>(dim),UGraph<T,S>(dim),TGraph<T,S>(dim){}
   //cto creating a graph with the different dimension for node and edge weight vectors
   GNGModulGraph(const int& dimNode,const int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge),UGraph<T,S>(dimNode,dimEdge),TGraph<T,S>(dimNode,dimEdge){}
   // std dto
   ~GNGModulGraph(){}
   // function sets the max random value for the random initializiation of the context vector 
   void       setMaxRandomValue(const int& value){maxRandomValue = value;}
   // adds an edge between the nodes given by their indeces by calling corresponding upper class method of UGraph
   virtual void inline  addEdge(const int& x,const int& y){if (x!=y) UGraph<T,S>::addEdge(x,y); }
   // removes an edge between the nodes given by their indeces if there exists one by calling corresponding upper class method of UGraph    
   virtual void inline  rmEdge(const int& x,const int& y){if (x!=y) UGraph<T,S>::rmEdge(x,y);}

 protected:
   // returns the maximal random value
   int        getMaxRandomValue() const{return maxRandomValue;}

 private:
   // maxRandomValue max value for the random initializiation of the context vector
   int                                maxRandomValue;
};

#endif
