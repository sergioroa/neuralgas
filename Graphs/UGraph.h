/** 
* \class UGraph
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/


#ifndef UGRAPH_H
#define UGRAPH_H

#include "Base_Graph.h"


/** \brief UGraph provides a structure for an undirected graph with n-dim weighted nodes and m-dim weighted edges
*
* The Base_Graph class is extended such that edge handling like adding and removing 
* is implemented such that the underlying data structure represents an undirected graph.
* The first template parameter is the type of the node weigth vectors and the second template
* parameter is the type of the edge weigth vectors.
* The following paragraph is important for derived classes. 
* It is possible to use a user defined type of node for the graph.
* In this case that newly defined node has to be derived from the struct Base_Node<T,S>
* and the virtual function newNode() has to be overloaded in the derived graph class such that 
* it returns a pointer to that newly defined but derived node data type. 
* The same holds for the edges.
*
* 
*/

template<typename T,typename S> class UGraph : virtual public Base_Graph<T,S>
{
  public:
    //cto creating a graph with the same dimension for node and edge weight vectors
                        UGraph(const int& dim) : Base_Graph<T,S>(dim){}
    //cto creating a graph with the different dimension for node and edge weight vectors
                        UGraph(const int& dimNode,const int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge){}  
    //std dto
                        ~UGraph(){}
    //returns whether the two nodes are connected by an edge
    inline bool         areConnected(const int&,const int&)const;
    // adds an edge between the nodes given by their indeces      
    void                addEdge(const int&,const int&);
    // removes an edge between the nodes given by their indeces if there exists one     
    void                rmEdge(const int&,const int&);

};  

/** \brief returns whether the two nodes are connected by an edge
*/
template<typename T,typename S> inline bool UGraph<T,S>::areConnected(const int& x,const int& y)const
{
  return ( this->_nodes[x]->edges[y] != NULL );
} 


/** \brief Adds an edge between the nodes given by their indeces 
*
*   If there is no edge then it is added and the corresponding number of connections 
*   in the ingoing and outgoing node are updated.
*   The edge from x to y is the same as the edge from y to x.
* 
*   \param x outgoing node, plays no role in the undirected graph
*   \param y ingoing node, plays no role in the undirected graph
*/

template<typename T,typename S> void UGraph<T,S>::addEdge(const int& x,const int& y)
{
  if ( this->_nodes[x]->edges[y] == NULL ) //self edges are allowed
  {
    Base_Edge<S,T>* new_edge = this->newEdge();
    new_edge->in             = this->_nodes[x]; //sets the ingoing node
    new_edge->out            = this->_nodes[y]; //sets the outgoing node
    new_edge->weight.resize(this->_dimEdge);    //resizes the dim of the edge weight vec
    
    this->_nodes[x]->edges[y]    =  new_edge;   //shares the same edge since its undirected
    this->_nodes[y]->edges[x]    =  new_edge;
    
    this->_nodes[x]->num_connections++;      // number of connections are increased
    this->_nodes[y]->num_connections++;  
  }      
  
}  

/** \brief Removes an edge between the nodes given by their indeces if there exists one 
*
*   If there is an edge then it is removed and the corresponding number of connections 
*   in the ingoing and outgoing node are updated.
*   The edge from x to y is the same as the edge from y to x.
* 
*   \param x outgoing node, plays no role in the undirected graph
*   \param y ingoing node, plays no role in the undirected graph
*/
template<typename T,typename S> void UGraph<T,S>::rmEdge(const int& x,const int& y)
{

  if ( this->_nodes[x]->edges[y] != NULL ) //self edges are allowed
  {
    delete this->_nodes[x]->edges[y];
    this->_nodes[x]->edges[y]     =  NULL;
    this->_nodes[y]->edges[x]     =  NULL;
    this->_nodes[x]->num_connections--;
    this->_nodes[y]->num_connections--;  
  }  
}

#endif
