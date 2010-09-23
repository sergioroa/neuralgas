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
* Teh same holds for the edges.
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


//returns whether the two nodes are connected by an edge
template<typename T,typename S> inline bool UGraph<T,S>::areConnected(const int& x,const int& y)const
{
  return ( (*_nodes[x])[y] != NULL );
} 

template<typename T,typename S> void UGraph<T,S>::addEdge(const int& x,const int& y)
{
  if ( (*_nodes[x])[y] == NULL ) //self edges are allowed
  {
    Base_Edge<S,T>* new_edge =  newEdge();
    new_edge->in             = _nodes[x];
    new_edge->out            = _nodes[y];
    new_edge->weight.resize(_dimEdge);
    (*_nodes[x]).edges[y]    =  new_edge;
    
    (*_nodes[y]).edges[x]    =  new_edge;
    (*_nodes[x]).num_connections++;
    (*_nodes[y]).num_connections++;  
  }      
  
}  

template<typename T,typename S> void UGraph<T,S>::rmEdge(const int& x,const int& y)
{
  if ( (*_nodes[x])[y] != NULL ) //self edges are allowed
  {
    delete (*_nodes[x])[y];
    (*_nodes[x]).edges[y]     =  NULL;
    (*_nodes[y]).edges[x]     =  NULL;
    (*_nodes[x]).num_connections--;
    (*_nodes[y]).num_connections--;  
  }  
}

#endif
