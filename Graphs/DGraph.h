/** 
* \class DGraph
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/


#ifndef DGRAPH_H
#define DGRAPH_H

#include "Base_Graph.h"

namespace neuralgas {

//TODO doxygen funcs & mark the dangerous upcast within the code

/** \brief The derived node for directed graph
*
*   The node permits having different in and outgoing edges which is essential
*   for a directed graph.
*   Since this was already implied in the Base_Node but not the possibility
*   to reflect the distinction of in and outgoing edges within the variable 
*   num_connections containing the number of connections.Therefore two new variables
*   num_in_edges and num_out_edges are added representing the number of 
*   ingoing respectively outgoing edges.
*
*   \param num_in_edges number of ingoing edges
*   \param num_out_edges number of outgoing edges
*/

template < typename T , typename S > struct DNode : Base_Node<T,S>
{
 DNode()
 {num_in_edges=0; num_out_edges=0;}  
 //increases the number of ingoing edges
 void incInEdges(){++num_in_edges;} 
 //increases the number of outgoing edges
 void incOutEdges(){++num_out_edges;std::cout << num_out_edges<<std::endl;}
 //decreases the number of ingoing edges
 void decInEdges(){--num_in_edges;} 
 //decreases the number of outgoing edges
 void decOutEdges(){--num_out_edges;}
 //number of ingoing edges
 int num_in_edges;
 //number of outgoing edges
 int num_out_edges;   
};


/** \brief DGraph provides a structure for a directed graph with n-dim weighted nodes and m-dim weighted edges
*
* The Base_Graph class is extended such that edge handling like adding and removing 
* is implemented such that the underlying data structure represents an directed graph.
* The first template parameter is the type of the node weigth vectors and the second template
* parameter is the type of the edge weigth vectors.
* The following paragraph is important for derived classes. 
* It is possible to use a user defined type of node for the graph.
* In this case that newly defined node has to be derived from the struct Base_Node<T,S>
* and the virtual function newNode() has to be overloaded in the derived graph class such that 
* it returns a pointer to that newly defined but derived node data type. 
* The same holds for the edges.
*
*/


template<typename T,typename S> class DGraph : public virtual Base_Graph<T,S>
{
  public:
    //cto creating a graph with the same dimension for node and edge weight vectors
                               DGraph(const int& dim) : Base_Graph<T,S>(dim){}
    //cto creating a graph with the different dimension for node and edge weight vectors
                               DGraph(const int& dimNode,const int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge){}  
    //std dto
                               ~DGraph(){}
    //returns whether the two nodes are connected by an edge
    inline bool                areConnected(const int&,const int&)const;
    // adds an edge between the nodes given by their indeces      
    void                       addEdge(const int&,const int&);
    // removes an edge between the nodes given by their indeces if there exists one     
    void                       rmEdge(const int&,const int&);
    // returns a pointer to a edge of a type that is currently used by the graph
    virtual DNode<T,S>*        newNode();
    void                       getID();

};  


/** \brief returns whether the two nodes are connected by an edge either from x to y or from y to x
*/
template<typename T,typename S> inline bool DGraph<T,S>::areConnected(const int& x,const int& y)const
{
  return ( this->_nodes[x]->edges[y] != NULL && this->_nodes[y]->edges[x] != NULL);
} 

/** \brief Adds an edge between the nodes given by their indeces 
*
*   If there is no edge from x to y then it is added and the corresponding number 
*   of connections in the ingoing and outgoing node are updated.
*   No edge from y to x is added.
* 
*   \param x outgoing node
*   \param y ingoing node
*/
template<typename T,typename S> void DGraph<T,S>::addEdge(const int& x,const int& y)
{

  if ( this->_nodes[x]->edges[y] == NULL ) //self edges are allowed
  {
    Base_Edge<S,T>* new_edge = this->newEdge();
    new_edge->in             = this->_nodes[x];
    new_edge->out            = this->_nodes[y];
    new_edge->weight.resize(this->_dimEdge);
    this->_nodes[x]->edges[y]    =  new_edge;   
    // DANGER DownCast is performed via dynamic_cast
    (dynamic_cast< DNode<T,S>* >(this->_nodes[x]))->incOutEdges(); 
    (dynamic_cast< DNode<T,S>* >(this->_nodes[y]))->incInEdges(); 
  }      
  
}  

/** \brief Removes an edge between the nodes given by their indeces if there exists one 
*
*   If there is an edge going from x to y then it is removed and the corresponding 
*   number of connections in the ingoing and outgoing node are updated.
*   It does not remove an edge from y to x.
* 
*   \param x outgoing node
*   \param y ingoing node
*/
template<typename T,typename S> void DGraph<T,S>::rmEdge(const int& x,const int& y)
{
  if ( this->_nodes[x]->edges[y] != NULL ) //self edges are allowed
  {
    delete this->_nodes[x]->edges[y];
    this->_nodes[x]->edges[y]     =  NULL;
    // DANGER DownCast is performed via dynamic_cast  
    (dynamic_cast< DNode<T,S>* >(this->_nodes[x]))->decOutEdges(); 
    (dynamic_cast< DNode<T,S>* >(this->_nodes[y]))->decInEdges();     
  }  
}

/** \brief returns a pointer to a edge of a type that is currently used by the graph
*
* This virtual function allows the use of a (sub)class specific edge type.
* If subclasses want to use another type of edge, then this edge has to be derived from
* the struct Base_Edge<S,T> and the function newEdge() has to be reimplemented in the subclass
* such that the reimplemented function returns a pointer to that new defined and derived 
* edge type.
* This function has to be called by the function addEdge() and inheritance rules guarantee
* that the function of the subclass and not of the superclass is called within the function
* addEdge(),resulting in the use of the used defined edge type as edge for the graph structure.
*/

template<typename T,typename S> DNode<T,S>* DGraph<T,S>::newNode()
{
  DNode<T,S>* n = new DNode<T,S>;
  return n;
}

template<typename T,typename S> inline void DGraph<T,S>::getID()
{
  std::cout << typeid((*this)[0]).name()<<std::endl;
}  

} // namespace neuralgas

#endif
