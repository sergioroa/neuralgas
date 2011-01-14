/** 
* \class TGraph
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

//todo decAge below zero check and following edge removal

#ifndef TGRAPH_H
#define TGRAPH_H

#include "Base_Graph.h"
#include <typeinfo>

namespace neuralgas {

/** \brief The derived edge for the time / age graph
*
*   The edge permits giving each edge a time / age which can be
*   decreased, increased, set and get directly by correspondig functions.
*
* \param age represents the age of the edge which is initialised to 0.0
*/

template<typename A,typename B> struct TEdge : Base_Edge<A,B>
{
  TEdge()
  {age=0.0;}

  //represents the age of the edge which is initialised to 0.0
  float age;
  
  // sets the age of the edge
  void setAge(const float new_age){age=new_age;}
  // increases the age by one
  void incAge(){++age;}
  // decreases the age by one
  void decAge(){--age;}
  // returns the age
  float getAge() const {return age;}

};  



/** \brief TGraph provides a structure for a time or aging graph with n-dim weighted 
*   nodes and m-dim weighted edges
*
* The Base_Graph class is extended such that edge aging is implemented 
* representing an aging ot timed graph.
* The class is abstract, therefore it cannot be instantiated and has to be used as
* a superclass, perhaps in combination with another graph like DGraph or UGraph.
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

template<typename T,typename S> class TGraph : public virtual Base_Graph<T,S>
{ 
  public:
    //cto creating a graph with the same dimension for node and edge weight vectors
    TGraph(const int& dim):Base_Graph<T,S>(dim){}
    //cto creating a graph with the different dimension for node and edge weight vectors
    TGraph(const int& dimNode,const int& dimEdge) : Base_Graph<T,S>(dimNode,dimEdge){}  
    //get the age of the edge going from the first parameter to the second
    float                            getAge(const int&, const int&) const;
    //set the age of the edge going from the first parameter to the second
    void                             setAge(const int&, const int&,const float&);
    //inc the age of the edge going from the first parameter to the second
    void                             incAge(const int&,const int&);
    //dec the age of the edge going from the first parameter to the second
    void                             decAge(const int&,const int&);
    //returns a pointer to an edge of a type that is currently used by the graph
    virtual TEdge<S,T>*              newEdge();
    //abstract addEdge  
    virtual void                     addEdge(const int&, const int&)=0;
    //abstract rmEdge
    virtual void                     rmEdge(const int&, const int&)=0;
    void                             getID(const int&, const int&);
};  
/*
template<typename T,typename S> void TGraph<T,S>::addEdge(const int& x,const int& y)
{
 if ( this->_nodes[x]->edges[y] == NULL ) //self edges are allowed
  {
    TEdge<S,T>* new_edge = this->newEdge();
    new_edge->in             = this->_nodes[x];
    new_edge->out            = this->_nodes[y];
    new_edge->weight.resize(this->_dimEdge);
    this->_nodes[x]->edges[y]    =  new_edge;
  }
}  

template<typename T,typename S> void TGraph<T,S>::rmEdge(const int& a,const int& b)
{
} 
*/
/** Get the age of the edge going from x to y if it does exist otherwise returns -1
*
*   \param x source of the edge
*   \param y destination of the edge
*/
template<typename T,typename S> inline float TGraph<T,S>::getAge(const int& x, const int& y) const
{
// DANGER DownCast is performed via dynamic_cast
 if ( this->_nodes[x]->edges[y] !=NULL && 0 <= x && x < this->size() && 0 <= y && y < this->size() )
     return (dynamic_cast< TEdge<S,T>* >(this->_nodes[x]->edges[y]) )->getAge();
 else
     return -1;
}

/** Set the age of the edge going from x to y if it does exists
*
*   \param x source of the edge
*   \param y destination of the edge
*/
template<typename T,typename S> void TGraph<T,S>::setAge(const int& x, const int& y,const float& value)
{
 // DANGER DownCast is performed via dynamic_cast  
  if ( 0 <= x && x < this->size() && 0 <= y && y < this->size() && value >=0)
  {
   addEdge(x,y);
   (dynamic_cast< TEdge<S,T>* >(this->_nodes[x]->edges[y]) )->setAge(value); 
  }
  else rmEdge(x,y);
}

/** Dec the age of the edge going from x to y if it does exists
*
*   \param x source of the edge
*   \param y destination of the edge
*/
template<typename T,typename S> void TGraph<T,S>::decAge(const int& x,const int& y)
{
 // DANGER DownCast is performed via dynamic_cast
  addEdge(x,y);
  if ( 0 <= x && x < this->size() && 0 <= y && y < this->size() )
     (dynamic_cast< TEdge<S,T>* >(this->_nodes[x]->edges[y]) )->decAge();
}

/** Inc the age of the edge going from x to y if it does exists
*
*   \param x source of the edge
*   \param y destination of the edge
*/
template<typename T,typename S> void TGraph<T,S>::incAge(const int& x,const int& y)
{
 // DANGER DownCast is performed via dynamic_cast
  addEdge(x,y);
  if ( 0 <= x && x < this->size() && 0 <= y && y < this->size() )
  {
     (dynamic_cast< TEdge<S,T>* >(this->_nodes[x]->edges[y]) )->incAge();
  }
}

/** \brief returns a pointer to a edge of a type that is currently used by the graph
*
* This virtual function allows the use of a (sub)class specific edge type.
* If subclasses want to use another type of edge, then this edge has to be derived from
* the struct TEdge<S,T> and the function newEdge() has to be reimplemented in the subclass
* such that the reimplemented function returns a pointer to that new defined and derived 
* edge type.
* This function has to be called by the function addEdge() and inheritance rules guarantee
* that the function of the subclass and not of the superclass is called within the function
* addEdge(),resulting in the use of the used defined edge type as edge for the graph structure.
*/
template<typename T,typename S> TEdge<S,T>* TGraph<T,S>::newEdge()
{
  TEdge<S,T>* edge = new TEdge<S,T>;
  return edge;
}  

template<typename T,typename S> inline void TGraph<T,S>::getID(const int& x, const int& y)
{
  std::cout << typeid( (*this->_nodes[x]->edges[y])).name()<<std::endl;
}  

} // namespace neuralgas

#endif

