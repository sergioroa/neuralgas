/** 
* \file Base_Graph.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/


#ifndef BASE_GRAPH_H
#define BASE_GRAPH_H

#include <vector>
#include <map>
#include <fstream>
#include "Vector.h"
#include <iostream>
#include <cstdlib>
#include <tools/metrics.h>
#include <tools/helpers.h>

namespace neuralgas {

// TODO check edge deleting for derived classes (or better use shared_ptr)

template < typename T , typename S > struct Base_Node;
template < typename T , typename S> struct Base_Edge;

template< typename T, typename S > class NeuralGas;

/** \class Base_Node
 *  \brief The base class for nodes
 *
 * \param num_connections number of connected nodes to this node
 * \param weight vectors as data vector
 * \param func user definded function called in the function update() for changing the node's values
 * \param edges a vector of ptrs to the node's edges 
 */
template < typename T , typename S > struct Base_Node
{
  Base_Node()
  {
   num_connections=0;
   func=NULL;  
  }

  virtual ~Base_Node()
  { 
    int esize = edges.size(); 
    for(int i = 0; i < esize; i++)
        if (edges[i]!=NULL)
        {
		delete edges[i];
                // edges[i] = NULL;
        }  
    
  }
  /** \brief Calls the user defined function and applies it on this node
  */
  virtual void update(const unsigned int& index=0){(*func)(this);}  
  //user defined function called in update for changing the node's values
  void (*func) (Base_Node<T,S>*);
  // number of connected nodes to this node
  int num_connections;
  // weight vectors as data vector
  Vector<T> weight;
  // edge vector containing the edges of the node  
  std::vector< Base_Edge<S,T>* > edges;

};  

/** \class Base_Edge
 *  \brief The base class for edges
 *
 * \param weight vectors as data vector
 * \param in ptr to the ingoing node 
 * \param out ptr to the outgoing node
 */
template < typename T , typename S> struct Base_Edge
{
  
  Base_Edge(){in = NULL; out=NULL;}
  ~Base_Edge(){//in = NULL; out=NULL;
                    }

  /** \brief Calls the user defined function and applies it on this edge
  */
  virtual void update(){(*func)(this);}
  //user definded function called in update for changing the edges's values
  void (*func) (Base_Edge<T,S>*);
  // weight vectors as data vector
  Vector<T> weight;
  //ptr to the ingoing node
  Base_Node<S,T>* in;
  //ptr to the outgoing node
  Base_Node<S,T>* out;
};  

/** \class Base_Graph
 *  \brief Base_Graph provides a basic structure for a graph with n-dim weighted nodes and m-dim weighted edges
 *
 * The graph structure comes up with adding, removing nodes; 
 * initializing the graph with a given number of nodes where the weight vectors have 
 * random values. Further functions are checking whether a node has neighbors and
 * getting the neighbors function. 
 * Kind of iterator functions allow to apply a given user defined function to all
 * nodes or to the neighbor nodes of a node. 
 * Edge handling like adding and removing has to be treated in the derived classes in order to
 * ensure the users requirements to implement either an undirected or a directed graph.
 * It is recommended in the case of an undirected graph to use for two connected nodes exactly 
 * one edge that is refered to from the pointers within both nodes.
 * And in the case for directed graphs it is expected to use up to two edges, each edge referring 
 * to the in node and corresponding out node. 
 * The first template parameter is the type of the node weigth vectors and the second template
 * parameter is the type of the edge weigth vectors.
 * The following paragraph is important for derived classes. 
 * It is possible to use a user defined type of node for the graph.
 * In this case that newly defined node has to be derived from the struct Base_Node<T,S>
 * and the virtual function newNode() has to be overloaded in the derived graph class such that 
 * it returns a pointer to that newly defined but derived node data type. 
 * The same holds for the edges.
 *
 * \param _nodes array of pointer to the nodes of the graph
 * \param _dimNode dimension of the node's weight vectors
 * \param _dimEdge dimension of the edge's weight vectors
 */
template<typename T, typename S> class Base_Graph
{
      public:
	     typedef T (Base_Graph::*Metric)(const Vector<T>&,const Vector<T>&) const;
             //cto creating a graph with the same dimension for node and edge weight vectors
             Base_Graph(const unsigned int&);
             //dummy copy constructor. Copying procedures should be done in derived classes
             Base_Graph (const Base_Graph&) {}
                         
             //cto creating a graph with the different dimension for node and edge weight vectors
             Base_Graph(const unsigned int&,const unsigned int&);
             
             //std dto
             virtual ~Base_Graph();            
             //returning a reference to the node indexed by the given index
             Base_Node<T,S>&                     operator[](const unsigned int&);
             //returning a const reference to the node indexed by the given index
             Base_Node<T,S>&                     operator[](const unsigned int&) const;
             //inits the graph with the given number of nodes with random valued weight vectors
	     void                                initRandomGraph(const unsigned int&,const Vector<T>&, const Vector<T>&);
	//adds a new uninitialized, edgeless node into the graph
             void                                addNode(void);
             // removes the node given by the index, removes its edges and updates the number of connections of its neighbors
             virtual void                        rmNode(const unsigned int&); 
             // adds an edge between the nodes given by their indeces      
             virtual void                        addEdge(const unsigned int&,const unsigned int&)=0;
             // removes an edge between the nodes given by their indeces if there exists one     
             virtual void                        rmEdge(const unsigned int&,const unsigned int&)=0;
             // saves the nodes weight in a file
             bool                                save(const char*, bool t = false);
             // loads nodes 
             void                                setNodes( std::vector < Base_Node<T, S>* >* nodes);
             //returns a vector of ints representing the indices of the neighboring nodes
             std::vector<unsigned int>           getNeighbors(const unsigned int&) const;
             // Returns neighbors set cardinality for some node
             int                                 getNeighborsSize(const unsigned int&) const;
             //returns whether the node of the given index has neighbors or not
             inline bool                         isConnected(const unsigned int&) const;
             //applies the given func to all nodes
             inline void                         applyFunc2AllNodes(void (*func)(Base_Node<T,S>*,const float&),const float&);
             //applies the given func to the neighboring nodes 
             inline void                         applyFunc2Neighbors(const unsigned int&,void (*func)(Base_Node<T,S>*,const float&),const float&);
             //calls the update function declared within the node struct
             inline void                         update();
             //returns the number of nodes currently in the graph
             inline unsigned int                 size(void) const;
	     //sets a user defined metric, used as distance of reference and data vector 
             inline void                         setMetric(Metric); 
             // pre-specified metric is the standard L2 euclidean metric
             virtual T                           metric(const Vector<T>&, const Vector<T>&) const;

             // Returns distances from a node to a set of nodes
             std::map<unsigned int, T>           get1toMDistances (const unsigned int&, std::vector<unsigned int>&);

             void                                showGraph();

	     // get nodes vector
	     inline std::vector< Base_Node<T,S>* >* getNodes ();

      protected:
             //returns a pointer to a node of a type that is currently used by the graph
             virtual Base_Node<T,S>*             newNode();
             //returns a pointer to an edge of a type that is currently used by the graph
             virtual Base_Edge<S,T>*             newEdge();
             // array of pointer to the nodes of the graph
             std::vector< Base_Node<T,S>* >      _nodes;            
            
             // dimension of the node's weight vectors
             unsigned int                        _dimNode;
             // dimension of the edge's weight vectors
             unsigned int                        _dimEdge;
             //user specified metric
             Metric _metric_to_use;
	//private:

};

template < typename T , typename S> void Base_Graph<T,S>::showGraph()
{
 std::cout << "Graphsize " << _nodes.size() << std::endl;

 for(unsigned int i=0; i < _nodes.size(); i++)
 {        
          for(unsigned int j=0; j < _nodes[i]->edges.size(); j++)
                  std::cout <<( ( _nodes[i]->edges[j]==NULL) ? " " : "*");
          std::cout<<std::endl;
 }
 std::cout << std::endl;
}
/*
template<typename T,typename S> void Base_Graph<T,S>::addEdge(const unsigned int& a,const unsigned int& b)
{
}  

template<typename T,typename S> void Base_Graph<T,S>::rmEdge(const unsigned int& a,const unsigned int& b)
{
}  
*/
/** \brief cto creating a graph with the given dimension of the weight vectors where 
* the values of the vectors are chosen randomly
*
* \param dim is the dimension of the weight vectors of the nodes in the graph
*/
template<typename T,typename S> Base_Graph<T,S>::Base_Graph(const unsigned int& dim)
{ 
  _dimNode = dim;
  _dimEdge = dim;
  _metric_to_use  = NULL;   // no user defined metric is used, std metric is used
  ::srand( (unsigned)time( NULL ) );                    //inits the random function 
}

/** \brief cto creating a graph with the given dimension of the node weight vectors where 
* the values of the vectors are chosen randomly and given possibly different dimension of the 
* edge vectors.
*
* \param dimNode is the dimension of the weight vectors of the nodes in the graph
* \param dimEdge is the dimension of the weight vectors of the edges in the graph
*/
template<typename T,typename S> Base_Graph<T,S>::Base_Graph(const unsigned int& dimNode,const unsigned int& dimEdge)
{
  _dimNode = dimNode;
  _dimEdge = dimEdge; 
  ::srand( (unsigned)time( NULL ) );                    //inits the random function 
}

/** \brief std dto
*/
template<typename T,typename S> Base_Graph<T,S>::~Base_Graph()
{
  int nsize  = size();
  for(int i = 0; i < nsize; i++)
  {
    delete  _nodes[i];                                // delete ptrs to the nodes
    //_nodes[i] = NULL;            
    //_nodes.pop_back();
  }
  _nodes.clear();
  
}

/** \brief returning a reference to the node indexed by the given index
*
* \param index of the node that shall be returned 
*/
template<typename T,typename S> Base_Node<T,S>& Base_Graph<T,S>::operator[](const unsigned int& index)
{
  assert (index < size () );
  return *(_nodes[index]);
}  

/** \brief returning a const reference to the node indexed by the given index
*
* \param index of the node that shall be returned 
*/
template<typename T,typename S> Base_Node<T,S>& Base_Graph<T,S>::operator[](const unsigned int& index) const
{
  assert (index < size());
  return *(_nodes[index]);
}



/** \brief Saves the nodes weight to a file
*
*   The functions saves the weights of the nodes to a file where
*   for n nodes Xi with d dimensional weights the format looks as follows
*   X11 X12 X13 ... X1d
*   .
*   .
*   .
*   Xn1 Xn2 Xn3 ... Xnd
*
* which means that they are separated by spaces.
*
* \param filename is the name of the file where to store the data to
* \param text save the file in text format if true
*/
template<typename T,typename S> bool Base_Graph<T,S>::save( const char* filename, bool text)
{
  if (text)
    return neuralgas::saveNodesText (filename, &_nodes);
  else
    return neuralgas::saveNodes (filename, &_nodes);
}

/** \brief Loads nodes
    \param nodes
 */
template<typename T, typename S> void Base_Graph<T,S>::setNodes( std::vector < Base_Node<T, S>* >* nodes)
{
  _nodes = *nodes;
}



/** \brief Inits the graph with the given number of nodes with random valued weight vectors.
* If a graph already existed it is going to be erased.
* 
* \param num_of_nodes is the number of nodes with random valued weight vectors the graph is going to be initialized with
*/
  
template<typename T,typename S> void Base_Graph<T,S>::initRandomGraph(const unsigned int& num_of_nodes, const Vector<T>& low_limits, const Vector<T>& high_limits)
{    
  unsigned int nsize  = size();
  for(unsigned int i = 0; i < nsize; i++)
  {
    delete  _nodes[i];                                // delete ptrs to the nodes
    _nodes[i] = NULL;            
    _nodes.pop_back();                                // rm ptr from the node array 
  }  
 
  for (unsigned int i = 0; i < num_of_nodes; i++)
  {
     addNode();                                     // adds a new node
     (_nodes[i])->weight.resize(_dimNode);          // sets dimension of the weight vector
     
     for(unsigned int j = 0; j < _dimNode; j++)          
     {
	     (_nodes[i])->weight[j] = (T) ((rand() / (static_cast<T>(RAND_MAX) + 1.0)) * (high_limits[j] - low_limits[j]) + low_limits[j] );   //sets the value of the weights to random values
     }
  }  
}   


/** \brief returns a pointer to a node of a type that is currently used by the graph
*
* This virtual function allows the use of a (sub)class specific node type.
* If subclasses want to use another type of node, then this node has to be derived from
* the struct Base_Node<T,S> and the function newNode() has to be reimplemented in the subclass
* such that the reimplemented function returns a pointer to that new defined and derived 
* node type.
* This function is called by the function addNode() and inheritance rules guarantee
* that the function of the subclass and not of the superclass is called within the function
* addNode(),resulting in the use of the used defined node type as node for the graph structure.
*/

template<typename T,typename S> Base_Node<T,S>* Base_Graph<T,S>::newNode(void)
{
  Base_Node<T,S>* n = new Base_Node<T,S>;
  return n;
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

template<typename T,typename S> Base_Edge<S,T>* Base_Graph<T,S>::newEdge(void)
{
  Base_Edge<S,T>* e = new Base_Edge<S,T>;
  return e;
} 


/** \brief adds a new uninitialized, edgeless node into the graph
* The function uses the newNode function to get the currently (sub) class dependent 
* Node type and adds this node to the graph.
* There are no edges and the weight vectors are uninitialized, solely the dimension
* of the weight vectors is set.
*/
template <typename T,typename S> void Base_Graph<T,S>::addNode(void)
{
  // function does as follows
  // creates a new node of (sub)class specific type and adds it to the _node array
  // a new slot for a possible edge is added to each node
  unsigned int nsize=size();  
  
  Base_Node<T,S>*    n   =   newNode();  // gets the actual used node type 
  n->weight.resize(_dimNode);
  for(unsigned int i=0; i < nsize; i++)
  {
   n->edges.push_back(NULL);    
  }
  _nodes.push_back(n);                  // add the new node
 
  nsize++;
    
  for(unsigned int i=0; i < nsize; i++)
  {
   _nodes[i]->edges.push_back(NULL);    
  }
 
}

/** \brief Removes the node given by the index, removes its edges and updates the number 
* of connections of its neighbors
*
* \param index is the node that shall be deleted
*/
template<typename T,typename S> void Base_Graph<T,S>::rmNode(const unsigned int& index)
{
  // function does as follows
  // checks whether the neighbors have an (directed) edge to index and deletes them
  // deletes node in the _nodes array  
 
  unsigned int     nsize=size();
  assert ( index < nsize );
  std::vector<unsigned int> neighbors = getNeighbors(index);
  
  for(unsigned int i=0; i < neighbors.size(); i++)
    if ( _nodes[ neighbors[i] ]->edges[index]!=NULL ) // edge from i to index
    {
      delete  _nodes[  neighbors[i] ]->edges[index];
      _nodes[ neighbors[i] ]->edges[index]=NULL;   
      _nodes[ neighbors[i] ]->num_connections--;
    }  
  
  for(unsigned int i=0; i < nsize; i++)
  {
    _nodes[ i ]->edges.erase( _nodes[ i ]->edges.begin() + index ); 
  }
  
  delete _nodes[index];                                 // delete ptrs to the nodes
  _nodes[index] = NULL;             
  
  _nodes.erase(_nodes.begin() + index);                 // erase node within the array 
 
}


/** \brief returns a vector of ints representing the indeces of the neighboring nodes
*
* \param index of the homie whose neighbors we are looking for
*/
template<typename T,typename S> std::vector<unsigned int> Base_Graph<T,S>::getNeighbors(const unsigned int& index) const
{
  unsigned int nsize = size();
  assert (index < nsize );

  std::vector<unsigned int> result_v;
  int num_of_neighbors = _nodes[index]->num_connections;
  int found_neighbors = 0;
  
  for (unsigned int i=0; i < nsize /*&& found_neighbors < num_of_neighbors*/; i++)
    if( _nodes[index]->edges[i] !=NULL && index!=i)
    {
        result_v.push_back(i);   
        found_neighbors++;
    }
  //if the following is false, something must be wrong
  assert (found_neighbors <= num_of_neighbors);
  return result_v;
}

/** \brief Returns neighbors set cardinality for some node
 * \param index
 */
template<typename T,typename S> int Base_Graph<T,S>::getNeighborsSize(const unsigned int& index) const
{
  assert ( index < size() );
  int num_of_neighbors = _nodes[index]->num_connections;
  int found_neighbors = 0;
	
  for (unsigned int i=0; i < size() /*&& found_neighbors < num_of_neighbors*/; i++)
    if( _nodes[index]->edges[i] !=NULL && index!=i)
      found_neighbors++;
  //if the following is false, something must be wrong
  assert (found_neighbors <= num_of_neighbors);

  return found_neighbors;

}

/** \brief Returns distances from a node to a set of nodes
 *  \param index index of the node
 *  \param indices indices of the set of nodes
 */
template<typename T, typename S> std::map<unsigned int, T> Base_Graph<T,S>::get1toMDistances (const unsigned int& index, std::vector<unsigned int>& indices)
{
  assert ( index < size() );
  for (unsigned int i=0; i < indices.size(); i++ )
    assert ( indices[i] < size() );

  std::map<unsigned int, T> distances;
  for (unsigned int i=0; i < indices.size(); i++)
    distances[indices[i]] = metric (_nodes[index]->weight, _nodes[indices[i]]->weight);
  return distances;
}

/** \brief Returns whether the node of the given index has neighbors
*
* The function returns true if the node of the given index has neighbors, 
* edges respectively
* and returns false if the node is solely without neighbors, edges respectively
*
* \param index is the index of the node in question for connectedness
*/
template<typename T,typename S> inline bool Base_Graph<T,S>::isConnected(const unsigned int& index) const
{
 if ( index < size() ) 
     return ( _nodes[index]->num_connections > 0 ) ? true : false;
 else
     return false; 
}


/** \brief applies the given func to all nodes
*
* An internal function is provided since an external access by the graph opertor[]
* is slower due to the range check.
* \param func is a pointer to a function that shall be applied to all nodes
*/
template<typename T,typename S> inline void Base_Graph<T,S>::applyFunc2AllNodes(void (*func)(Base_Node<T,S>*,const float&),const float& value)
{
  unsigned int nsize = size();
  for(unsigned int i = 0; i < nsize; i++)
    (*func)(_nodes[i],value);
}  

/** \brief applies the given func to the neighboring nodes 
*
* An internal function is provided since an external access by the graph opertor[]
* is slower due to the range check.
* \param func is a pointer to a function that shall be applied to all neighboring nodes
*/
template<typename T,typename S> inline void Base_Graph<T,S>::applyFunc2Neighbors(const unsigned int& index,void (*func)(Base_Node<T,S>*,const float&),const float& value)
{
  unsigned int nsize = size();
 
  if ( index < nsize )
  {   
    for (unsigned int i=0; i < nsize; i++)
       if ( _nodes[index]->edges[i]!=NULL )
          (*func)(_nodes[i],value);   //applies the func to the node since it is a neighbor          
  
  }
  
}  

/** \brief calls the update function declared within the node struct
*
* The Base_Node<T> has an update function and therefore all derived structs as well.
* It is possible to define for each node an own unique update function reflecting
* different update rules.
* ATTENTION ! Since no check is performed whether an update rule for a node was declared
* one has to take care before calling this function that to each node was
* an update rule declared. Otherwise that will provoce unexpected progam behavior.
*/
template<typename T,typename S> inline void Base_Graph<T,S>::update()
{
  unsigned int nsize = size();
  for(unsigned int i = 0; i < nsize; i++)
    (_nodes[i])->update();
}              

/** \brief returns the number of nodes currently in the graph
*
*/
template<typename T, typename S > inline unsigned int Base_Graph<T,S>::size(void) const
{
 return _nodes.size();
}

//!
/*! 
  \return the vector of nodes
*/
template< typename T, typename S> inline std::vector< Base_Node<T,S>* >* Base_Graph<T,S>::getNodes ()
{
	return &_nodes;
}


/** \brief standard is the pre-specified L2 euclidean metric
*
* \param x first vector
* \param y second vector
*/
template< typename T, typename S> T Base_Graph<T,S>::metric(const Vector<T>& x, const Vector<T>& y) const
{
 if (_metric_to_use==NULL) // is a non-standard metric set ?
   return euclidean<T,S> (x, y);
 else
   return (this->*_metric_to_use)(x,y);   // use the non-standard user defined metric
}  

/** \brief Sets an user defined metric, used as distance of reference and data vector.
*
* The user can define an own metric to measure the distance or distortion of the reference vectors and the data vectors. 
* The metric has to be of the form
*
* T user_metric(const Vector<T>&,const Vector<T>&)
*
* Implicitly this defines a topology.
* The default value of the function pointer is NULL, this allows after an user defined metric
* was given to call the function with no parameter which effects that the default implemented
* metric is going to be used.
* \param *metric_to_use is function ptr to an user defined metric
*/
template < typename T, typename S > inline void Base_Graph<T,S>::setMetric(Metric metric_to_use=NULL )
{
  _metric_to_use=metric_to_use;
}

} //namespace neuralgas

#endif
