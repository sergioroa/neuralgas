#ifndef MGNGGRAPH_H
#define MGNGGRAPH_H

#include "Graph.h"

/** \brief The specific node derived from struct Base_Node<T> for the MGNG graph derived 
* from class Graph<T>
*
* \param context
* \param counter
*/
template <typename T> struct MGNGNode : Base_Node<T>
{ 
 vektor<T> context;
 float counter;  
};

template<typename T> class MGNGGraph : public Graph<T>
{
  public:
                            MGNGGraph(const int&);
                            ~MGNGGraph();   
  inline int                getAge(const int&, const int&) const;
  void                      incAge(const int&,const int&);
};

template<typename T> MGNGGraph<T>::MGNGGraph(const int& dim) : Graph(dim)
{}

template<typename T> MGNGGraph<T>::~MGNGGraph()
{}

template<typename T> inline int MGNGGraph<T>::getAge(const int& x, const int& y) const
{
 if ( 0 <= x && x < size() && 0 <= y && y < size() )
     return (x < y) ? (*_graph)[x][y] : (*_graph)[y][x];
 else
     return -1;
}
template<typename T> void MGNGGraph<T>::incAge(const int& x,const int& y)
{
 if ( 0 <= x && x < size() && 0 <= y && y < size() )
    (x < y) ? ((*_graph)[x][y])++ : ((*_graph)[y][x])++;
      
}

#endif
