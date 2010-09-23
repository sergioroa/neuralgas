#ifndef TGRAPH_H
#define TGRAPH_H

#include "Base_Graph.h"


template<typename A,typename B> struct TEdge;

template<typename A,typename B> struct TEdge : Base_Edge<A,B>
{
  TEdge(){h=0;}
  int h;
  void update()
  {h++;
  std::cout << h <<std::endl;}
  
};  

template<typename T,typename S> class TGraph : virtual public Base_Graph<T,S>
{ 
  public:
    TGraph():Base_Graph<T,S>(10){}
    virtual inline int               getAge(const int&, const int&) const;
    virtual void                     incAge(const int&,const int&);
    TEdge<S,T>*                      getEdge();  
    void                             addEdge(const int&, const int&);
    void                             rmEdge(const int&, const int&);
    void                             getID();
};  

template<typename T,typename S> void TGraph<T,S>::addEdge(const int& a,const int& b)
{
}  

template<typename T,typename S> void TGraph<T,S>::rmEdge(const int& a,const int& b)
{
} 

template<typename T,typename S> inline int TGraph<T,S>::getAge(const int& x, const int& y) const
{
 if ( 0 <= x && x < size() && 0 <= y && y < size() )
 //    return (*this)[x].edge[y].h;
 //else
     return -1;
}
template<typename T,typename S> void TGraph<T,S>::incAge(const int& x,const int& y)
{
  if ( 0 <= x && x < size() && 0 <= y && y < size() )
     (*this)[x][y]->update();
}

template<typename T,typename S> inline TEdge<S,T>* TGraph<T,S>::getEdge()
{
  TGraph<S,T>* edge = new TGraph<S,T>;
  
  return edge;
}  

template<typename T,typename S> inline void TGraph<T,S>::getID()
{
  std::cout << typeid((*this)[0]).name()<<std::endl;
}  
#endif

