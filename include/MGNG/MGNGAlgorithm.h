#ifndef MGNGALGORITHM_H
#define MGNGALGORITHM_H

#include <vector>
#include <math.h>

#include "vektor.h"
#include "Graph.h"

template<typename T> struct mypair
{
 T x;
 T y;
};

template<typename T> void sort(std::vector< mypair<T> >& to_sort)
{
 
}

template<typename T> class MGNGAlgorithm
{
 public:
                                      MGNGAlgorithm(const std::vector< vektor<T> >&,const float&,
                                                    const float&,const int&,const float&,const float&,
                                                    const float&,const float&,const int&,const int&,
                                                    const int&);
                                      ~MGNGAlgorithm();
         std::vector< vektor<T> >     getData() const;
         float                        getAlpha() const;
         float                        getBeta() const;
         int                          getGamma() const;
         float                        getDelta() const;
         float                        getEpsilon_n() const;
         float                        getEpsilon_b() const;
         float                        getEta() const;
         int                          getTheta() const;
         int                          getLambda() const;
         int                          getNeighbors() const;

         void                         setData(const std::vector< vektor<T> >&);
         void                         setAlpha(const float&);
         void                         setBeta(const float&);
         void                         setGamma(const int&);
         void                         setDelta(const float&);
         void                         setEpsilon_n(const float&);
         void                         setEpsilon_b(const float&);
         void                         setEta(const float&);
         void                         setTheta(const int&);
         void                         setLambda(const int&);
         void                         setNeighbors(const int&);
         
         void                         run();
         void                         run(T (*metric_to_use)(const vektor<T>&,const vektor<T>&));
                         
 private:
         T                            metric(const vektor<T>&,const vektor<T>&);
         Graph<T>*                    _graph;
         std::vector< vektor<T> >     _data;
         float                        _alpha;
         float                        _beta;
         int                          _gamma;
         float                        _delta;
         float                        _epsilon_n;
         float                        _epsilon_b;
         float                        _eta;
         int                          _theta;
         int                          _lambda;
         int                          _neighbors;

};

template<typename T> MGNGAlgorithm<T>::MGNGAlgorithm(const std::vector< vektor<T> >& data,const float& alpha, 
                                const float& beta, const int& gamma,const float& delta,
                                const float& epsilon_n, const float& epsilon_b,
                                const float& eta, const int& theta,const int& lambda,
                                const int& neighbors=2)
{
 int size=data.size();
 _data.resize(size);
 for(int i=0; i < size; i++)
         _data[i]=data[i];
 
 _alpha                    =     alpha;
 _beta                     =     beta;
 _gamma                    =     gamma;
 _delta                    =     delta;
 _epsilon_n                =     epsilon_n;
 _epsilon_b                =     epsilon_b;
 _eta                      =     eta;
 _theta                    =     theta;
 _lambda                   =     lambda;
 _neighbors                =     neighbors;
 _graph                    =     NULL;
}

template<typename T> MGNGAlgorithm<T>::~MGNGAlgorithm()
{
 if (_graph != NULL)
    delete _graph;
}

template<typename T> std::vector<vektor<T> > MGNGAlgorithm<T>::getData() const
{
 return _data;
}

template<typename T> float MGNGAlgorithm<T>::getAlpha() const
{
 return _alpha;
}

template<typename T> float MGNGAlgorithm<T>::getBeta() const
{
 return _beta;
}

template<typename T> int MGNGAlgorithm<T>::getGamma() const
{
 return _gamma;
}

template<typename T> float MGNGAlgorithm<T>::getDelta() const
{
 return _delta;
}

template<typename T> float MGNGAlgorithm<T>::getEpsilon_n() const
{
 return _epsilon_n;
}

template<typename T> float MGNGAlgorithm<T>::getEpsilon_b() const
{
 return _epsilon_b;
}

template<typename T> float MGNGAlgorithm<T>::getEta() const
{
 return _eta;
}

template<typename T> int MGNGAlgorithm<T>::getTheta() const
{
 return _theta;
}

template<typename T> int MGNGAlgorithm<T>::getLambda() const
{
 return _lambda;
}

template<typename T> int MGNGAlgorithm<T>::getNeighbors() const
{
 return _neighbors;
}

template<typename T> void MGNGAlgorithm<T>::setData(const std::vector< vektor<T> >& data)
{
 _data.resize(data.size());
 for(int i=0; i < data.size(); i++)
         _data[i]=data[i];
}

template<typename T> void MGNGAlgorithm<T>::setAlpha(const float& alpha)
{
 _alpha=alpha;
}

template<typename T> void MGNGAlgorithm<T>::setBeta(const float& beta)
{
 _beta=beta;
}

template<typename T> void MGNGAlgorithm<T>::setGamma(const int& gamma)
{
 _gamma=gamma;
}

template<typename T> void MGNGAlgorithm<T>::setDelta(const float& delta)
{
 _delta=delta;
}

template<typename T> void MGNGAlgorithm<T>::setEpsilon_n(const float& epsilon_n)
{
 _epsilon_n=epsilon_n;
}

template<typename T> void MGNGAlgorithm<T>::setEpsilon_b(const float& epsilon_b)
{
 _epsilon_b=epsilon_b;
}

template<typename T> void MGNGAlgorithm<T>::setEta(const float& eta)
{
 _eta=eta;
}

template<typename T> void MGNGAlgorithm<T>::setTheta(const int& theta)
{
 _theta=theta;
}

template<typename T> void MGNGAlgorithm<T>::setLambda(const int& lambda)
{
 _lambda=lambda;
}

template<typename T> void MGNGAlgorithm<T>::setNeighbors(const int& neighbors)
{
 _neighbors=neighbors;
}

template<typename T> void MGNGAlgorithm<T>::run(T (*metric_to_use)(const vektor<T>& x,const vektor<T>& y) )
{
 int dimension = _data[0].size();
 
 if (dimension>0)
 {
  // line 1 - 3
  _graph                          =     new Graph<T>(dimension);
  
  // number of total steps for the algorithm
  int tsize                       =     _data.size();
  // creating a zero of type T
  T zero
                            =     _data[0][0]-_data[0][0];
  // line 4
  vektor<T> globalContextV;
  globalContextV.resize(dimension,zero);
  
  //line 5-6
  for(int t = 0; t < tsize; t++)
  {
   int first_winner               =     0;
   int second_winner              =     0;
   T distance                     =     zero;
   T best_distance                =     zero;
   
   for (int j = 0; j < (*_graph).size(); j++)
   {
    
    distance= (1 - _alpha) * pow((*metric_to_use)(_data[t],(*_graph)[j].weight_v) ,2);
    distance+= _alpha*pow((*metric_to_use)(globalContextV,(*_graph)[j].context_v) ,2);
    
    if (distance <= best_distance)
    {
     second_winner              =       first_winner;
     first_winner               =       j;
     best_distance              =       distance;
    }       
   }
   
   // line 7
   (*_graph)[first_winner].counter++;
   
   // line 10
   std::vector<int> first_winner_neighbors = (*_graph).getNeighbors(first_winner);
   for(int j=0; j < first_winner_neighbors.size();j++)
           (*_graph).incAge(first_winner, first_winner_neighbors[j] ); 
   
   // line 8 - 9
   // line 10 is realised before 8-9 for efficiency reasons, in line 10 we had otherwise
   // to check in each for-loop with an if statement whether the neighbor is the second_winner
   // in order to not increment the edge age
   // therefore we possibly (if the second_winner was already a neighbor) increment the age first
   // and set it then to 0 by creating an edge ( if there was already an edge, nothing is changed 
   // by the command )
   
   (*_graph).addEdge(first_winner,second_winner);
   
   // line 11
   for (int j = 0; j < (*_graph).size(); j++)
   {
    for (int k = 0; k < j; k++)
        if ( (*_graph).getAge(k,j) > _gamma )
           (*_graph).rmEdge(k,j);
   }
   
   // line 12
   for (int j = 0; j < (*_graph).size(); j++)
       if (!(*_graph).isConnected(j)) 
          (*_graph).deleteNode(j);
   
   // line 13
   (*_graph)[first_winner].weight_v  += 
               _epsilon_b * ( _data[t]-(*_graph)[first_winner].weight_v);
   (*_graph)[first_winner].context_v += 
               _epsilon_b * ( globalContextV - (*_graph)[first_winner].context_v);
   for(int j=0; j < first_winner_neighbors.size();j++)
   {
      (*_graph)[first_winner_neighbors[j]].weight_v  += 
               _epsilon_n * (globalContextV - (*_graph)[first_winner_neighbors[j]].weight_v);
      (*_graph)[first_winner_neighbors[j]].context_v +=
               _epsilon_n * (globalContextV - (*_graph)[first_winner_neighbors[j]].context_v);
   } 
   
   // line 14
   globalContextV=( 1 - _beta)*(*_graph)[first_winner].weight_v + _beta * (*_graph)[first_winner].context_v;
 
   // line 15
   if( t % _lambda == 0 && (*_graph).size() < _theta)
   {
       float max_counter                   = 0.0;
       int max_counter_index               = 0;
       for (int i=0; i < (*_graph).size(); i++)
       {
           if ((*_graph)[i].counter > max_counter )
           {
            max_counter             = (*_graph)[i].counter;
            max_counter_index       = i;
           }
       }
       std::vector<int> neighbors   = (*_graph).getNeighbors(max_counter_index);
       float max_counter_n                 = 0.0;
       int max_counter_index_n             = 0;

       for (int i=0; i < neighbors.size(); i++)
       {
           if ((*_graph)[neighbors[i]].counter > max_counter_n )
           {
            max_counter_n             = (*_graph)[neighbors[i]].counter;
            max_counter_index_n       = i;
           }
       }
       
       (*_graph).newNode();
       int gsize                      = (*_graph).size();
       (*_graph)[gsize].weight_v      = (*_graph)[max_counter_index_n].weight_v
                                       +(*_graph)[max_counter_index].weight_v;
       (*_graph)[gsize].context_v    = (*_graph)[max_counter_index_n].context_v
                                       +(*_graph)[max_counter_index].context_v;
                                    
   } 
   
   // line 16
   for (int i = 0; i < (*_graph).size(); i++)
   {
       (*_graph)[i].counter*=_eta;
   }
  }
 }
                  
}

template<typename T> void MGNGAlgorithm<T>::run()
{
 int dimension = _data[0].size();
 
 if (dimension>0)
 {
  // line 1 - 3
  _graph                          =     new Graph<T>(dimension);
  
  // number of total steps for the algorithm
  int tsize                       =     _data.size();
  // creating a zero of type T
  T zero
                            =     _data[0][0]-_data[0][0];
  // line 4
  vektor<T> globalContextV;
  globalContextV.resize(dimension,zero);
  
  //line 5-6
  for(int t = 0; t < tsize; t++)
  {
   int first_winner               =     0;
   int second_winner              =     0;
   T distance                     =     zero;
   T best_distance                =     zero;
   
   for (int j = 0; j < (*_graph).size(); j++)
   {
    
    distance= (1 - _alpha) * pow(metric(_data[t],(*_graph)[j].weight_v) ,2);
    distance+= _alpha*pow(metric(globalContextV,(*_graph)[j].context_v) ,2);
    
    if (distance <= best_distance)
    {
     second_winner              =       first_winner;
     first_winner               =       j;
     best_distance              =       distance;
    }       
   }
   
   // line 7
   (*_graph)[first_winner].counter++;
   
   // line 10
   std::vector<int> first_winner_neighbors = (*_graph).getNeighbors(first_winner);
   for(int j=0; j < first_winner_neighbors.size();j++)
           (*_graph).incAge(first_winner, first_winner_neighbors[j] ); 
   
   // line 8 - 9
   // line 10 is realised before 8-9 for efficiency reasons, in line 10 we had otherwise
   // to check in each for-loop with an if statement whether the neighbor is the second_winner
   // in order to not increment the edge age
   // therefore we possibly (if the second_winner was already a neighbor) increment the age first
   // and set it then to 0 by creating an edge ( if there was already an edge, nothing is changed 
   // by the command )
   
   (*_graph).addEdge(first_winner,second_winner);
   
   // line 11
   for (int j = 0; j < (*_graph).size(); j++)
   {
    for (int k = 0; k < j; k++)
        if ( (*_graph).getAge(k,j) > _gamma )
           (*_graph).rmEdge(k,j);
   }
   
   // line 12
   for (int j = 0; j < (*_graph).size(); j++)
       if (!(*_graph).isConnected(j)) 
          (*_graph).deleteNode(j);
   
   // line 13
   (*_graph)[first_winner].weight_v  += 
               _epsilon_b * ( _data[t]-(*_graph)[first_winner].weight_v);
   (*_graph)[first_winner].context_v += 
               _epsilon_b * ( globalContextV - (*_graph)[first_winner].context_v);
   for(int j=0; j < first_winner_neighbors.size();j++)
   {
      (*_graph)[first_winner_neighbors[j]].weight_v  += 
               _epsilon_n * (globalContextV - (*_graph)[first_winner_neighbors[j]].weight_v);
      (*_graph)[first_winner_neighbors[j]].context_v +=
               _epsilon_n * (globalContextV - (*_graph)[first_winner_neighbors[j]].context_v);
   } 
   
   // line 14
   globalContextV=( 1 - _beta)*(*_graph)[first_winner].weight_v + _beta * (*_graph)[first_winner].context_v;
 
   // line 15
   if( t % _lambda == 0 && (*_graph).size() < _theta)
   {
       float max_counter                   = 0.0;
       int max_counter_index               = 0;
       for (int i=0; i < (*_graph).size(); i++)
       {
           if ((*_graph)[i].counter > max_counter )
           {
            max_counter             = (*_graph)[i].counter;
            max_counter_index       = i;
           }
       }
       std::vector<int> neighbors   = (*_graph).getNeighbors(max_counter_index);
       float max_counter_n                 = 0.0;
       int max_counter_index_n             = 0;

       for (int i=0; i < neighbors.size(); i++)
       {
           if ((*_graph)[neighbors[i]].counter > max_counter_n )
           {
            max_counter_n             = (*_graph)[neighbors[i]].counter;
            max_counter_index_n       = i;
           }
       }
       
       (*_graph).newNode();
       int gsize                      = (*_graph).size();
       (*_graph)[gsize].weight_v      = (*_graph)[max_counter_index_n].weight_v
                                       +(*_graph)[max_counter_index].weight_v;
       (*_graph)[gsize].context_v    = (*_graph)[max_counter_index_n].context_v
                                       +(*_graph)[max_counter_index].context_v;
                                    
   } 
   
   // line 16
   for (int i = 0; i < (*_graph).size(); i++)
   {
       (*_graph)[i].counter*=_eta;
   }
  }
 }
}

template<typename T> T MGNGAlgorithm<T>::metric(const vektor<T>& x,const vektor<T>& y)
{
 T result;
 T value;
 
 result =  x[0] - x[0];
 
 for (int i=0; i < x.size(); i++)
 {
     value = (x[i]*x[i] - y[i]*y[i]);
     result+= ( value>0 ) ? value : -value;
 }
    
 return sqrt(result); 

}

#endif
