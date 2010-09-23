/** 
* \class NeuralGas
* \author Manuel Noll
* 
*  Copyright(c) 20010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/



#ifndef NEURALGAS_H
#define NEURALGAS_H

#include <math.h>
#include <vector>
#include "Graph.h"

/** \brief Class declares the basic operations for the Neural Gas Algorithm and is intended as a super class for inheritance.
*
* The Neural Gas Algortihm as proposed by Martinetz and Schulten, "A Neural Gas Network Learns Topologies"
* served to build a general algorithmic framework.
* \param *graphptr a pointer to the underlying data structure representing the graph
* \param _dimension is dimension of the input data and therefore the dimension of the graph weight vectors
* \param *_data is a pointer to the given input data; a pointer is used, since for large data a duplicate set could be to memory consuming
* \param *_metric_to_use is a function pointer which is NULL if not a particular metric is given via setMetric
*/

template < typename T > class NeuralGas
{
  public:
    //cto with size of dimension as input
                            NeuralGas(const int&);         
    //std dto
                            ~NeuralGas(void);               
    //sets the data to be processed
    inline void             setData(const std::vector<T>*);  
    //sets a user defined metric, used as distance of reference and data vector 
    inline void             setMetric(T (*)(const T& a,const T& b));  
    
    
  protected:
    // pre-specified metric is the standard L2 euclidean metric
    virtual T               metric(const Vector<T>&, const Vector<T>&);
    // sets the number of inital reference vectors
    void                    setRefVectors(const int&);    
    // applies a given function to all neighbors
    void                    applyFunc2Neighbors(const int&, void (*)(Base_Node<T>& n)); 
    // applies a given function to all nodes
    void                    applyFunc2AllNodes(void (*)(Base_Node<T>& n)); 
    // returns the dimension of the data vectors
    inline int              getDimension(void) const; 
    // ptr to the underlying graph
    Graph<T>*               graphptr; 
  
  private:
    //dimension of the vectors
    int                     _dimension; 
    // ptr to the input data
    std::vector<T>*         _data; 
    //user specified metric
    T                       (*_metric_to_use)(const T&,const T&); 
    
      
};  

/** \brief Cto with the dimension of the vectors. 
*
* Cto has as input the to-be-used dimension of the vectors. 
* The dimension has to be the dimension of the input data.
* The dimension cannot be changed after the creation of the object.
* \param dimension is dimension of the input data and therefore the dimension of the graph weight vectors
*/
template < typename T > NeuralGas<T>::NeuralGas(const int& dimension)
{
  _dimension      = dimension;
  _data           = NULL;
  _metric_to_use  = NULL;   // no user defined metric is used, std metric is used
}  

/** \brief Dto frees memory by deleting the underlying graph data structure.
*/
template < typename T > NeuralGas<T>::~NeuralGas(void)
{if (graphptr!=NULL) delete graphptr;}

/** \brief Sets the data that has to be processed in the next algorithmic run.
*
* The function sets the data that has to be processed in the next algorithmic run.
* If a run already took place and the function setRefVectors is not called until the next run
* then the underlying graph persists and will fit to this newly added data, allowing
* a succesiv adaption of the graph.
* \param *data is the given data
*/
template < typename T > inline void NeuralGas<T>::setData(const std::vector<T>* data)
{_data       = data;}

/** \brief Sets an user defined metric, used as distance of reference and data vector.
*
* The user can define an own metric to measure the distance or distorian of the reference 
* vectors and the data vectors. 
* The metric has to be of the form
*
* T user_metric(const T&,const T&)
*
* Implicitly this defines a topology.
* The default value of the function pointer is NULL, this allows after an user defined metric
* was given to call the function with no parameter which effects that the default implemented
* metric is going to be used.
* \param *metric_to_use is function ptr to an user defined metric
*/

template < typename T > inline void NeuralGas<T>::setMetric(T (*metric_to_use)(const T& a,const T& b)=NULL )
{_metric_to_use=metric_to_use;}  

/** \brief standard is the pre-specified L2 euclidean metric
*
* \param x first vector
* \param y second vector
*/
template< typename T> T NeuralGas<T>::metric(const Vector<T>& x, const Vector<T>& y)
{
 T result;
 T value;
 Vector<T> z = x - y;
 
 result =  z[0] - z[0];
 
 for (int i=0; i < x.size(); i++)
 {
     value  = (z[i]*z[i]);
     result+=  z[i];
 }
    
 return sqrt(result); 
}  

/** \brief Sets the number of initial reference vectors.
*
* The function sets the initial number of reference vectors and initializes those
* with random values.
* \param num_of_ref_vec is the number of initial reference vectors
*/
template < typename T > void NeuralGas<T>::setRefVectors(const int& num_of_ref_vec)
{
  if (graphptr!=NULL)
      delete graphptr;
  
  graphptr  = new Graph<T>(_dimension);
  (*graphptr).initRandomGraph(num_of_ref_vec);  // creates a Graph object with given size of the vectors and number of ref vectors initilized with random values
}  

/** \brief Applies a given function to all neighboring nodes.
*
* The function applies a given function of the form
* 
* void (node<T>& n)
*
* to all neigboring nodes of the given node which is addresed by its index.
* \param home is the homie of the neighbors in question 
* \param func function that has to be applied
*/
template < typename T > void NeuralGas<T>::applyFunc2Neighbors(const int& home, void (*func)(Base_Node<T>& n))
{
  (*graphptr).applyFunc2Neighbors(home,func);
}  

/** \brief Applies a given function to all nodes.
*
* The function applies a given function of the form
* 
* void (node<T>& n)
*
* to all nodes of the underlying graph.
* \param func function that has to be applied
*/
template < typename T > void NeuralGas<T>::applyFunc2AllNodes(void (*func)(Base_Node<T>& n))
{
  (*graphptr).applyFunc2AllNodes(func);
}  

/** \brief Function returns the vectors dimension.
*/
template < typename T > int NeuralGas<T>::getDimension(void) const 
{return _dimension;}  

#endif
