/** 
* \class NeuralGas
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

// TODO      reconsidering the handling of the data as pointers, in particular memory freeing

#ifndef NEURALGAS_H
#define NEURALGAS_H

#include <math.h>
#include <vector>
#include <limits>
#include <Graphs/Base_Graph.h>

/// Classes that follow the idea of Hebbian Learning proposed by the original Neural Gas algorithm
namespace neuralgas {

/// NUM_PARAM is the maximal number of used parameters, has to be changed if a larger number of
/// parameters has to be used in a derived class
#define NUM_PARAM 9

/**
 * \enum sampling_mode
 * \brief for using different mechanisms for sampling learning data
 */
enum _sampling_mode {sequential, /**< sample instances sequentially */
		     randomly /**< sample instances randomly */
};

/**
 * \enum stopping_criterion
 * \brief for using different mechanisms for stopping the algorithms
 */
enum _stopping_criterion {epochs, /**< nr of training epochs */
                          global_error /**< a global error measure */
};


/** \brief Class declares the basic operations for the Neural Gas Algorithm and is intended as a super class for inheritance.
*
* The Neural Gas Algortihm as proposed by Martinetz and Schulten, "A Neural Gas Network Learns Topologies"
* served to build a general algorithmic framework.
* The underlying graph structure is the Base_Graph which is an abstract graph class
* forcing the NeuralGas class to be an abstract class as well.
* \param *graphptr a pointer to the underlying data structure representing the graph
* \param _zero is the zero element of the unknown datatype T
* \param _unit is the unit element
* \param _dimension is dimension of the input data and therefore the dimension of the graph weight vectors
* \param *_data is a pointer to the given input data; a pointer is used, since for large data a duplicate set could be to memory consuming
* \param *_metric_to_use is a function pointer which is NULL if not a particular metric is given via setMetric
*/

template < typename T, typename S > class NeuralGas
{
public:
    typedef T (NeuralGas::*Metric)(const Vector<T>&,const Vector<T>&) const;
    //cto with size of dimension as input
    NeuralGas(const int&);         
    //std dto
    ~NeuralGas(void);               
    //sets the data to be processed
    inline void             setData(std::vector< Vector<T>* >*);  
    //adds a single datum
    inline void             addData(Vector<T>*);
    //adds an arbitrary number of data
    inline void             addData(std::vector< Vector<T>* >*);
    // sets the number of inital reference vectors
    virtual void            setRefVectors(const int&,const int&)=0;    
    //sets a user defined metric, used as distance of reference and data vector 
    //inline void             setMetric(T (*)(const Vector<T>& a,const Vector<T>& b));
    inline void             setMetric(Metric); 
    //assigns int depending functions to the parameters 
    void                    setFuncArray(float (*)(const int&),const int&); 
    //determines the maximal value within the given data set
    const int               maxRandomValue() const; 
    // saves the nodes weight in a file
    void                    save(const char*);
    // sets the sampling mode for learning instances
    virtual void            setSamplingMode (unsigned int);
    //sets a user defined stopping criterion
    virtual void            setStoppingCriterion (unsigned int);
    
    //an abstract run func
    virtual void            run()=0;

  protected:
    // pre-specified metric is the standard L2 euclidean metric
    virtual T               metric(const Vector<T>&, const Vector<T>&) const;
    // applies a given function to all neighbors
    inline void             applyFunc2Neighbors(const int&, void (*)(Base_Node<T,S>* n,const float&),const float&); 
    // applies a given function to all nodes
    inline void             applyFunc2AllNodes(void (*)(Base_Node<T,S>* n,const float&),const float&); 
    // returns the dimension of the data vectors
    inline int              getDimension(void) const; 
    // ptr to the underlying graph
    Base_Graph<T,S>*        graphptr;          
    // returns the number of data items currently stored
    inline int              size(void) const;
    // returns a reference to the indexed element
    Vector<T>&              operator[](const int&);
    // returns a const reference to the indexed element
    const Vector<T>&        operator[](const int&) const;
    // zero element of unknown datatype T, is set in constructor
    T                       _zero;
    // functions for assigning an integer depending function to the parameters
    float   (*_funcArray[NUM_PARAM])(const int&);
    // parameters for the algorithm
    float   params[NUM_PARAM];
    // sampling mode of learning instances
    unsigned int sampling_mode;
    // stopping criterion
    unsigned int stopping_criterion;
 
  private:
    //dimension of the vectors
    int                     _dimension; 
    // ptr to the input data
    std::vector< Vector<T>* >* _data; 
    //user specified metric
    //T                       (*_metric_to_use)(const Vector<T>&,const Vector<T>&);
    Metric _metric_to_use;
    
      
};  

/** \brief Assigns int depending functions to the parameters 
*
*   \param func is a func ptr pointing to an int depending function reflecting a variable parameter
*   \param index is the index of the parameter to be assigned by a function
*/
template<typename T,typename S> 
void NeuralGas<T,S>::setFuncArray(float (*func)(const int& time),const int& index)
{
 _funcArray[index]=func;
}

/** \brief Determines the maximal value within the given data set
*
*/
template<typename T,typename S> const int NeuralGas<T,S>::maxRandomValue() const
{
 T max_data_value = _zero;
 for ( unsigned int i = 0; i < _data->size(); i++)
  for ( unsigned int j = 0; j < (*_data)[0]->size(); j++)
      if (_data->operator[](i)->operator[](j) > max_data_value ) 
         max_data_value = abs( int(ceil(_data->operator[](i)->operator[](j))) );
 return int(max_data_value);
} 

/** \brief Cto with the dimension of the vectors. 
*
* Cto has as input the to-be-used dimension for the vectors, leading to an equal
* dimension for the weight vectors of the Node and the Edge struct. 
* The dimension has to be the dimension of the input data.
* The dimension cannot be changed after the creation of the object.
* \param dimension is dimension of the input data and therefore the dimension of the graph weight vectors
*/
template < typename T, typename S > NeuralGas<T,S>::NeuralGas(const int& dimension)
{
  T x = T();
  _zero           = x - x;
  _dimension      = dimension;
  _data           = new std::vector< Vector<T>* >;
  _data           = NULL;
  graphptr        = NULL;
  _metric_to_use  = NULL;   // no user defined metric is used, std metric is used
  for (int i=0;i < NUM_PARAM; i++)
      _funcArray[i] = NULL;
  sampling_mode = sequential;
  stopping_criterion = epochs;
}

/** \brief Dto frees memory by deleting the underlying graph data structure.
*/
template < typename T, typename S > NeuralGas<T,S>::~NeuralGas(void)
{
 //if (graphptr!=NULL) 
   // delete graphptr;
 graphptr=NULL;

 if (_data!=NULL)
    delete _data;

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
*
* \param filename is the name of the file where to store the data to
*/


template<typename T,typename S> void NeuralGas<T,S>::save(const char* filename)
{
 graphptr->save(filename);
}


/**
 * \brief Sets the sampling mode for selecting learning instances
 * \param mode sampling mode from \p _sampling_mode enum variables
 */
template<typename T,typename S> void NeuralGas<T,S>::setSamplingMode(unsigned int mode)
{
  sampling_mode = mode;
}

/**
 * \brief Sets the stopping criterion of the algorithm
 * \param criterion criterion from \p _stopping_criterion enum variab
 */
template<typename T,typename S> void NeuralGas<T,S>::setStoppingCriterion(unsigned int criterion)
{
  stopping_criterion = criterion;
}


/** \brief Sets the data that has to be processed in the next algorithmic run.
*
* The function sets the data that has to be processed in the next algorithmic run.
* If a run already took place and the function setRefVectors is not called until the next run
* then the underlying graph persists and will fit to this newly added data, allowing
* a succesiv adaption of the graph.
* \param *data is the given data
*/
template < typename T, typename S > inline void NeuralGas<T,S>::setData(std::vector< Vector<T>* >* data)
{
 _data       = data;
}

/** \brief Adds a single datum
*
*   \param to_add the datum to add
*/
template < typename T, typename S > inline void NeuralGas<T,S>::addData(Vector<T>* to_add)
{_data->push_back(to_add);}

/** \brief adds an arbitrary number of data
*
*   \param to_add the data to add
*/
template < typename T, typename S > inline void NeuralGas<T,S>::addData(std::vector< Vector<T>* >* to_add)
{
 for(int i=0;i < to_add.size(); i++)
         _data->push_back( (*to_add)[i]);
}


/** \brief Sets an user defined metric, used as distance of reference and data vector.
*
* The user can define an own metric to measure the distance or distorian of the reference 
* vectors and the data vectors. 
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

//template < typename T, typename S > inline void NeuralGas<T,S>::setMetric(T (*metric_to_use)(const Vector<T>& a,const Vector<T>& b)=NULL )
//{_metric_to_use=metric_to_use;}  

template < typename T, typename S > inline void NeuralGas<T,S>::setMetric(Metric metric_to_use=NULL )
{_metric_to_use=metric_to_use;}  


/** \brief standard is the pre-specified L2 euclidean metric
*
* \param x first vector
* \param y second vector
*/
template< typename T, typename S> T NeuralGas<T,S>::metric(const Vector<T>& x, const Vector<T>& y) const
{
 if (_metric_to_use==NULL) // is a non-standard metric set ?
 {
     T result;
     T value;
     Vector<T> z = x - y;
     
     //result =  _zero;
     result = 0;
     int tsize = z.size();
     
     for (int i=0; i < tsize; i++)
     {
         value  = (z[i]*z[i]);
         result+=  value;
     }
     return T(sqrt(result));
 }
 else
	 return (this->*_metric_to_use)(x,y);   // use the non-standard user defined metric
}  

/** \brief Sets the number of initial reference vectors.
*
* The function sets the initial number of reference vectors and initializes those
* with random values.
* \param num_of_ref_vec is the number of initial reference vectors
*/
/*template < typename T, typename S > void NeuralGas<T,S>::setRefVectors(const int& num_of_ref_vec,const int& max_value)
{
  if (graphptr!=NULL)
      delete graphptr;

  graphptr  = new Base_Graph<T,S>(_dimension);
  graphptr->initRandomGraph(num_of_ref_vec,max_value);  // creates a Graph object with given size of the vectors and number of ref vectors initilized with random values
} */ 

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
template < typename T, typename S > void NeuralGas<T,S>::applyFunc2Neighbors(const int& home, void (*func)(Base_Node<T,S>* n,const float&),const float& value)
{
  graphptr->applyFunc2Neighbors(home,func,value);
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
template < typename T, typename S > void NeuralGas<T,S>::applyFunc2AllNodes(void (*func)(Base_Node<T,S>* n,const float&),const float& value)
{
  graphptr->applyFunc2AllNodes(func,value);
}  

/** \brief Function returns the vectors dimension.
*/
template < typename T, typename S > int NeuralGas<T,S>::getDimension(void) const 
{return _dimension;}  

/** returns the number of data items currently stored
*/
template < typename T, typename S > int NeuralGas<T,S>::size(void) const
{return _data->size();}

/** operator[] returns a reference to the indexed data element element
*
* \param index
*/
template <typename T, typename S > Vector<T>& NeuralGas<T,S>::operator[](const int& index)
{return *(*_data)[index];}

/** const operator[] returns a const reference to the indexed data element element
*
* \param index
*/
template <typename T, typename S > const Vector<T>& NeuralGas<T,S>::operator[](const int& index) const
{return *(*_data)[index];}

} // namespace neuralgas

#endif
