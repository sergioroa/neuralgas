/*
 *   This file is part of NeuralGas.
 *
 *   NeuralGas is free software: you can redistribute it and/or modify it
 *   under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   NeuralGas is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with NeuralGas.  If not, see <http://www.gnu.org/licenses/>.
 */

/** 
* \file NeuralGas.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  Copyright(c) 2011 Sergio Roa
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
#include <tools/math_helpers.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <sstream>

/** \namespace neuralgas
    \brief Classes that follow the idea of Hebbian Learning proposed by the original Neural Gas algorithm
*/
namespace neuralgas {

/// NUM_PARAM is the maximal number of used parameters, has to be changed if a larger number of
/// parameters has to be used in a derived class
#define NUM_PARAM 9

/**
 * \enum _sampling_mode
 * \brief for using different mechanisms for sampling learning data
 */
enum _sampling_mode {sequential, /**< sample instances sequentially */
		     randomly /**< sample instances randomly */
};

/**
 * \enum _stopping_criterion
 * \brief for using different mechanisms for stopping the algorithms
 */
enum _stopping_criterion {epochs, /**< nr of training epochs */
                          global_error, /**< a global error measure */
			  stability /**< a network stability measure */
};

typedef boost::iostreams::tee_device<std::stringstream, std::ofstream> TeeDev;
typedef boost::iostreams::stream<TeeDev> TeeStream;

/** \class NeuralGas
 * \brief This class defines the basic operations for the Neural Gas Algorithm and is intended as a super class for inheritance.
 *
 * The Neural Gas Algorithm as proposed by Martinetz and Schulten, "A Neural Gas Network Learns Topologies"
 * served to build a general algorithmic framework.
 * The underlying graph structure is the Base_Graph which is an abstract graph class
 * forcing the NeuralGas class to be an abstract class as well.
 * \param *graphptr a pointer to the underlying data structure representing the graph
 * \param _zero is the zero element of the unknown datatype T
 * \param _dimension is dimension of the input data and therefore the dimension of the graph weight vectors
 * \param *_data is a pointer to the given input data; a pointer is used, since for large data a duplicate set could be to memory consuming
 * \param *_metric_to_use is a function pointer which is NULL if not a particular metric is given via setMetric
 */
template < typename T, typename S > class NeuralGas
{
    friend class boost::serialization::access;
public:
    typedef T (NeuralGas::*Metric)(const Vector<T>&,const Vector<T>&) const;
    //cto with size of dimension as input
    NeuralGas(const unsigned int&);
    // copy cto
    NeuralGas(const NeuralGas&);
    // default cto for deserialization
    NeuralGas ();
    //std dto
    ~NeuralGas(void);
    //sets the data to be processed
    inline void             setData(std::vector< Vector<T>* >*);
	
    //adds a single datum
    inline void             addData(Vector<T>*);
    //adds an arbitrary number of data
    inline virtual void     addData(std::vector< Vector<T>* >*);
    // returns the number of data items currently stored
    inline unsigned int     size(void) const;
    // sets the number of inital reference vectors
    virtual void            setRefVectors(const unsigned int&)=0;
    //sets a user defined metric, used as distance of reference and data vector 
    //inline void             setMetric(T (*)(const Vector<T>& a,const Vector<T>& b));
    inline void             setMetric(Metric); 
    //assigns int depending functions to the parameters 
    void                    setFuncArray(T (*)(const unsigned int&),const unsigned int&);
    //determines the maximal value within the given data set
    const T                 maxValue() const;
    //determines the minimal value within the given data set
    const T                 minValue() const;
    //determines the maximal values for each dim within the given data set
    Vector<T>               maxValues() const;
    //determines the minimal values for each dim within the given data set
    Vector<T>               minValues() const;
    // saves the nodes weight in a file
    void                    save(const char*, bool t = false);
    // loads nodes
    void                    setNodes( std::vector < Base_Node<T, S>* >* nodes);
    // sets the sampling mode for learning instances
    virtual void            setSamplingMode (unsigned int);
    //sets a user defined stopping criterion
    virtual void            setStoppingCriterion (unsigned int);
    
    //an abstract run func
    virtual void            run()=0;

    //output redirection
    void                    redirectOutput (std::string file);
protected:
    // pre-specified metric is the standard L2 euclidean metric
    virtual T               metric(const Vector<T>&, const Vector<T>&) const;
    // applies a given function to all neighbors
    inline void             applyFunc2Neighbors(const unsigned int&, void (*)(Base_Node<T,S>* n,const float&),const float&); 
    // applies a given function to all nodes
    inline void             applyFunc2AllNodes(void (*)(Base_Node<T,S>* n,const float&),const float&); 
    // returns the dimension of the data vectors
    inline unsigned int     getDimension(void) const; 
    // ptr to the underlying graph
    Base_Graph<T,S>*        graphptr;          
    // returns a reference to the indexed element
    Vector<T>&              operator[](const unsigned int&);
    // returns a const reference to the indexed element
    const Vector<T>&        operator[](const unsigned int&) const;
    // zero element of unknown datatype T, is set in constructor
    T                       _zero;
    // functions for assigning an integer depending function to the parameters
    T (*_funcArray[NUM_PARAM])(const unsigned int&);
    // parameters for the algorithm
    float   params[NUM_PARAM];
    // sampling mode of learning instances
    unsigned int sampling_mode;
    // stopping criterion
    unsigned int stopping_criterion;
    /// streams for output redirection
    std::ostream* out;
    std::ofstream* logout;
    TeeDev* tdev;
    TeeStream* tout;  
    std::stringstream tss;
    //dimension of the vectors
    unsigned int _dimension; 
    // ptr to the input data
    std::vector< Vector<T>* >* _data; 
    //user specified metric
    //T                       (*_metric_to_use)(const Vector<T>&,const Vector<T>&);
    Metric _metric_to_use;


private:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);


};  

/** \brief Cto with the dimension of the vectors. 
*
* Cto has as input the to-be-used dimension for the vectors, leading to an equal
* dimension for the weight vectors of the Node and the Edge struct. 
* The dimension has to be the dimension of the input data.
* The dimension cannot be changed after the creation of the object.
* \param dimension is dimension of the input data and therefore the dimension of the graph weight vectors
*/
template < typename T, typename S > NeuralGas<T,S>::NeuralGas(const unsigned int& dimension)
  : out (&std::cout)
{
  T x = T();
  _zero           = x - x;
  _dimension      = dimension;
  // _data           = new std::vector< Vector<T>* >;
  _data           = NULL;
  graphptr        = NULL;
  _metric_to_use  = NULL;   // no user defined metric is used, std metric is used
  for (int i=0;i < NUM_PARAM; i++)
      _funcArray[i] = NULL;
  sampling_mode = sequential;
  stopping_criterion = epochs;

  logout = 0;
  tdev = 0;
  tout = 0;
}

/** \brief Copy constructor
 */
template < typename T, typename S > NeuralGas<T,S>::NeuralGas (const NeuralGas& n) :
  _zero (n._zero),
  sampling_mode (n.sampling_mode),
  stopping_criterion (n.stopping_criterion),
  out (&std::cout),
  logout (0),
  tdev (0),
  tout (0),
  _dimension (n._dimension),
  _data (0),
  _metric_to_use (n._metric_to_use)
{
}

/** \brief cto only used for deserialization
 */
template < typename T, typename S > NeuralGas<T,S>::NeuralGas () :
  graphptr (0),
  _zero (0),
  sampling_mode (sequential),
  stopping_criterion (epochs),
  out (&std::cout),
  logout (0),
  tdev (0),
  tout (0),
  _dimension (0),
  _data (0),
  _metric_to_use (0)
{
}

/** \brief Dto frees memory by deleting the underlying graph data structure.
*/
template < typename T, typename S > NeuralGas<T,S>::~NeuralGas(void)
{
 //if (graphptr!=NULL) 
   // delete graphptr;
  graphptr=NULL;
  
  if (_data!=NULL)
  {
    for(unsigned int i=0; i < _data->size(); i++)
      delete (*_data)[i];
    _data->clear();
    delete _data;
    
  }
  if (tout)
  {
    tout->close();
    delete tout;
    if (tdev)
	    delete tdev;
    if (logout)
	    delete logout;
  }
}

/** \brief Assigns int depending functions to the parameters 
*
*   \param func is a func ptr pointing to an int depending function reflecting a variable parameter
*   \param index is the index of the parameter to be assigned by a function
*/
template<typename T,typename S> 
void NeuralGas<T,S>::setFuncArray(T (*func)(const unsigned int& time),const unsigned int& index)
{
 _funcArray[index]=func;
}

/** \brief Determines the maximal value within the given data set
 *
 */
template<typename T,typename S> const T NeuralGas<T,S>::maxValue() const
{
    return neuralgas::maxValue (_data);
}
/** \brief Determines the minimal value within the given data set
 *
 */
template<typename T,typename S> const T NeuralGas<T,S>::minValue() const
{
    return neuralgas::minValue (_data);
}
/** \brief Determines the minimal values in each dim within the given data set
 *
 */
template<typename T,typename S> Vector<T> NeuralGas<T,S>::minValues() const
{
    return neuralgas::minValues (_data);
}
/** \brief Determines the maximal values in each dim within the given data set
 *
 */
template<typename T,typename S> Vector<T> NeuralGas<T,S>::maxValues() const
{
    return neuralgas::maxValues (_data);
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
* \param text save the file in text format if true
*/


template<typename T,typename S> void NeuralGas<T,S>::save(const char* filename, bool text)
{
  graphptr->save(filename, text);
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
* A copy of the data pointer is done. Thus, please delete properly the memory if needed.
* \param *data is the given data
*/
template < typename T, typename S > inline void NeuralGas<T,S>::setData(std::vector< Vector<T>* >* data)
{
	// _data       = data;
  if (_data != NULL) {
    for(unsigned int i = 0; i < _data->size(); i++)
      delete  _data->at(i);                                // delete ptrs to the nodes
    _data->clear();
  }
  else
    _data = new std::vector< Vector<T>* >;

  for (unsigned int i=0; i < data->size(); i++) {
    Vector<T> *new_item = new Vector<T>(*data->at(i));
    _data->push_back (new_item);
  }

}

/** \brief Adds a single item
*
*   \param to_add the item to add
*/
template < typename T, typename S > inline void NeuralGas<T,S>::addData(Vector<T>* to_add)
{_data->push_back(to_add);}

/** \brief adds an arbitrary number of data
*
*   \param to_add the data to add
*/
template < typename T, typename S > inline void NeuralGas<T,S>::addData(std::vector< Vector<T>* >* to_add)
{
  for(unsigned int i=0;i < to_add->size(); i++)
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
{
  _metric_to_use=metric_to_use;
  graphptr->setMetric (metric_to_use);

}  


/** \brief standard is the pre-specified L2 euclidean metric
*
* \param x first vector
* \param y second vector
*/
template< typename T, typename S> T NeuralGas<T,S>::metric(const Vector<T>& x, const Vector<T>& y) const
{
 if (_metric_to_use==NULL) // is a non-standard metric set ?
 {
     return euclidean<T,S> (x, y);
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
/*template < typename T, typename S > void NeuralGas<T,S>::setRefVectors(const unsigned int& num_of_ref_vec,const int& max_value)
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
* \param index is the nodex index 
* \param func function that has to be applied
* \param value value to be applied
*/
template < typename T, typename S > void NeuralGas<T,S>::applyFunc2Neighbors(const unsigned int& index, void (*func)(Base_Node<T,S>* n,const float&),const float& value)
{
  graphptr->applyFunc2Neighbors(index,func,value);
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
template < typename T, typename S > unsigned int NeuralGas<T,S>::getDimension(void) const 
{return _dimension;}  

/** returns the number of data items currently stored
*/
template < typename T, typename S > unsigned int NeuralGas<T,S>::size(void) const
{return _data->size();}

/** operator[] returns a reference to the indexed data element element
*
* \param index
*/
template <typename T, typename S > Vector<T>& NeuralGas<T,S>::operator[](const unsigned int& index)
{return *(*_data)[index];}

/** const operator[] returns a const reference to the indexed data element element
*
* \param index
*/
template <typename T, typename S > const Vector<T>& NeuralGas<T,S>::operator[](const unsigned int& index) const
{return *(*_data)[index];}

/** \brief Loads the nodes weight from a file
    \param filename
    \param text if text file then true
 */
template<typename T, typename S> void NeuralGas<T,S>::setNodes( std::vector < Base_Node<T, S>* >* nodes)
{
  graphptr->setNodes (nodes);
}

/** \brief Redirects output to a file */
template<typename T, typename S>
void NeuralGas<T,S>::redirectOutput (std::string logname)
{
  logout = new std::ofstream(logname.c_str());
  if(!logout->is_open())
    std::cerr << "can't open log file " << logname << std::endl;
  tdev = new TeeDev(tss, *logout);
  tout = new TeeStream(*tdev);
  out = tout;
  std::cout << "writing to log file " << logname << std::endl;
  
  
}

template<typename T, typename S>
template<class Archive>
void 
NeuralGas<T,S>::serialize(Archive & ar, const unsigned int /* file_version */)
{
  ar & BOOST_SERIALIZATION_NVP(sampling_mode);
  ar & BOOST_SERIALIZATION_NVP(stopping_criterion);
  ar & BOOST_SERIALIZATION_NVP(_dimension);
}


} // namespace neuralgas

#endif
