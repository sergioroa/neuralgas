/** 
* \file ErrorTesting.h
* \author Manuel Noll
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef ERRORTESTING_H
#define ERRORTESTING_H

#include <vector>
#include <GrowingNeuralGas/GNGModul.h>

namespace neuralgas {

/** \class ErrorTesting
 *  \brief This class allows to compute the quantization error of NeuralGas Algorithms.
 *
 *   By giving a ptr of type GNGModul or of a thereof derived class to the class constructor
 *   the ErrorTesting class has full access to all members of that GNGModul object since
 *   the ErrorTesting class is defined as friend class of GNGModul.
 *   The class uses the fitted GNG network and the within the GNG object contained data 
 *   to compute the error either by an arbitrary number of data vectors beginning from the end
 *   or by a number of randomly chosen data vectors.
 *   As error metric the distance used within the GNG object is used as default. But it is
 *   also possible to define an own error metric and set it by the func setErrorMetric.
 *
 * \param _gngptr ptr to the given GNG object
 * \param _pastTimeSteps defines the number of data items that shall be used, correspoding to the number of returned error points
 * \param _random bool that says whether the data items are selected randomly, default is false
 * \param _errors a vector containing the number of errors
 * \param _error_metric_to_use is the user defined error metric,as default the metric of the given derived GNGModul object is used
 *
 */   
template < typename T, typename S > class ErrorTesting
{
 public:
    ErrorTesting();
 
    // cto taking a GNGModul ptr 
    ErrorTesting(GNGModul<T,S>*);
    // std dto
    ~ErrorTesting();

    void            setGNGObject(GNGModul<T,S>*);
    // returns the errors for given number of data items that are either selected at random or not
    std::vector<T>  getErrors(const unsigned int &,const bool&);
    //sets a user defined metric, used for error measuring 
    inline void     setErrorMetric(T (*)(const Vector<T>& a,const Vector<T>& b));  
            
 private:
    // returns the distance of the node and the datum vector
    T               getErrorDistance(const unsigned int&, const unsigned int&);
    // sets whether to choose the data randomly, default is false
    void            chooseRandomly(const bool&);
    // returns a random index within the range of the data size
    const unsigned int       getRandomIndex() const;
    // returns the shortest distance for a given datum
    const T         getShortestDistance(const unsigned int&);
    // calculates the errors 
    void            calcErrors();
    // sets the number of data items that shall be used
    void            setPastTimeSteps(const unsigned int&);
    // ptr to the given GNG object
    GNGModul<T,S> *  _gngptr;
    // defines the number of data items that shall be used, correspoding to the number of returned error points
    unsigned int     _pastTimeSteps;
    // bool that says whether the data items are selected randomly, default is false
    bool             _random;
    // returns the errors for given number of data items that are either selected at random or not    
    std::vector<T>   _errors;
    // is the user defined error metric
    T               (*_error_metric_to_use)(const Vector<T>&,const Vector<T>&); 

};


template < typename T, typename S > ErrorTesting<T,S>::ErrorTesting()
{
 _gngptr              = NULL;
 _pastTimeSteps       = 1;
 _random              = false;
 _errors.resize(_pastTimeSteps);
 _error_metric_to_use = NULL;

}

/** \brief cto taking a GNGModul ptr 
]*
*   \param gngptr a ptr to a GNGModul or a thereof derived class for which the error shall be computed
*/
template < typename T, typename S > ErrorTesting<T,S>::ErrorTesting(GNGModul<T,S>* gngptr)
{
 _gngptr              = gngptr;
 _pastTimeSteps       = 1;
 _random              = false;
 _errors.resize(_pastTimeSteps);
 _error_metric_to_use = NULL;
}

/** \brief std dto
*/
template < typename T, typename S > ErrorTesting<T,S>::~ErrorTesting()
{_gngptr  = NULL;}

//! \brief set the GNG object (algorithm) pointer to be tested
/*! 
  \param gngptr 
*/
template < typename T, typename S > void ErrorTesting<T,S>::setGNGObject(GNGModul<T,S>* gngptr)
{ _gngptr              = gngptr;}

/* \brief Sets the number of data items that shall be used
*
*  The number of items is defined for which the error calculation is done.
*  If the random flag is set then they are chosen randomly, otherwise the last
*  pastTimeSteps data items are used.
*
*  \param pastTimeSteps
*/
template < typename T, typename S > void ErrorTesting<T,S>::setPastTimeSteps(const unsigned int& pastTimeSteps)
{        
 if( pastTimeSteps <= _gngptr->size() ) //check whether the desired number of steps is larger than the available data
  _pastTimeSteps = pastTimeSteps;
 else
  _pastTimeSteps = _gngptr->size();
}

/** \brief sets whether to choose the data randomly, default is false
*/ 
template < typename T, typename S > void ErrorTesting<T,S>::chooseRandomly(const bool& random)
{_random=random;}

/** \brief returns the errors for given number of data items that are either selected at random or not
*
*
* \param steps number of data items that shall be taken to compute the error
* \param random determines whether the data is taken randomly or whether the last number of data items is used
*/
template < typename T, typename S > std::vector<T> ErrorTesting<T,S>::getErrors(const unsigned int& steps,const bool& random=false) 
{
 chooseRandomly(random);
 setPastTimeSteps(steps);
 calcErrors();
 return _errors;
}

/** \brief Sets an user defined metric, used as distance for error measuring.
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
* \param *errormetric_to_use is function ptr to an user defined metric
*/
template < typename T, typename S > inline void ErrorTesting<T,S>::setErrorMetric(T (*error_metric_to_use)(const Vector<T>& a,const Vector<T>& b)=NULL )
{_error_metric_to_use=error_metric_to_use;}  

// returns the distance of the node and the item vector
template<typename T,typename S> T ErrorTesting<T,S>::getErrorDistance(const unsigned int& node_index, const unsigned int& item_index)
{
     return _error_metric_to_use( (*_gngptr)[item_index], _gngptr->graphptr->operator[](node_index).weight);
}

// returns a random index within the range of the data size
template < typename T, typename S > const unsigned int ErrorTesting<T,S>::getRandomIndex() const 
{return (rand() % _gngptr->size()); }

// returns the shortest distance for a given datum
template < typename T, typename S > const T ErrorTesting<T,S>::getShortestDistance(const unsigned int& index)
{
    // init with zero
   T distance      = _gngptr->_zero;
   T best_distance = _gngptr->_zero; 

   // best_distance set to "infinity"
   best_distance = std::numeric_limits<T>::max();

   if (_error_metric_to_use == NULL)
   {
       for (unsigned int j = 0; j < _gngptr->graphptr->size(); j++)
       {
	distance = dynamic_cast<GNGModulGraph<T,S>* >(_gngptr->graphptr)->getDistance((*_gngptr)[index],j);
    
        if (distance < best_distance)
        {
         best_distance              =       distance;
        }       
       }
   }
   else
   {
       for (unsigned int j = 0; j < _gngptr->graphptr->size(); j++)
       {
        distance = getErrorDistance(index,j);
    
        if (distance < best_distance)
        {
         best_distance              =       distance;
        }       
       }
   }
   return best_distance;
}
   // calculates the errors 

template < typename T, typename S > void ErrorTesting<T,S>::calcErrors()
{
 _errors.resize(_pastTimeSteps);

 if (_random)
 {
  for(unsigned int i = 0; i < _pastTimeSteps; i++)
          _errors[i]=getShortestDistance(getRandomIndex());         
 }
 else
 {
  for(unsigned int i = 0; i < _pastTimeSteps; i++)
          _errors[i]=getShortestDistance(_gngptr->size() - _pastTimeSteps + i);         
 }
}

} // namespace neuralgas

#endif
