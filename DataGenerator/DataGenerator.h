/** 
* \file DataGenerator.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  Copyright(c) 2010 Manuel Noll - All rights reserved
*  \version 1.0
*  \date    2010
*/

#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H

#include <Graphs/Vector.h>
#include <tools/helpers.h>
#include <fstream>

namespace neuralgas {

/** \class DataGenerator
 *  \brief The class is an abstract class that provides the basic operations for data generation
 *   for the NeuralGas algorithms.
 *
 *   The class is an abstract class that provides the data generation process for the NeuralGas
 *   Algorithms. It is intended to be derived by classes that either generate the data
 *   by their own e.g. via a mathematical funcition or represent a file format
 *   that allows reading the data from a file or to stream in data online successively.
 *
 *   \param _data contains the generated data
 */
template<typename T> class DataGenerator
{
 public:
        // std cto with dimension of the data
                                                DataGenerator(const int&);
        // std dto
                                                ~DataGenerator();
         // erases the data
        void                                    reset();
        // returns the next item
        virtual   Vector<T>*                    next()=0;
        // returns the data
        std::vector< Vector<T>* >*              getData();
        // generates a given number of data
        virtual void                            generate(const int&);
	// saves dataset to a text file
	bool                                    save(const char* filename, bool t = false);

 protected:
	// generate an item
        virtual   Vector<T>*                    generate()=0;
        // dimension of the data
        int                                    _dim;
        // contains the generated data
        std::vector< Vector<T>* >*             _data;
        // zero element of unknown datatype T, is set in constructor
        T                       _zero;   
};

/** \brief std cto with the dimension of the data
*
*  \param dim is the dimension of the data
*/
template<typename T> DataGenerator<T>::DataGenerator(const int& dim)
{
 T x = T ();
 _zero            = x - x;
 _data            =          new std::vector< Vector<T>* >;                             
 _dim             =          dim;
}

/** \brief std dto
*/
template<typename T> DataGenerator<T>::~DataGenerator()
{

 if (_data != NULL)
 {
	 for(unsigned int i=0; i < _data->size(); i++)
		 delete (*_data)[i];
	 _data->clear ();
	 delete _data;
 }

}


/** \brief Erases the data generated so far
*/
template<typename T> void DataGenerator<T>::reset()
{
 for(unsigned int i=0; i < _data->size(); i++)
         delete (*_data)[i];
 _data->clear();
 _data=NULL;
}

/** \brief Returns the next item
*
*   The function is intended for implementing successive item returning,
*   either in the sense that step by step a new item is generated and returned
*   or is read in from a file line by line or waits for instreaming online data.
*   It may be necessary to implement the func in such a way that first all the data
*   within the vector _data is returned and afterwards newly generated data is returned.
*   Dependent on the requirements this newly generated data may be added to the data vector.
*/
template<typename T> Vector<T>* DataGenerator<T>::next(){}

/** \brief Returns the data
*/
template<typename T> std::vector< Vector<T>* >* DataGenerator<T>::getData()
{return _data;}

/** \brief Generates a given number of data
*
* \param number is the number of data items that shall be generated
*/
template<typename T> void DataGenerator<T>::generate(const int& number)
{
 for(int i=0; i < number; i++)
 {
  _data->push_back(generate());
 }
}

/** \brief Generates a single item
 *
 *   This function represents the heart of the class. It defines how to "generate"
 *   item, meaning whether it should be read from a file, generated via a function
 *   or online streamed.
 *   
 */
template<typename T> Vector<T>* DataGenerator<T>::generate()=0;

/** \brief Saves the dataset to a text file
*
*   The functions saves each vector in a dataset to a file where
*   for n elements Xi of d dimensional vectors the format looks as follows
*   X11 X12 X13 ... X1d
*   .
*   .
*   .
*   Xn1 Xn2 Xn3 ... Xnd
*
* which means that they are separated by spaces.
*
* \param filename is the name of the file where to store the data to
*/
template<typename T> bool DataGenerator<T>::save(const char* filename, bool text)
{
	if (text)
		return neuralgas::saveDataText (filename, _data);
	else
		return neuralgas::saveData (filename, _data);
}

} // namespace neuralgas

#endif
