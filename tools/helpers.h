/** 
* \file helpers.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef HELPERS_H
#define HELPERS_H

#include <Graphs/Vector.h>
#include <fstream>
#include <iostream>

namespace neuralgas {

template < typename T , typename S > struct Base_Node;

template<typename T>
bool readDataText(const char* filename, std::vector < Vector<T>* >* _data)
{
    std::ifstream myfile (filename);
    assert ( _data != NULL );

    if (myfile.is_open())
    {
	char comment;
	int size, dim;

	myfile >> comment;
	myfile >> size;
	myfile >> comment;
	myfile >> dim;
	for (int i=0; i < size; i++)
	{
	    Vector<T>* item = new Vector<double>(dim);
	    for (int j=0; j<dim; j++)
		myfile >> item->at(j);
	    _data->push_back (item);
	}
	
        myfile.close();
	return true;
    }
    return false;
}

template<typename T>
bool readData(const char* filename, std::vector< Vector<T>* >* _data)
{
    std::ifstream myfile (filename, std::ios::in | std::ios::binary);
    assert ( _data != NULL );
    
    if (myfile.is_open ())
    {
	int size, dim;
	myfile.read ((char *)&size, sizeof (size));
	myfile.read ((char *)&dim, sizeof (dim));

	for (int i=0; i<size; i++)
	{
	    Vector<T>* item = new Vector<T>(dim);
	    for (int j=0; j<dim; j++)
		myfile.read ((char *)&item->at(j), sizeof (item->at(j)));
	    _data->push_back (item);
	}
        myfile.close();
	return true;
    }
    return false;
}

template<typename T>
bool saveData(const char* filename, std::vector < Vector<T>* >* _data)
{
    std::ofstream myfile (filename, std::ios::out | std::ios::binary);
    assert ( _data->size() );
    
    if (myfile.is_open ())
    {
	int size = _data->size();
	int dim = _data->at(0)->size();
	
	myfile.write ((const char*)&size, sizeof (size));
	myfile.write ((const char*)&dim, sizeof (dim));
	for (int i=0; i<size; i++)
	    for (int j=0; j<dim; j++)
		myfile.write((const char*)&_data->at(i)->at(j), sizeof(_data->at(i)->at(j)));
    }
    return false;

    
}

template<typename T>
bool saveDataText(const char* filename, std::vector < Vector<T>* >* _data)
{
    std::ofstream myfile (filename);
    assert (_data->size());
   
    if (myfile.is_open())
    {
	int size = _data->size();
	int dim = _data->at(0)->size();
	char comment = '#';
	myfile << comment << " ";
	myfile << size << std::endl; 
	myfile << comment << " ";
	myfile << dim << std::endl;
	for(int i = 0; i < size; i++)
	{      
	    for(int j=0; j < dim; j++)
		myfile << _data->at(i)->at(j) << " ";
	    if (i != size-1)
		myfile << std::endl;
	}
	
	myfile.close();
	return true;
    }
    return false;
}



template<typename T>
void showData(std::vector < Vector<T>* >* _data)
{
    for(unsigned int i=0; i < _data->size(); i++)
    {
	for (unsigned int j=0; j < _data->at(i)->size(); j++)
	    std::cout << (*(*_data)[i])[j] << " ";
	std::cout << std::endl;
    }
}

template<typename T, typename S>
void showNodes (std::vector < Base_Node<T, S>* >* _nodes)
{
    for(unsigned int i=0; i < _nodes->size(); i++)
    {
	for (unsigned int j=0; j < _nodes->at(i)->weight.size(); j++)
	    std::cout << (*(*_nodes)[i]).weight[j] << " ";
	std::cout << std::endl;
    }
}


template<typename T, typename S>
bool readNodesText(const char* filename, std::vector < Base_Node<T, S>* >* _nodes)
{
    assert ( _nodes != NULL );
    std::ifstream myfile (filename);
    
    if (myfile.is_open())
    {
	char comment;
	int size, dim;

	myfile >> comment;
	myfile >> size;
	myfile >> comment;
	myfile >> dim;
	for (int i=0; i<size; i++)
	{
	    Base_Node<T, S>* node = new Base_Node<T, S>;
	    node->weight.resize (dim);
	    for (int j=0; j<dim; j++)
		myfile >> node->weight[j];
	    _nodes->push_back (node);
	}

	
        myfile.close();
	return true;
    }
    return false;

}


template<typename T, typename S>
bool readNodes(const char* filename, std::vector < Base_Node<T, S>* >* _nodes)
{
    assert ( _nodes != NULL );
    std::ifstream myfile (filename, std::ios::in | std::ios::binary );
    
    if (myfile.is_open())
    {
	int size, dim;
	myfile.read ((char *)&size, sizeof (size));
	myfile.read ((char *)&dim, sizeof (dim));

	for (int i=0; i<size; i++)
	{
	    Base_Node<T, S>* node = new Base_Node<T, S>;
	    node->weight.resize (dim);
	    for (int j=0; j<dim; j++)
		myfile.read ((char *)&node->weight[j], sizeof (node->weight[j]));
	    _nodes->push_back (node);
	}
	myfile.close();
	return true;
    }
    return false;
}

template<typename T, typename S>
bool saveNodes(const char* filename, std::vector < Base_Node<T, S>* >* _nodes)
{
    std::ofstream myfile (filename, std::ios::out | std::ios::binary);
    assert (_nodes->size());

    if (myfile.is_open())
    {
	int size = _nodes->size();
	int dim = _nodes->at(0)->weight.size();
	myfile.write ((const char*)&size, sizeof (size));
	myfile.write ((const char*)&dim, sizeof (dim));
	for (int i=0; i<size; i++)
	    for(int j=0; j < dim; j++)
		myfile.write ((const char*)&_nodes->at(i)->weight[j], sizeof(_nodes->at(i)->weight[j]));

	myfile.close ();
	return true;
	    
    }
    return false;
}


template<typename T, typename S>
bool saveNodesText(const char* filename, std::vector < Base_Node<T, S>* >* _nodes)
{
    std::ofstream myfile (filename);
    assert (_nodes->size());
    
    if (myfile.is_open())
    {
	char comment = '#';
	unsigned int size = _nodes->size();
	unsigned int dim = _nodes->at(0)->weight.size();
	myfile << comment << " ";
	myfile << size << std::endl;
	myfile << comment << " ";
	myfile << dim << std::endl;
	for(unsigned int i = 0; i < size; i++)
	{      
	    for(unsigned int j=0; j < dim; j++)
		myfile << _nodes->at(i)->weight[j] << " ";
	    if (i != size-1)
		myfile << std::endl;
	}
       
	myfile.close();
	return true;
    }
    return false;

}


} // namespace neuralgas

#endif // HELPERS_H
