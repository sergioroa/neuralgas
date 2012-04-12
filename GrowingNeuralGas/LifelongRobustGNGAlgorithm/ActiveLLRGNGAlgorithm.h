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
* \file LLRGNGAlgorithm.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef ACTIVELLRGNGALGORITHM_H
#define ACTIVELLRGNGALGORITHM_H

#include <GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGAlgorithm.h>

namespace neuralgas {

/** \class ActiveLLRGNGAlgorithm
 *  \brief A version of LLRGNGAlgorithm that actively selects samples for learning
*/
template<typename T, typename S>
class ActiveLLRGNGAlgorithm : virtual public LLRGNGAlgorithm<T,S>
{
public:
	// cto class initialization
	ActiveLLRGNGAlgorithm (const unsigned int& dim, const unsigned int& window = 81) :
		LLRGNGAlgorithm<T,S> (dim, window) { }
	// std dto
	~ActiveLLRGNGAlgorithm ();
	// run the algorithm
	void run();
	// add new data to the current data set
	inline virtual void addData (std::vector< Vector<T>* >*);
	// initialize data for visualization
	void initializeDataVisualization ();
protected:
	/// data indices used for training
	std::map<unsigned int, unsigned int> selected_indices;
	// add new data indices from receptive fields of winner nodes */
	inline void addNewDataIndices (std::vector<unsigned int>&);

};

template<typename T, typename S>
ActiveLLRGNGAlgorithm<T,S>::~ActiveLLRGNGAlgorithm () {
	if (this->mdl_history != NULL)
		this->closeMDLHistory ();
}


/** \brief Runs an algorithm based on different implementations of GNG and own ideas
*/
template<typename T, typename S>
void ActiveLLRGNGAlgorithm<T,S>::run()
{
	assert (this->getDimension() > 0);
	//setRefVectors() must be called before
	//setData and/or addData must be called before
	assert (this->_graphptr);
	assert (this->insertion_rate);

	this->calculateValueRange ();
	if (this->sampling_mode == randomly && this->stopping_criterion == stability)
	{
		::srand( (unsigned)time( NULL ) );
		this->epoch = 0;
		this->stable_graph = false;
		this->min_mdl_graphptr = NULL;
		this->last_epoch_mdl_reduction = 0;
		for (unsigned int i=0; i<this->_graphptr->size(); i++)
		{
			this->_graphptr->setLastEpochImprovement (i, 0);
			LLRGNGNode<T,S>* node = static_cast<LLRGNGNode<T,S>* > (&(*this->_graphptr)[i]);
			node->last_avgerror = node->prev_avgerror = 0;
			node->errors.clear();
			node->dim_errors.clear();
		}
		do {
			for(unsigned int t = 0; t < this->insertion_rate; t++)
			{
				unsigned int selected_index = ::rand() % selected_indices.size();
				std::vector<unsigned int> winners = this->learning_loop (selected_indices[selected_index], t);
				/*if (this->_graphptr->getNeighborsSize(selected_index) > winners.size())
					this->updateMinimalGraphMDL();*/
				addNewDataIndices (winners);
				if (this->stable_graph)
					break;
			}
			if (this->mdl_history != NULL)
				this->saveMDLHistory ();
			if (this->stable_graph)
				break;
			this->epoch++;
				
		}
		while (true);
		std::cout << "Finished!" << std::endl;

	}
	else
		std::cerr << "This algorithm only works with *stability* stopping critterion and random sampling" << std::endl;
}

/** \brief Adds new data to the current data set to be processed in the next algorithmic run */
template<typename T, typename S>
inline void ActiveLLRGNGAlgorithm<T,S>::addData (std::vector< Vector<T>* >* data)
{
	if (this->_data == NULL)
		this->_data = new std::vector< Vector<T>* >;
	selected_indices.clear ();
  
	for (unsigned int i=0; i < data->size(); i++) {
		Vector<T> *new_item = new Vector<T>(*data->at(i));
		this->_data->push_back (new_item);
		selected_indices[selected_indices.size()] = this->size() - 1;
	}

}

/** \brief add new data indices from receptive fields of winner nodes */
template<typename T, typename S>
inline void ActiveLLRGNGAlgorithm<T,S>::addNewDataIndices (std::vector<unsigned int>& winners)
{
	for (unsigned int i=0; i<winners.size(); i++)
	{
		if (winners[i] < this->_graphptr->size())
		{
			LLRGNGNode<T,S>* winner_node = static_cast<LLRGNGNode<T,S>* > (&(*this->_graphptr)[winners[i]]);

			for (unsigned int i=0; i<winner_node->data.size(); i++)
				if (selected_indices.find (winner_node->data[i]) == selected_indices.end())
					selected_indices[selected_indices.size()] = winner_node->data[i];
		}
	}
	
}

/** \brief Initialize data (only for visualization in Qt Widget) */
template<typename T, typename S> void ActiveLLRGNGAlgorithm<T,S>::initializeDataVisualization ()
{
	if (this->visualizing)
		emit initializeData (this->_data, this->_graphptr->getNodes());
}

} // namespace neuralgas

#endif
