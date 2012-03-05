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

#include <cstdlib>
#include <iostream>
#include "ErrorTesting.h"
#include <GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGAlgorithm.h>
#include "tools/helpers.h"

using namespace std;
using namespace neuralgas;

int main(int argc, char *argv[])
{

    std::vector < Base_Node<double, int>* > nodes;
    if (!readNodes ("nodes.dat", &nodes))
    {
	std::cerr << "unable to read file nodes.dat" << std::endl;
	return 1;
    }
    
    assert (nodes.size());
    unsigned int dim = nodes[0]->weight.size ();

    std::vector<Vector<double>* > data;
    if (!readData ("data.dat", &data))
    {
	std::cerr << "unable to read file data.dat" << std::endl;
	return 1;
    }
    assert (data.size());
    assert (data.at(0)->size() == dim);
    
    LLRGNGAlgorithm<double, int> gng (dim);
    gng.setNodes ( &nodes );
    gng.setData ( &data );
    
    gng.showGraph ();

    ErrorTesting<double,int> et(&gng);

    std::vector<double> errors = et.getErrors(data.size());
    double total_error=0.0;
    for(unsigned int j=0; j < errors.size(); j++)
	total_error+=errors[j];
    std::cout << "Avg error: "<< total_error / data.size() <<std::endl;
    std::cout << std::endl;
    
    return EXIT_SUCCESS;
}
