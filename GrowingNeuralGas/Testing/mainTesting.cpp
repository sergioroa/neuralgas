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

    std::vector<Vector<double>* >* data = new std::vector<Vector<double>* >;
    if (!readData ("data.dat", data))
    {
	std::cerr << "unable to read file nodes.dat" << std::endl;
	return 1;
    }
    assert (data->size());
    assert (data->at(0)->size() == dim);
    
    LLRGNGAlgorithm<double, int> gng (dim);
    gng.setNodes ( &nodes );
    gng.setData ( data );
    
    gng.showGraph ();

    ErrorTesting<double,int> et(&gng);

    std::vector<double> errors = et.getErrors(data->size());
    double total_error=0.0;
    for(unsigned int j=0; j < errors.size(); j++)
	total_error+=errors[j];
    std::cout << "Avg error: "<< total_error / data->size() <<std::endl;
    std::cout << std::endl;
    
    return EXIT_SUCCESS;
}
