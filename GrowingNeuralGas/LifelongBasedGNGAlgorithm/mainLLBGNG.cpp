#include <iostream>
#include "LLBGNGAlgorithm.h"
#include <DataGenerator/NoisyAutomata.h>


using namespace std;
using namespace neuralgas;

int main (int argc, char* argv[])
{
	LLBGNGAlgorithm<float, int>* llbgng = new LLBGNGAlgorithm<float, int>(2);

	int size;
	float sigma, transProb;
	if (argc == 4) {
		size = atoi(argv[1]);
		sigma = atof(argv[2]);
		transProb = atof(argv[3]);
	}
	else {
		cerr << "Usage: " << argv[0] << " size sigma transition_prob" << endl;
		return 1;
	}
    
	// noisy automata testing
	NoisyAutomata na;
	na.setSigma(sigma);
	na.setTransProb(transProb);
	na.generate(size);
	na.save("data.txt");
	llbgng->setData(na.getData());
	float min = llbgng->minValue();
	float max = llbgng->maxValue();
	llbgng->setRefVectors(2,min,max);

	llbgng->setTimeWindows (20, 100, 100);
	llbgng->setLearningRates (0.1, 0.01, 0.1);
	llbgng->setInsertionRate (10);
	llbgng->setAdaptationThreshold (0.05);
	llbgng->setInsertionTolerance (0.2);
	llbgng->setDeletionThreshold (0.1);
	llbgng->setMinimalNodeAge (0.01);
	llbgng->setMaximalEdgeAge (50);
	llbgng->setStabilization (0.01);

	llbgng->setSamplingMode (randomly);
	llbgng->setStoppingCriterion (epochs);
	llbgng->setMaxEpochs (50);

	llbgng->run(); 
	llbgng->save("nodes.txt");
	return EXIT_SUCCESS;

}
