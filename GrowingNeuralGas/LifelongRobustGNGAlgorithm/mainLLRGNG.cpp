#include <iostream>
#include "LLRGNGAlgorithm.h"
#include <DataGenerator/NoisyAutomata.h>
#include "GrowingNeuralGas/Testing/ErrorTesting.h"


using namespace std;
using namespace neuralgas;

int main (int argc, char* argv[])
{
	LLRGNGAlgorithm<double, int>* llrgng = new LLRGNGAlgorithm<double, int>(2);

	int size;
	double sigma, transProb;
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
	na.setInputFormat (vectorial);
	na.openCrySSMExFile ();
	na.generate(size);
	na.save("data.txt");
	llrgng->setData(na.getData());
	Vector<double> mins = llrgng->minValues();
	Vector<double> maxs = llrgng->maxValues();
	// double min = llrgng->minValue();
	// double max = llrgng->maxValue();
	for (unsigned int i=0; i< mins.size(); i++)
	{
		cout << "min: " << mins[i] << endl;
		cout << "max: " << maxs[i] << " for dim " << i << endl;
	}
	// cout << "min: " << min << endl;
	// cout << "max: " << max << endl;
	
	llrgng->setRefVectors(2,mins,maxs);
	// llrgng->setRefVectors(2,mins,maxs);

	// llrgng->setTimeWindows (20, 100, 100);
	llrgng->setTimeWindows (100, 60, 80*size);
	llrgng->setLearningRates (0.1, 0.001);
	llrgng->setInsertionRate (size);
	llrgng->setAdaptationThreshold (0.0);
	llrgng->setMaximalEdgeAge (100);
	llrgng->setDataAccuracy (0.001);
	// llrgng->setDataAccuracy (0.00000000001);
	// llrgng->setMaxNodes (3);
	llrgng->setMaxEpochsErrorReduction (5);
	llrgng->setMaxEpochsMDLReduction (60);

	llrgng->setSamplingMode (randomly);
	//llrgng->setStoppingCriterion (epochs);
	//llrgng->setMaxEpochs (100);
	llrgng->setStoppingCriterion (stability);

	llrgng->run(); 
	llrgng->save("nodes.txt");
	llrgng->showGraph ();

	ErrorTesting<double,int> et(llrgng);
	double total_error=0.0;
	std::vector<double> errors = et.getErrors(size-1);
	for (unsigned int i =0; i < errors.size(); i++)
		total_error += errors[i];
	
	std::cout << total_error / size<<std::endl;

	return EXIT_SUCCESS;

}
