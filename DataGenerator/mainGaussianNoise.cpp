#include <iostream>
#include <GrowingNeuralGas/LifelongBasedGNGAlgorithm/LLBGNGAlgorithm.h>
#include <DataGenerator/GaussianNoise.h>
#include "GrowingNeuralGas/Testing/ErrorTesting.h"

using namespace std;
using namespace neuralgas;


int main(int argc, char *argv[])
{
	
	LLBGNGAlgorithm<float, int>* llbgng = new LLBGNGAlgorithm<float, int>(2);
	int size;
	float whitenoise_prob = -1.0;
	if (argc >= 2) {
		size = atoi(argv[1]);
	}
	if (argc >= 3) {
		whitenoise_prob = atof(argv[2]);
		if (whitenoise_prob < 0.0 || whitenoise_prob > 1.0)
		{
			cerr << "Please enter a valid probability value!" << endl;
			return 1;
		}
	}
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " size [whitenoise_prob]" << endl;
		return 1;
	}


	GaussianNoise gn;
	//gn.setCanonicalDataset ();
	gn.setCustomizedDataset ();
	if (whitenoise_prob != -1.0)
		gn.setWhiteNoiseProb (whitenoise_prob);
	gn.generate(size);
	vector<Vector<float>*>* data = gn.getData();
	
	for (unsigned int i=0; i < data->size(); i++)
		std::cout <<data->operator[](i)->operator[](0)<<" "<<data->operator[](i)->operator[](1)<<std::endl;

	gn.save("data.txt");
	llbgng->setData(gn.getData());
	Vector<float> mins = llbgng->minValues();
	Vector<float> maxs = llbgng->maxValues();
	// float min = llbgng->minValue();
	// float max = llbgng->maxValue();
	for (unsigned int i=0; i< mins.size(); i++)
	{
		cout << "min: " << mins[i] << endl;
		cout << "max: " << maxs[i] << " for dim " << i << endl;
	}
	// cout << "min: " << min << endl;
	// cout << "max: " << max << endl;
	
	llbgng->setRefVectors(2,mins,maxs);
	// llbgng->setRefVectors(2,mins,maxs);

	llbgng->setTimeWindows (20, 100, 100);
	llbgng->setLearningRates (0.1, 0.001, 0.01);
	llbgng->setInsertionRate (10);
	llbgng->setAdaptationThreshold (0.05);
	llbgng->setInsertionTolerance (0.1);
	llbgng->setDeletionThreshold (0.5);
	llbgng->setMinimalNodeAge (0.001);
	llbgng->setMaximalEdgeAge (50);
	//llbgng->setMaxNodes (5);
	llbgng->setStabilization (1.001);

	llbgng->setSamplingMode (randomly);
	//llbgng->setStoppingCriterion (epochs);
	//llbgng->setMaxEpochs (1);
	llbgng->setStoppingCriterion (stability);

	llbgng->run(); 
	llbgng->save("nodes.txt");
	llbgng->showGraph ();

	ErrorTesting<float,int> et(llbgng);
	float total_error=0.0;
	std::vector<float> errors = et.getErrors(size-1);
	for (unsigned int i =0; i < errors.size(); i++)
		total_error += errors[i];
	
	std::cout << total_error / size <<std::endl;

	
	return EXIT_SUCCESS;
}
