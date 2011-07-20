#include <iostream>
#include <GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGAlgorithm.h>
#include <DataGenerator/GaussianNoise.h>
#include "GrowingNeuralGas/Testing/ErrorTesting.h"

using namespace std;
using namespace neuralgas;


int main(int argc, char *argv[])
{
	
	LLRGNGAlgorithm<double, int>* llrgng = new LLRGNGAlgorithm<double, int>(2);
	int size;
	double whitenoise_prob = -1.0;
	string dataset;
	if (argc >= 3) {
		size = atoi(argv[2]);
		dataset = string (argv[1]);
	}
	if (argc >= 4) {
		whitenoise_prob = atof(argv[3]);
		if (whitenoise_prob < 0.0 || whitenoise_prob > 1.0)
		{
			cerr << "Please enter a valid probability value!" << endl;
			return 1;
		}
	}
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " 1/2/file size [whitenoise_prob]" << endl;
		cerr << "where 1/2/file means canonical, customized or read dataset (file)" << endl;
		return 1;
	}


	GaussianNoise gn;
	if (dataset == "1")
		gn.setCanonicalDataset ();
	else if (dataset == "2")
		gn.setCustomizedDataset ();
	else
	{
		const char* filename = dataset.c_str();
		if (!gn.readCustomizedDataset (filename))
		{
			cerr << "Error opening file..." << endl;
			return 1;
		}
	}
	if (whitenoise_prob != -1.0)
		gn.setWhiteNoiseProb (whitenoise_prob);
	gn.generate(size);
	vector<Vector<double>*>* data = gn.getData();
	
	for (unsigned int i=0; i < data->size(); i++)
		std::cout <<data->operator[](i)->operator[](0)<<" "<<data->operator[](i)->operator[](1)<<std::endl;

	gn.save("data.txt");
	llrgng->setData(gn.getData());
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
	llrgng->setTimeWindows (100, 60, size);
	llrgng->setLearningRates (0.1, 0.001);
	llrgng->setInsertionRate (size);
	llrgng->setAdaptationThreshold (0.0);
	llrgng->setMaximalEdgeAge (50);
	llrgng->setDataAccuracy (0.001);
	// llrgng->setDataAccuracy (0.00000000001);
	//llrgng->setMaxNodes (5);
	llrgng->setMaxEpochsErrorReduction (5);
	llrgng->setMaxEpochsMDLReduction (80);
	
	llrgng->setSamplingMode (randomly);
	//llrgng->setStoppingCriterion (epochs);
	//llrgng->setMaxEpochs (1);
	llrgng->setStoppingCriterion (stability);

	llrgng->run(); 
	llrgng->save("nodes.txt");
	llrgng->showGraph ();

	ErrorTesting<double,int> et(llrgng);
	double total_error=0.0;
	std::vector<double> errors = et.getErrors(size-1);
	for (unsigned int i =0; i < errors.size(); i++)
		total_error += errors[i];
	
	std::cout << total_error / size <<std::endl;

	
	return EXIT_SUCCESS;
}
