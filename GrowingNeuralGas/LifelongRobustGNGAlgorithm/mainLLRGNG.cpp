#include <iostream>
#include "LLRGNGAlgorithm.h"
#include <DataGenerator/NoisyAutomata.h>
#include "GrowingNeuralGas/Testing/ErrorTesting.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;


using namespace std;
using namespace neuralgas;

int main (int argc, char* argv[])
{
	int size;
	double sigma, transProb;

	po::options_description desc("Allowed parameters:");
	desc.add_options()
		("help,h", "produce help message")
		("debug,d", "debug using a Qt Widget for Voronoi visualization")
		("size", po::value (&size)->default_value (1000), "dataset size")
		("sigma", po::value (&sigma)->default_value (0.1), "sigma std dev parameter")
		("transprob,p", po::value (&transProb)->default_value (0.1), "transition probability parameter");
    
	// Declare an options description instance which will include
	// all the options
	po::options_description all("Basic usage:");
	all.add(desc);
  
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
		  options(all).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}

	LLRGNGAlgorithm<double, int>* llrgng = new LLRGNGAlgorithm<double, int>(2);
	if (vm.count("debug"))
		llrgng->allocGUI (argc, argv);

	cout << sigma << "," << transProb << endl;
  
	// noisy automata testing
	NoisyAutomata na;
	na.setSigma(sigma);
	na.setTransProb(transProb);
	na.setInputFormat (vectorial);
	na.setOutputFormat (vectorial);
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
	// llrgng->setRefVectors(2,min,max);

	// llrgng->setTimeWindows (20, 100, 100);
	llrgng->setTimeWindows (100, 60, 80*size);
	llrgng->setLearningRates (0.1, 0.001);
	llrgng->setInsertionRate (size);
	llrgng->setAdaptationThreshold (0.0);
	llrgng->setMaximalEdgeAge (50);
	llrgng->setDataAccuracy (0.001);
	// llrgng->setDataAccuracy (0.00000000001);
	// llrgng->setMaxNodes (3);
	llrgng->setMaxEpochsErrorReduction (5);
	llrgng->setMaxEpochsMDLReduction (60);
	llrgng->setModelEfficiencyConst (1);

	llrgng->setSamplingMode (randomly);
	//llrgng->setStoppingCriterion (epochs);
	//llrgng->setMaxEpochs (100);
	llrgng->setStoppingCriterion (stability);

	llrgng->begin();

	if (vm.count("debug"))
	{
		llrgng->getVWindow()->show();
		llrgng->getApp()->exec();
	}
	llrgng->wait ();
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
