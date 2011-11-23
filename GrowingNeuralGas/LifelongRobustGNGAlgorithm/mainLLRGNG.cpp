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
		("transprob,p", po::value (&transProb)->default_value (0.1), "transition probability parameter")
		("symbinput,i", "store CrySSMEx dataset input in symbolic format")
		("symboutput,o", "store CrySSMEx dataset output in symbolic format");
    
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

	QApplication *a;
	VoronoiMainWindow *vWindow;
	LLRGNGAlgorithm<double, int>* llrgng = new LLRGNGAlgorithm<double, int>(2);
	if (vm.count("debug"))
	{
		a = new QApplication (argc, argv);
		vWindow = new VoronoiMainWindow;
		// llrgng->allocGUI (/*argc, argv*/);
		vWindow->setMutex (llrgng->getMutex ());
		vWindow->setWaitCondition (llrgng->getWaitCondition ());
		llrgng->setDebugging (true);
		// vWindow->moveToThread(QApplication::instance()->thread());
		llrgng->connect (llrgng, SIGNAL(updateData(SeqNeurons*)), vWindow, SLOT (updateData(SeqNeurons*)));
		llrgng->connect (llrgng, SIGNAL(initializeData(SeqData*, SeqNeurons*, unsigned int)), vWindow, SLOT (initializeData(SeqData*, SeqNeurons*, unsigned int)));
	}
	
	cout << sigma << "," << transProb << endl;
	
	// noisy automata testing
	NoisyAutomata na;
	na.setSigma(sigma);
	na.setTransProb(transProb);
	if (vm.count("symbinput"))
		na.setInputFormat (symbolic);
	else
		na.setInputFormat (vectorial);
	if (vm.count("symboutput"))
		na.setOutputFormat (symbolic);
	else
		na.setOutputFormat (vectorial);
	na.openCrySSMExFile ();
	na.generate(size);
	na.save("data.dat");
	llrgng->setData(na.getData());
	Vector<double> mins = minValues(na.getData());
	Vector<double> maxs = maxValues(na.getData());
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
	if (vm.count("debug"))
		vWindow->show();

	// llrgng->setTimeWindows (20, 100, 100);
	llrgng->setTimeWindows (100, 60, 80*size);
	llrgng->setLearningRates (0.1, 0.001);
	llrgng->setInsertionRate (size);
	llrgng->setAdaptationThreshold (0.0);
	llrgng->setMaximalEdgeAge (50);
	llrgng->setDataAccuracy (0.0001);
	// llrgng->setDataAccuracy (0.00000000001);
	// llrgng->setMaxNodes (3);
	llrgng->setMaxEpochsErrorReduction (5);
	llrgng->setMaxEpochsMDLReduction (400);
	llrgng->setModelEfficiencyConst (1);

	llrgng->setSamplingMode (randomly);
	//llrgng->setStoppingCriterion (epochs);
	//llrgng->setMaxEpochs (100);
	llrgng->setStoppingCriterion (stability);

	
	llrgng->begin();

	if (vm.count("debug"))
		a->exec();
	
	llrgng->wait ();
	llrgng->save("nodes.dat");
	llrgng->showGraph ();

	ErrorTesting<double,int> et(llrgng);
	double total_error=0.0;
	std::vector<double> errors = et.getErrors(size-1);
	for (unsigned int i =0; i < errors.size(); i++)
		total_error += errors[i];
	
	std::cout << total_error / size<<std::endl;
	na.reset();
	delete llrgng;

	return EXIT_SUCCESS;

}
