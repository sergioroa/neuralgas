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

#include <iostream>
#include <GrowingNeuralGas/LifelongRobustGNGAlgorithm/ActiveLLRGNGAlgorithm.h>
#include <DataGenerator/GaussianNoise.h>
// #include "GrowingNeuralGas/Testing/ErrorTesting.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace neuralgas;


int main(int argc, char *argv[])
{
	int size;
	// int insertion_rate;
	double whitenoise_prob;
	string dataset;
	unsigned int max_epochs_mdl_reduction;
	unsigned int max_epochs_error_reduction;
	double model_efficiency;
	double data_accuracy;

	po::options_description desc("Allowed parameters:");
	desc.add_options()
		("help,h", "produce help message")
		("visualize,v", "use a Qt Widget for Voronoi visualization")
		("size,s", po::value (&size)->default_value (1000), "dataset size")
		// ("insertionrate,i", po::value (&insertion_rate)->default_value (1000), "insertion rate")
		("maxepochsmdl,m", po::value (&max_epochs_mdl_reduction)->default_value (300), "Set maximum nr of epochs after mdl reduction is expected")
		("maxepochserror,e", po::value (&max_epochs_error_reduction)->default_value (3), "Set maximum nr of epochs after error reduction is expected")
		("datafile,f", po::value (&dataset)->default_value ("1"), "dataset file: use 1 for canonical and 2 for customized input from console")
		("whitenoise_prob,w", po::value (&whitenoise_prob), "white noise probability parameter")
		("mdl,o", "save MDL history to a file mdl.txt")
		("modelefficiency,c", po::value (&model_efficiency)->default_value (1), "model efficiency constant")
		("dataaccuracy,y", po::value (&data_accuracy)->default_value (0.001), "data accuracy constant");

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
	
	QApplication *a = NULL;
	VoronoiMainWindow *vWindow = NULL;
	ActiveLLRGNGAlgorithm<double, int>* llrgng = new ActiveLLRGNGAlgorithm<double, int>(2);
	if (vm.count("visualize"))
	{
		a = new QApplication (argc, argv);
		vWindow = new VoronoiMainWindow;
		vWindow->setMutex (llrgng->getMutex ());
		vWindow->setWaitCondition (llrgng->getWaitCondition ());
		llrgng->setVisualizing (true);
		// vWindow->moveToThread(QApplication::instance()->thread());
		llrgng->connect (llrgng, SIGNAL(updateData(SeqNeurons*)), vWindow, SLOT (updateData(SeqNeurons*)));
		llrgng->connect (llrgng, SIGNAL(initializeData(SeqData*, SeqNeurons*, unsigned int)), vWindow, SLOT (initializeData(SeqData*, SeqNeurons*, unsigned int)));
	}

	if (vm.count("whitenoise_prob")) {
		if (whitenoise_prob < 0.0 || whitenoise_prob > 1.0)
		{
			cerr << "Please enter a valid probability value!" << endl;
			return 1;
		}
	}
		
	GaussianNoise* gn = new GaussianNoise;
	if (dataset == "1")
		gn->setCanonicalDataset ();
	else if (dataset == "2")
		gn->setCustomizedDataset ();
	else
	{
		const char* filename = dataset.c_str();
		if (!gn->readCustomizedDataset (filename))
		{
			cerr << "Error opening file..." << endl;
			return 1;
		}
	}
	if (vm.count("whitenoise_prob"))
		gn->setWhiteNoiseProb (whitenoise_prob);
	gn->generate(size);
	vector<Vector<double>*>* data = gn->getData();
	
	for (unsigned int i=0; i < data->size(); i++)
		std::cout <<data->operator[](i)->operator[](0)<<" "<<data->operator[](i)->operator[](1)<<std::endl;

	gn->save("data.dat");

	// unsigned int subset_size = gn->getData()->size() / 10;
	unsigned int subset_size = 50;

	// for (unsigned int i=0; i<10; i++)
	for (unsigned int i=0; i<gn->getData()->size() / subset_size; i++)
	{
		vector<Vector<double>* >* subset = new vector<Vector<double>* >;
		subset->insert (subset->end(), data->begin()+i*subset_size, data->begin()+(i+1)*subset_size);
		// copy(data->begin()+i*subset_size,data->begin()+(i+1)*subset_size,subset->begin());

		llrgng->addData (subset);
		Vector<double> mins = llrgng->minValues();
		Vector<double> maxs = llrgng->maxValues();
		for (unsigned int j=0; j< mins.size(); j++)
		{
			cout << "min: " << mins[i] << endl;
			cout << "max: " << maxs[i] << " for dim " << i << endl;
		}
		
		if (i == 0)
		{
			llrgng->setRefVectors(2,mins,maxs);
			// llrgng->setTimeWindows (20, 100, 100);
			// llrgng->setTimeWindows (100, 60, size);
			llrgng->setTimeWindows (20, 12, size);
			llrgng->setLearningRates (0.8, 0.001);
			llrgng->setInsertionRate (subset_size);
			llrgng->setAdaptationThreshold (0.0);
			llrgng->setMaximalEdgeAge (50);
			llrgng->setDataAccuracy (data_accuracy);
			llrgng->setMaxEpochsErrorReduction (max_epochs_error_reduction);
			llrgng->setMaxEpochsMDLReduction (max_epochs_mdl_reduction);
			if (dataset == "1")
				llrgng->setModelEfficiencyConst (1);
			else
				llrgng->setModelEfficiencyConst (model_efficiency);
			
			llrgng->setSamplingMode (randomly);
			llrgng->setStoppingCriterion (stability);
		        // llrgng->setMeanDistanceMode (arithmetic);
			if (vm.count("mdl"))
				llrgng->saveMDLHistory ("mdl.txt");
			
		}
		else llrgng->initializeDataVisualization ();
		
		if (vm.count("visualize"))
		{
			vWindow->show ();
		}

		llrgng->begin();
		
		if (vm.count("visualize"))
			a->exec();

		llrgng->wait ();
		
		delete subset;
	}

	llrgng->save("nodes.dat");
	llrgng->showGraph ();
	// ErrorTesting<double,int> et(llrgng);
	// double total_error=0.0;
	// std::vector<double> errors = et.getErrors(size-1);
	// for (unsigned int i =0; i < errors.size(); i++)
	// 	total_error += errors[i];
	
	// std::cout << total_error / size <<std::endl;
	delete gn;
	delete llrgng;
	if (vm.count("visualize"))
	{
		delete vWindow;
		delete a;
	}
	
	return EXIT_SUCCESS;
}
