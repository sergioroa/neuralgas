#include <cstdlib>
#include <iostream>
#include "NeuralGasSuite.h"
#include "GrowingNeuralGas/Testing/ErrorTesting.h"
#include "GrowingNeuralGas/MergeGrowingNeuralGas/MGNGAlgorithm.h"
#include "GrowingNeuralGas/MergeGrowingNeuralGas/CDNAlgorithm.h"
#include "GrowingNeuralGas/GNGAlgorithm.h"
#include "GrowingNeuralGas/ErrorBasedGNGAlgorithm/EBGNGAlgorithm.h"
#include "GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGAlgorithm.h"
#include "DataGenerator/MackeyGlass.h"
#include "DataGenerator/NoisyAutomata.h"
#include <math.h>

using namespace std;
using namespace neuralgas;

double func(const unsigned int& time)
{return 0.5 ;}


double constalpha(const unsigned int& time)
{return 0.5;}

double constbeta(const unsigned int& time)
{return 0.75;}

double constgamma(const unsigned int& time)
{return 88;}

double constepsilonw(const unsigned int& time)
{return 0.05;}

double constepsilonn(const unsigned int& time)
{return 0.0006;}

double constdelta(const unsigned int& time)
{return 0.5;}

double consteta(const unsigned int& time)
{return 0.99995;}

double consttheta(const unsigned int& time)
{return 100;}

double constlambda(const unsigned int& time)
{return 600 ;}

double funcalpha(const unsigned int& time)
{
      double random_alpha = double(rand()) / double(50000);
      return random_alpha;
}


double funcbeta(const unsigned int& time)
{
     return 0.75 ;
}

double funcgamma(const unsigned int& time)
{return 1+sqrt(time);}

double funclambda(const unsigned int& time)
{return 3 ;}

double functheta(const unsigned int& time)
{
      return 100;
}

enum _algorithm { _gng, _ebgng, _mgng, _cdn, _llbgng };
enum _dataset {_mg, _na };

void display_error_msg (char *argv[])
{
      cerr << "Usage: " << argv[0] << " gng/ebgng/mgng/cdn/llrgng mg/na size" << endl;
      exit(1);
}

int main(int argc, char *argv[])
{

    if (argc < 4)
	    display_error_msg (argv);

    unsigned int dataset;
    if (string(argv[2]) == "mg")
	    dataset = _mg;
    else if (string(argv[2]) == "na")
	    dataset = _na;
    else
	    display_error_msg (argv);
    
    GNGModul<double, int> *gng;
    unsigned int algorithm;
    //MGNGAlgorithm<double,int>* mgng    = new MGNGAlgorithm<double,int>(1);
    //GNGAlgorithm<double,int>* gng      = new GNGAlgorithm<double,int>(1);
    //EBGNGAlgorithm<double,int>* ebgng  = new EBGNGAlgorithm<double,int>(2);
    //CDNAlgorithm<double,int>* cdn      = new CDNAlgorithm<double,int>(1);
    if (string(argv[1]) == "gng")
    {
	    gng = new GNGAlgorithm<double,int>(dataset + 1);
	    algorithm = _gng;
    }
    else if (string(argv[1]) == "ebgng")
    {
	    gng = new EBGNGAlgorithm<double,int>(dataset + 1);
	    algorithm = _ebgng;
    }
    else if (string(argv[1]) == "mgng")
    {
	    gng = new MGNGAlgorithm<double,int>(dataset + 1);
	    algorithm = _mgng;
    }
    else if (string(argv[1]) == "cdn")
    {
	    gng = new CDNAlgorithm<double,int>(dataset + 1);
	    algorithm = _cdn;
    }
    else if (string(argv[1]) == "llrgng")
    {
	    gng = new LLRGNGAlgorithm<double,int>(dataset + 1);
	    algorithm = _llbgng;
    }
    else
	    display_error_msg (argv);
       
    int sizeofdata=atoi (argv[3]);
    
    NeuralGasSuite<double,int> ng;
    DataGenerator<double> *dG;

    if (dataset == _mg)
    {
	    dG = new MackeyGlass;
	    static_cast<MackeyGlass*>(dG)->setPastTimeSteps(17);
	    static_cast<MackeyGlass*>(dG)->setBoundary(0.4);
	    static_cast<MackeyGlass*>(dG)->setPower(10);
    }
    else if (dataset == _na)
    {
	    dG = new NoisyAutomata;
	    double sigma, transprob;
	    cout << "Sigma: ";
	    cin >> sigma;
	    cout << "Transition prob.: ";
	    cin >> transprob;
	    static_cast<NoisyAutomata*>(dG)->setSigma (sigma);
	    static_cast<NoisyAutomata*>(dG)->setTransProb (transprob);
    }
    
    
    ng.setDataGenerator(dG);
    //ng.setDataGenerator(na);
    
    
    dG->generate(sizeofdata);
    //na->generate(sizeofdata);
    //vector<Vector<double>*>* data = na->getData();
    vector<Vector<double>*>* data = dG->getData();
    
    for (unsigned int i=0; i < data->size(); i++)
    {
	    std::cout <<data->operator[](i)->operator[](0);
	    if (dataset == _na)
		    std::cout << " "<< data->operator[](i)->operator[](1);
	    std::cout << std::endl;
    }

    ng.add(gng);
    ng.setData();
    ng.setRefVectors(2);
   
    for(int i=0; i < NUM_PARAM; i++)
    {
        gng->setFuncArray(func,i);
        //ebgng->setFuncArray(func,i);
        //cdn->setFuncArray(func,i);
    }
    
    // mgng->setFuncArray(constalpha,0);
    // mgng->setFuncArray(constbeta,1);
    // mgng->setFuncArray(constgamma,2);
    // mgng->setFuncArray(constdelta,3);
    // mgng->setFuncArray(constepsilonw,4);
    // mgng->setFuncArray(constepsilonn,5);
    // mgng->setFuncArray(consttheta,6);
    // mgng->setFuncArray(consteta,7);
    // mgng->setFuncArray(constlambda,8);
    if (algorithm == _mgng)
    {
        gng->setFuncArray(constalpha,0);
        gng->setFuncArray(constbeta,1);
        gng->setFuncArray(constgamma,2);
        gng->setFuncArray(constdelta,3);
        gng->setFuncArray(constepsilonw,4);
        gng->setFuncArray(constepsilonn,5);
        gng->setFuncArray(consttheta,6);
        gng->setFuncArray(consteta,7);
        gng->setFuncArray(constlambda,8);
    }
    else if (algorithm == _cdn)
    {
        
    // cdn->setFuncArray(constalpha,0);
    // cdn->setFuncArray(constbeta,1);
    // cdn->setFuncArray(constgamma,2);
    // cdn->setFuncArray(constdelta,3);
    // cdn->setFuncArray(constepsilonw,4);
    // cdn->setFuncArray(constepsilonn,5);
    // cdn->setFuncArray(consttheta,6);
    // cdn->setFuncArray(consteta,7);
    // cdn->setFuncArray(constlambda,8);

        gng->setFuncArray(constalpha,0);
        gng->setFuncArray(constbeta,1);
        gng->setFuncArray(constgamma,2);
        gng->setFuncArray(constdelta,3);
        gng->setFuncArray(constepsilonw,4);
        gng->setFuncArray(constepsilonn,5);
        gng->setFuncArray(consttheta,6);
        gng->setFuncArray(consteta,7);
        gng->setFuncArray(constlambda,8);
    }
    else if (algorithm == _gng)
    {
        gng->setFuncArray(constgamma,3);
        gng->setFuncArray(functheta,7);
        gng->setFuncArray(funclambda,6);
        //gng->setStoppingCriterion (global_error);
        gng->setStoppingCriterion (epochs);
        gng->setMaxEpochs (100);
        gng->setSamplingMode(randomly);
        //(static_cast<GNGAlgorithm<double,int>*>(gng))->setMinGlobalError (0.01);
    }
    else if (algorithm == _ebgng)
    {
    // ebgng->setFuncArray(constgamma,3);
    // ebgng->setFuncArray(functheta,7);
    // ebgng->setFuncArray(funclambda,8);
	static_cast<EBGNGAlgorithm<double,int>*>(gng)->setErrorThreshold(0.03);
        gng->setFuncArray(constgamma,3);
        gng->setFuncArray(functheta,7);
        gng->setFuncArray(funclambda,8);
        gng->setSamplingMode(randomly);
        // gng->setStoppingCriterion (epochs);
	gng->setStoppingCriterion (stability);
        // gng->setMaxEpochs (1000000000);
    }
    else if (algorithm == _llbgng)
    {
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setTimeWindows (100, 60, sizeofdata);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setLearningRates (0.1, 0.001);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setInsertionRate (sizeofdata);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setAdaptationThreshold (0.0);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setMaximalEdgeAge (50);
	    // static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setDataAccuracy (0.001);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setDataAccuracy (0.00000000001);
	    if (dataset == _na)
		    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setModelEfficiencyConst (1.0);
	    else if (dataset == _mg)
		    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setModelEfficiencyConst (1.3);		    
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setMaxEpochsErrorReduction (5);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setMaxEpochsMDLReduction (80);
	    static_cast<LLRGNGAlgorithm<double,int>*>(gng)->setSamplingMode (randomly);
	    //gng->setStoppingCriterion (epochs);
	    //gng->setMaxEpochs (100);
	    //(static_cast<LLBGNGAlgorithm<double,int>*>(gng))->setMaxNodes(5);
	    gng->setStoppingCriterion (stability);

    }
   
    //(static_cast< CDNAlgorithm<double,int>* > (ng[1]))->setEnergy(0.1);
    //na->save("data.txt");
    dG->save("data.txt");
    ng.run();
    std::vector<double> errors;
    for (int i=0; i < ng.size(); i++)
    {
        std::cout << "Algorithm " << i << std::endl;
        //errors = ng.getErrors(i,500);
        errors = ng.getErrors(i,sizeofdata-1);
        double total_error=0.0;
        ng[i]->showGraph();
	ng[i]->save ("nodes.txt");
        for(unsigned int j=0; j < errors.size(); j++)
        {
         //       std::cout << errors [j] << " ";
                total_error+=errors[j];
        }
        std::cout << "Total error: "<< total_error / sizeofdata <<std::endl;
    
        std::cout << std::endl;
    }

    std::cout << std::endl;
    delete dG;
    //delete na;
    
    return EXIT_SUCCESS;
}
