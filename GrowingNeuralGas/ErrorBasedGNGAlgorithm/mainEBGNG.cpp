#include <cstdlib>
#include <iostream>

#include "EBGNGAlgorithm.h"
// #include <GrowingNeuralGas/ErrorBasedGNGAlgorithm/EBGNGAlgorithm.h>
#include <DataGenerator/MackeyGlass.h>
#include <DataGenerator/BinaryAutomata.h>
#include <DataGenerator/NoisyAutomata.h>
#include <GrowingNeuralGas/Testing/ErrorTesting.h>


using namespace std;
using namespace neuralgas;


float func(const int& time)
{return 0.5 ;}


float constalpha(const int& time)
{return 0.5;}

float constbeta(const int& time)
{return 0.75;}

float constgamma(const int& time)
{return 88;}

float constepsilonw(const int& time)
{return 0.05;}

float constepsilonn(const int& time)
{return 0.0006;}

float constdelta(const int& time)
{return 0.5;}

float consteta(const int& time)
{return 0.99995;}

float consttheta(const int& time)
{return 75;}

float constlambda(const int& time)
{return 600 ;}

float funcalpha(const int& time)
{
      float random_alpha = float(rand()) / float(50000);
      return random_alpha;
}


float funcbeta(const int& time)
{
     return 0.75 ;
}

float funcgamma(const int& time)
{return 1+sqrt(time);}

float funclambda(const int& time)
{return 3 ;}

float functheta(const int& time)
{
      return 75;
}

int main(int argc, char *argv[])
{
    EBGNGAlgorithm<float,int>* eb= new EBGNGAlgorithm<float,int>(2);
    
    eb->setFuncArray(constgamma,3);
    eb->setFuncArray(func,4);
    eb->setFuncArray(functheta,7);
    eb->setRefVectors(2,10);

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
    eb->setData(na.getData());   
    //binary automata testing
  /*  BinaryAutomata ba;
    ba.generate(size);
    eb->setData(ba.getData());
   */ 
    //Mackey Glass testing
   /* 
   
    ErrorTesting<float,int> et(eb);
    MackeyGlass mg;
    mg.setPastTimeSteps(17);
    mg.setBoundary(0.4);
    mg.setPower(10);
    
    mg.generate(size);
    
        
    eb->setData(mg.getData());
    eb->run();
    float total_error=0.0;
    std::vector<float> errors = et.getErrors(50);
    for (int i =0; i < errors.size(); i++)
        total_error += errors[i];
    
    std::cout << total_error <<std::endl;
   */
    //eb->setSamplingMode(randomly);
    //eb->setStoppingCriterion (epochs);
    //eb->setMaxEpochs (5);

    eb->run(); 
    eb->save("nodes.txt");
    return EXIT_SUCCESS;
}
