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


float func(const unsigned int& time)
{return 0.5 ;}


float constalpha(const unsigned int& time)
{return 0.5;}

float constbeta(const unsigned int& time)
{return 0.75;}

float constgamma(const unsigned int& time)
{return 88;}

float constepsilonw(const unsigned int& time)
{return 0.05;}

float constepsilonn(const unsigned int& time)
{return 0.0006;}

float constdelta(const unsigned int& time)
{return 0.5;}

float consteta(const unsigned int& time)
{return 0.99995;}

float consttheta(const unsigned int& time)
{return 75;}

float constlambda(const unsigned int& time)
{return 600 ;}

float funcalpha(const unsigned int& time)
{
      float random_alpha = float(rand()) / float(50000);
      return random_alpha;
}


float funcbeta(const unsigned int& time)
{
     return 0.75 ;
}

float funcgamma(const unsigned int& time)
{return 1+sqrt(time);}

float funclambda(const unsigned int& time)
{return 3 ;}

float functheta(const unsigned int& time)
{
      return 75;
}

int main(int argc, char *argv[])
{
    EBGNGAlgorithm<float,int>* eb= new EBGNGAlgorithm<float,int>(/*2*/1);
    
    eb->setFuncArray(constgamma,3);
    eb->setFuncArray(func,4);
    eb->setFuncArray(functheta,7);

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
    /*NoisyAutomata na;
    na.setSigma(sigma);
    na.setTransProb(transProb);
    na.generate(size);
    na.save("data.txt");
    eb->setData(na.getData());*/
    //binary automata testing
  /*  BinaryAutomata ba;
    ba.generate(size);
    eb->setData(ba.getData());
   */ 
    //Mackey Glass testing
    
   
    ErrorTesting<float,int> et(eb);
    MackeyGlass *mg = new MackeyGlass;
    mg->setPastTimeSteps(17);
    mg->setBoundary(0.4);
    mg->setPower(10);
    
    mg->generate(size);
    
        
    eb->setData(mg->getData());

    vector<Vector<float>*>* data = mg->getData();
    
    for (unsigned int i=0; i < data->size(); i++)
	    std::cout <<data->operator[](i)->operator[](0)<<" "<<data->operator[](i)->operator[](1)<<std::endl;

    Vector<float> mins = eb->minValues();
    Vector<float> maxs = eb->maxValues();
    // float min = eb->minValue();
    // float max = eb->maxValue();
    // eb->setRefVectors(2,min,max);
    eb->setRefVectors(2,mins,maxs);
    eb->setErrorThreshold(0.01);
   
    eb->setSamplingMode(randomly);
    // eb->setStoppingCriterion (epochs);
    // eb->setMaxEpochs (5);
    eb->setStoppingCriterion (stability);

    eb->run(); 
    // eb->save("nodes.txt");
    float total_error=0.0;
    std::vector<float> errors = et.getErrors(50);
    for (unsigned int i =0; i < errors.size(); i++)
        total_error += errors[i];
    
    std::cout << total_error <<std::endl;
    return EXIT_SUCCESS;
}
