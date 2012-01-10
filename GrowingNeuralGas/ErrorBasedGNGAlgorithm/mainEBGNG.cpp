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

#include <cstdlib>
#include <iostream>

#include "EBGNGAlgorithm.h"
#include <DataGenerator/MackeyGlass.h>
#include <DataGenerator/BinaryAutomata.h>
#include <DataGenerator/NoisyAutomata.h>
#include <GrowingNeuralGas/Testing/ErrorTesting.h>


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
{return 3;}

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
      return 3;
}

int main(int argc, char *argv[])
{
    EBGNGAlgorithm<double,int>* eb= new EBGNGAlgorithm<double,int>(/*2*/1);
    
    eb->setFuncArray(constgamma,3);
    eb->setFuncArray(func,4);
    eb->setFuncArray(functheta,7);

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
    
   
    ErrorTesting<double,int> et(eb);
    MackeyGlass *mg = new MackeyGlass;
    mg->setPastTimeSteps(17);
    mg->setBoundary(0.4);
    mg->setPower(10);
    
    mg->generate(size);
    
        
    eb->setData(mg->getData());

    vector<Vector<double>*>* data = mg->getData();
    
    for (unsigned int i=0; i < data->size(); i++)
	    std::cout <<data->operator[](i)->operator[](0)<<" "<<data->operator[](i)->operator[](1)<<std::endl;

    Vector<double> mins = minValues(data);
    Vector<double> maxs = maxValues(data);
    // double min = eb->minValue();
    // double max = eb->maxValue();
    // eb->setRefVectors(2,min,max);
    eb->setRefVectors(2,mins,maxs);
    eb->setErrorThreshold(0.01);
   
    eb->setSamplingMode(randomly);
    // eb->setStoppingCriterion (epochs);
    // eb->setMaxEpochs (5);
    eb->setStoppingCriterion (stability);

    eb->run(); 
    // eb->save("nodes.txt");
    double total_error=0.0;
    std::vector<double> errors = et.getErrors(50);
    for (unsigned int i =0; i < errors.size(); i++)
        total_error += errors[i];
    
    std::cout << total_error <<std::endl;
    return EXIT_SUCCESS;
}
