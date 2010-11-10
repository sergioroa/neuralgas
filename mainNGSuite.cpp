#include <cstdlib>
#include <iostream>
#include "NeuralGasSuite.h"
#include "GrowingNeuralGas/Testing/ErrorTesting.h"
#include "GrowingNeuralGas/MergeGrowingNeuralGas/MGNGAlgorithm.h"
#include "GrowingNeuralGas/MergeGrowingNeuralGas/CDNAlgorithm.h"
#include "GrowingNeuralGas/GNGAlgorithm.h"
#include "GrowingNeuralGas/ErrorBasedGNGAlgorithm/EBGNGAlgorithm.h"
#include "DataGenerator/MackeyGlass.h"
#include <math.h>

using namespace std;


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
{return 100;}

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
      return 100;
}

int main(int argc, char *argv[])
{
    int sizeofdata=50000;
    
    NeuralGasSuite<float,int> ng;
    
    MackeyGlass* mg = new MackeyGlass;
    mg->setPastTimeSteps(17);
    mg->setBoundary(0.4);
    mg->setPower(10);
    
    
    ng.setDataGenerator(mg);
    
    
    mg->generate(sizeofdata);
        
    MGNGAlgorithm<float,int>* mgng    = new MGNGAlgorithm<float,int>(1);
    GNGAlgorithm<float,int>* gng      = new GNGAlgorithm<float,int>(1);
    EBGNGAlgorithm<float,int>* ebgng  = new EBGNGAlgorithm<float,int>(1);
    CDNAlgorithm<float,int>* cdn      = new CDNAlgorithm<float,int>(1);     
    
    for(int i=0; i < NUM_PARAM; i++)
    {
            gng->setFuncArray(func,i);
            ebgng->setFuncArray(func,i);
            cdn->setFuncArray(func,i);
     }
    
    mgng->setFuncArray(constalpha,0);
    mgng->setFuncArray(constbeta,1);
    mgng->setFuncArray(constgamma,2);
    mgng->setFuncArray(constdelta,3);
    mgng->setFuncArray(constepsilonw,4);
    mgng->setFuncArray(constepsilonn,5);
    mgng->setFuncArray(consttheta,6);
    mgng->setFuncArray(consteta,7);
    mgng->setFuncArray(constlambda,8);
    
        
    cdn->setFuncArray(constalpha,0);
    cdn->setFuncArray(constbeta,1);
    cdn->setFuncArray(constgamma,2);
    cdn->setFuncArray(constdelta,3);
    cdn->setFuncArray(constepsilonw,4);
    cdn->setFuncArray(constepsilonn,5);
    cdn->setFuncArray(consttheta,6);
    cdn->setFuncArray(consteta,7);
    cdn->setFuncArray(constlambda,8);
        
    gng->setFuncArray(constgamma,3);
    gng->setFuncArray(functheta,7);
    gng->setFuncArray(funclambda,6);
    ebgng->setFuncArray(constgamma,3);
    ebgng->setFuncArray(functheta,7);
    ebgng->setFuncArray(funclambda,8);


   
    ng.add(mgng);
  //  ng.add(gng);
    ng.add(cdn);
    ng.setRefVectors(2);
    (dynamic_cast< CDNAlgorithm<float,int>* > (ng[1]))->setEnergy(0.1);
    ng.run(); 
    std::vector<float> errors;
    for (int i=0; i < ng.size(); i++)
    {
        std::cout << "Algorithm " << i << std::endl;
        errors = ng.getErrors(i,500);
        float total_error=0.0;
        //ng[i]->showGraph();   
        for(int j=0; j < errors.size(); j++)
        {
         //       std::cout << errors [j] << " ";
                total_error+=errors[j];
        }
        std::cout << "Gesamtfehler "<< total_error<<std::endl;
    
        std::cout << std::endl;
    }

    std::cout << std::endl;
    delete mg;
    
    system("PAUSE");
    return EXIT_SUCCESS;
}
