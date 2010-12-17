#include <cstdlib>
#include <iostream>
#include "ErrorTesting.h"
#include <GrowingNeuralGas/MergeGrowingNeuralGas/MGNGAlgorithm.h>

using namespace std;


float func(const int& time)
{return 0.5 ;}

float funclambda(const int& time)
{return 3 ;}

float functheta(const int& time)
{return 100;}


int main(int argc, char *argv[])
{

    std::vector< Vector<float>* >* v=new std::vector<Vector<float>* >;
    for( int i=0; i < 1500; i++)
    {
        Vector<float>* x;
        float odd = 3.7 * i*i -  i*5.2+2.3;
        float even =5.2 * i*i +  i*5.2+2.3;
        if (i % 2 || i%3 == 0)
           x = new Vector<float>(1,even);
        else
           x = new Vector<float>(1,odd);
        v->push_back(x);
       
    }
    MGNGAlgorithm<float,float>* mgng=new MGNGAlgorithm<float,float>(1);
    
    
    for(int i=0; i < NUM_PARAM; i++)
            mgng->setFuncArray(func,i);

    mgng->setFuncArray(functheta,6);
    mgng->setFuncArray(funclambda,8);

    mgng->setData(v);
    int max_value = mgng->maxRandomValue(); 
    mgng->setRefVectors(2,max_value);
    mgng->run();
    mgng->showGraph();
    ErrorTesting<float,float> et(mgng);

    std::vector<float> errors = et.getErrors(30,true);
    for(int i=0; i < errors.size(); i++)
            std::cout << errors [i] << " ";

    for( int i=0; i < 1500; i++)
         delete (*v)[i];
    
    //delete v;
    delete mgng;
    
    return EXIT_SUCCESS;
}
