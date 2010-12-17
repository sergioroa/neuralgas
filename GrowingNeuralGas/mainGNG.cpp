#include <cstdlib>
#include <iostream>
#include "GNGAlgorithm.h"

using namespace std;

float func(const int& time)
{return 0.5 ;}

float funclambda(const int& time)
{return 5 ;}

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
    GNGAlgorithm<float,float> mgng(1);
    
    
    for(int i=0; i < NUM_PARAM; i++)
            mgng.setFuncArray(func,i);

    mgng.setFuncArray(functheta,7);
    mgng.setFuncArray(funclambda,6);

    mgng.setData(v);
    int max_value = mgng.maxRandomValue(); 
    mgng.setRefVectors(2,max_value);
    mgng.run();
    mgng.showGraph();
    for( int i=0; i < 1500; i++)
         delete (*v)[i];
    
    //delete v;

    return EXIT_SUCCESS;
}
