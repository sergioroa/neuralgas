#include <cstdlib>
#include <iostream>
#include "MackeyGlass.h"

using namespace std;

static float* memory;

const float funcMackeyGlass(const int& time)
{

  
 if (time>17)
    return -0.1*memory[time]+(0.2*memory[time-17])/(1+pow(memory[time-17],10));
 else
    return -0.1*memory[time]+(0.2*memory[0])/(1+pow(memory[0],10));
 
}


int main(int argc, char *argv[])
{
    MackeyGlass mg;
    mg.setPastTimeSteps(17);
    mg.setBoundary(0.4);
    mg.setPower(10);
    int size = 100;
    for (int i =0; i < size; i++)
        std::cout << (*mg.next())[0] << " ";
    std::cout << std::endl;
    memory = new float[size];
    memory[0]=0.4;
    for(int i=0; i < size-1; i++)
    {
     float dx = funcMackeyGlass(i);
     memory[i+1]= memory[i] + dx;
     std::cout << memory[i]<<" ";
    }
    
    system("PAUSE");
    return EXIT_SUCCESS;
}
