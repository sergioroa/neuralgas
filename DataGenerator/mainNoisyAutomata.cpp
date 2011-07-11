#include <cstdlib>
#include <iostream>
#include <DataGenerator/NoisyAutomata.h>

using namespace std;
using namespace neuralgas;


int main(int argc, char *argv[])
{
    NoisyAutomata na;
    na.generate(1000);
    vector<Vector<double>*>* data = na.getData();
  
    for (unsigned int i=0; i < data->size(); i++)
        std::cout <<data->operator[](i)->operator[](0)<<" "<<data->operator[](i)->operator[](1)<<std::endl;
    return EXIT_SUCCESS;
}
