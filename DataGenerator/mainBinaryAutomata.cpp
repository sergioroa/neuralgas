#include <cstdlib>
#include <iostream>
#include <DataGenerator/BinaryAutomata.h>

using namespace std;
using namespace neuralgas;

int main(int argc, char *argv[])
{
    BinaryAutomata ba;
    ba.generate(1000);
    std::vector<Vector<double>* > * data = ba.getData();
    for(unsigned int i=0; i < data->size(); i++)
            std::cout << (*((*data)[i]))[0] << " ";
    return EXIT_SUCCESS;
}
