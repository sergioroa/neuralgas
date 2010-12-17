#include <cstdlib>
#include <iostream>
#include <DataGenerator/BinaryAutomata.h>

using namespace std;

int main(int argc, char *argv[])
{
    BinaryAutomata ba;
    ba.generate(1000);
    std::vector<Vector<float>* > * data = ba.getData();
    for(int i=0; i < data->size(); i++)
            std::cout << (*((*data)[i]))[0] << " ";
    return EXIT_SUCCESS;
}
