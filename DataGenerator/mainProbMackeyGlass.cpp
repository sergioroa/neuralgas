#include <cstdlib>
#include <iostream>
#include <DataGenerator/ProbMackeyGlass.h>

using namespace std;

int main(int argc, char *argv[])
{
    ProbMackeyGlass ph;
    ph.generate(1000);
    return EXIT_SUCCESS;
}
