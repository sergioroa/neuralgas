#include <NeuralGas/Vector.h>
#include <iostream>

int main (int argc, char* argv[])
{
	Vector<double> myVector;

	myVector = Vector<double>(3, 0.0);

	myVector[0] += 2.0;
	myVector[1] += 1.0;

	std::cout << std::endl << myVector[0] << "\t";
	std::cout << myVector[1] << "\t";
	std::cout << myVector[2] << std::endl;


}
