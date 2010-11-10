#include <cstdlib>
#include <iostream>
#include "Vector.h"
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    Vector<float>* x;
    x = new Vector<float>(10);
    std::cout << &x << std::endl;
    std::vector< Vector<float>* >* v=new std::vector< Vector<float>* >;
    std::cout << "1 Schritt"<<std::endl;
    
    for (int i=0; i < 15; i++)
    {
        v->push_back(x);
    }
    std::cout << x->getAddress()<<std::endl;
   delete x;
x=NULL;
   //delete (*v)[2];
    for (int i=0; i < 15; i++)
    {
        std::cout<<(*v)[i]->getAddress()<<std::endl;
    }
  //  std::cout << "x "<< x->getAddress()<<std::endl;    
    std::cout << "2 Schritt"<<std::endl;
    system("PAUSE");
    return EXIT_SUCCESS;
}
