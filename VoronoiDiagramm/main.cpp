#include <QtCore/QCoreApplication>
#include "Voronoi.h"

using namespace neuralgas;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    char* filename_data=new char[255];
    char* filename_neurons=new char[255];
    char* output=new char[255];
    strcpy(filename_data, "data.txt");
    strcpy(filename_neurons, "nodes.txt");
    strcpy(output, "voronoi.jpg");

    Voronoi v;
    v.getData(filename_data);
    v.showData();
    //v.setSize(500,1000);
    v.getNeurons(filename_neurons);
    v.setSizefromData(500);
    v.calcVoronoi();
    v.save(output);
    std::cout <<"run done"<<std::endl;
    return 0;
}
