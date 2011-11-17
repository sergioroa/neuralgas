#include <QtCore/QCoreApplication>
#include "Voronoi.h"

using namespace neuralgas;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    char* filename_data=new char[255];
    char* filename_neurons=new char[255];
    char* filename_voronoi=new char[255];
    char* output=new char[255];
    strcpy(filename_data, "data.txt");
    strcpy(filename_neurons, "nodes.txt");
    strcpy(filename_voronoi, "voronoi");
    strcpy(output, "voronoi.png");

    Voronoi v;
    if (argc > 1)
	    if (std::string(argv[1]) == "-w")
		    v.setWhiteBackground ();
    if (!v.getData(filename_data))
    {
	    std::cerr << "unable to read file data.txt" << std::endl;
	    return 1;
    }
    v.showData();
    //v.setSize(500,1000);
    if (!v.getNeurons(filename_neurons))
    {
	    std::cerr << "unable to read file nodes.txt" << std::endl;
	    return 1;
    }
    v.getMaxMinValue ();
    v.calcVoronoiGnuplot();
    v.saveVoronoiGnuplot(filename_voronoi, filename_data, filename_neurons);
    v.setSizefromData(500);
    v.discretizeData ();
    v.discretizeNeurons ();
    v.calcVoronoiImage();
    v.saveVoronoiImage(output);
    std::cout <<"run done"<<std::endl;
    return 0;
}
