#include <QtCore/QCoreApplication>
#include "Voronoi.h"

using namespace neuralgas;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    std::string filename_data="data.dat";
    std::string filename_neurons="nodes.dat";
    std::string filename_voronoi="voronoi";
    std::string output="voronoi.png";

    Voronoi v;
    if (argc > 1)
	    if (std::string(argv[1]) == "-w")
		    v.setWhiteBackground ();
    if (!v.readData(filename_data.c_str(), true))
    {
	    std::cerr << "unable to read file data.dat" << std::endl;
	    return 1;
    }
    v.saveData ("data.txt", true);
    // v.readData ("data.txt", true);
    std::cout << "data: " << std::endl;
    v.showData();
    //v.setSize(500,1000);
    if (!v.readNodes(filename_neurons.c_str(), true))
    {
	    std::cerr << "unable to read file nodes.dat" << std::endl;
	    return 1;
    }
    v.saveNodes ("nodes.txt", true);
    // v.readNodes ("nodes.txt", true);
    std::cout << "nodes: " << std::endl;
    v.showNodes ();
    v.getMaxMinValue ();
    v.calcVoronoiGnuplot();
    v.saveVoronoiGnuplot(filename_voronoi, "data.txt", "nodes.txt");
    v.setSizefromData(500);
    v.discretizeData ();
    v.discretizeNeurons ();
    v.calcVoronoiImage();
    v.saveVoronoiImage(output.c_str());
    std::cout <<"run done"<<std::endl;
    return 0;
}
