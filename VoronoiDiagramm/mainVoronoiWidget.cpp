#include <QtGui/QApplication>
#include "VoronoiWidget.h"

using namespace neuralgas;


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    VoronoiMainWindow m;

    std::string filename_data = "data.dat";
    std::string filename_neurons = "nodes.dat";
    if (!m.vw->voronoi->readData(filename_data.c_str()))
    {
	    std::cerr << "unable to read file data.txt" << std::endl;
	    return 1;
    }
    if (!m.vw->voronoi->readNodes(filename_neurons.c_str()))
    {
	    std::cerr << "unable to read file nodes.txt" << std::endl;
	    return 1;
    }
    m.vw->voronoi->getMaxMinValue ();
    m.vw->voronoi->setSizefromData(1000);
    m.vw->voronoi->discretizeData();
    m.vw->voronoi->discretizeNeurons();
    m.vw->voronoi->calcVoronoiImage();

    m.vw->setImageSize ();
    m.resize (m.vw->width(), m.vw->height());
    m.show();
    
    return a.exec();
}
