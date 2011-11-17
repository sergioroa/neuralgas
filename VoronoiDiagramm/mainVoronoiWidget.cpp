#include <QtGui/QApplication>
#include "VoronoiWidget.h"

using namespace neuralgas;


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    VoronoiMainWindow m;

    char* filename_data=new char[255];
    char* filename_neurons=new char[255];
    strcpy(filename_data, "data.txt");
    strcpy(filename_neurons, "nodes.txt");
    if (!m.vw->voronoi->getData(filename_data))
    {
	    std::cerr << "unable to read file data.txt" << std::endl;
	    return 1;
    }
    if (!m.vw->voronoi->getNeurons(filename_neurons))
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
