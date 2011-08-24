#include <QtGui/QApplication>
#include "VoronoiWidget.h"

using namespace neuralgas;


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow m;

    char* filename_data=new char[255];
    char* filename_neurons=new char[255];
    strcpy(filename_data, "data.txt");
    strcpy(filename_neurons, "nodes.txt");

    m.vw->voronoi->getData(filename_data);
    m.vw->voronoi->getNeurons(filename_neurons);
    m.vw->voronoi->setSizefromData(1600);
    m.vw->voronoi->calcVoronoi();

    m.vw->setImageSize ();
    m.show();

    return a.exec();
}
