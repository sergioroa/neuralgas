#ifndef VORONOI_H
#define VORONOI_H

#include <Qt/qimage.h>
#include <Qt/qpainter.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "VoronoiDiagramGenerator.h"
#include <Graphs/Vector.h>
#include <Graphs/Base_Graph.h>

namespace neuralgas {

class Voronoi
{
public:
    void getData(const char*);
    void getNeurons(const char*);
    void calcVoronoi();
    void save(const char*);
    void setSize(const int&, const int&);
    void setSizefromData(const int&);
    void showData();
    friend class VoronoiWidget;
protected:
    void addNeuron(const std::string&);
    void addData(const std::string&);
    void getMaxMinValue();
    void discretize();
    void setNeurons();
    void drawLine(double&, double&, double&, double&, QImage&);
private:
    VoronoiDiagramGenerator _vdg;
    double* _xValues;
    double* _yValues;
    int _height, _width;
    double _maxX, _maxY, _minX, _minY;
    std::vector < Vector<double>* >* _data; 
    std::vector < Base_Node<double, int>* >* _neurons;
};

} // namespace neuralgas

#endif // VORONOI_H
