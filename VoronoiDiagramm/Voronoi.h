#ifndef VORONOI_H
#define VORONOI_H

#include <Qt/qimage.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "VoronoiDiagramGenerator.h"

namespace neuralgas {

struct point
{
    float x;
    float y;
};

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
protected:
    void addNeuron(const std::string&);
    void addData(const std::string&);
    void getMaxMinValue();
    void discretize();
    void setNeurons();
    void drawLine(float&, float&, float&, float&, QImage&);
private:
    VoronoiDiagramGenerator _vdg;
    float* _xValues;
    float* _yValues;
    int _height, _width;
    float _maxX, _maxY, _minX, _minY;
    std::vector< point > _data;
    std::vector< point > _neurons;
};

} // namespace neuralgas

#endif // VORONOI_H
