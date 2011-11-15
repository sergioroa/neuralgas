/** 
* \file Voronoi.h
* \author Manuel Noll
* \author Sergio Roa
* 
*  \version 1.0
*  \date    2010
*/

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

typedef std::vector < Vector<double>* > SeqData;
typedef std::vector < Base_Node<double, int>* > SeqNeurons;

//! \class Voronoi
/*! \brief Encapsulates techniques for generating Voronoi Diagramms  
 */
class Voronoi
{
public:
    Voronoi ();
    ~Voronoi ();
    bool getData(const char*);
    bool getNeurons(const char*);
    void calcVoronoi();
    void save(const char*);
    void setSize(const int&, const int&);
    void setSizefromData(const int&);
    void showData();
    SeqData* getData () { return _data; }
    void setData (SeqData* d);
    void setNeurons (SeqNeurons* neurons);
    SeqNeurons* getNeurons () { return _neurons; }
    void getMaxMinValue();
    void discretizeData();
    void discretizeNeurons();
    void setWhiteBackground();
    friend class VoronoiWidget;
    friend class VoronoiMainWindow;
protected:
    void addNeuron(const std::string&);
    void addData(const std::string&);
    void setNeurons();
    void drawDiagram (QPainter&);
private:
    VoronoiDiagramGenerator _vdg;
    double* _xValues;
    double* _yValues;
    int _height, _width;
    double _maxX, _maxY, _minX, _minY, minX, minY;
    double rangeX, rangeY, factorX, factorY;
    SeqData* _data; 
    SeqNeurons* _neurons;
    QRgb dataColor;
    QRgb backgroundColor;
};

} // namespace neuralgas

#endif // VORONOI_H
