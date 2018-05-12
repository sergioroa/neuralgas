/*
 *   This file is part of NeuralGas.
 *
 *   NeuralGas is free software: you can redistribute it and/or modify it
 *   under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   NeuralGas is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with NeuralGas.  If not, see <http://www.gnu.org/licenses/>.
 */

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

#include "VoronoiDiagramGenerator.h"
#include <QImage>
#include <QPainter>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Graphs/Base_Graph.h>

namespace neuralgas {

typedef std::vector < Vector<double>* > SeqData;
typedef std::vector < Base_Node<double, int>*, _NGPoolAllocdoubleint_ > SeqNeurons;

//! \class Voronoi
/*! \brief Encapsulates techniques for generating Voronoi Diagramms  
 */
class Voronoi
{
public:
    Voronoi ();
    ~Voronoi ();
    bool readData(const char*, bool t=false );
    bool readNodes(const char*, bool t = false );
    bool saveData(const char*, bool t = false);
    bool saveNodes(const char*, bool t = false);
    void calcVoronoiImage();
    void calcVoronoiGnuplot();
    void saveVoronoiImage(const char*);
    void saveVoronoiGnuplot (std::string, std::string, std::string);
    void setSize(const int&, const int&);
    void setSizefromData(const int&);
    void showData();
    void showNodes();
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
    void setNeurons();
    void drawDiagram (QPainter&);
    void resetData();
    void resetNodes();
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
    std::stringstream gnuplotstream;

};

} // namespace neuralgas

#endif // VORONOI_H
