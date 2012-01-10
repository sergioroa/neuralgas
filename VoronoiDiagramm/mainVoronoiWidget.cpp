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
	    std::cerr << "unable to read file data.dat" << std::endl;
	    return 1;
    }
    if (!m.vw->voronoi->readNodes(filename_neurons.c_str()))
    {
	    std::cerr << "unable to read file nodes.dat" << std::endl;
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
