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
    if (!v.readData(filename_data.c_str()))
    {
	    std::cerr << "unable to read file data.dat" << std::endl;
	    return 1;
    }
    v.saveData ("data.txt", true);
    // v.readData ("data.txt", true);
    std::cout << "data: " << std::endl;
    v.showData();
    //v.setSize(500,1000);
    if (!v.readNodes(filename_neurons.c_str()))
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
