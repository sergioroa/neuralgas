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

#include <DataGenerator/GaussianNoise.h>
#include "Voronoi.h"
#include <tools/metrics.h>

using namespace neuralgas;

double getMSError (std::vector<Vector<double> >& centers, SeqNeurons& nodes)
{
	double MSError = 0.0;
	std::vector<double> errors;
	for (unsigned int i=0; i<centers.size(); i++)
	{
		double nearest_distance = std::numeric_limits<double>::max();
		double distance;
		for (unsigned int j=0; j<nodes.size(); j++)
		{
			distance = euclidean<double, int> (centers[i], nodes[j]->weight);
			if (distance < nearest_distance)
			{
				nearest_distance = distance;
			}
		}
		errors.push_back (nearest_distance);
	}
	for (unsigned int i =0; i < errors.size(); i++)
		MSError += errors[i];
	MSError /= errors.size();

	return MSError;
}

int main(int argc, char *argv[])
{

    std::string filename_neurons="nodes.dat";
    std::string dataset="1";

    if (argc > 1)
    {
	    dataset = std::string(argv[1]);
    }
    GaussianNoise* gn = new GaussianNoise;
    if (dataset == "1")
	    gn->setCanonicalDataset ();
    else if (dataset == "2")
	    gn->setCustomizedDataset ();
    else
    {
	    const char* filename = dataset.c_str();
	    if (!gn->readCustomizedDataset (filename))
	    {
		    std::cerr << "Error opening file..." << std::endl;
		    return 1;
	    }
    }
    std::vector<Vector<double> > centers_params = gn->getCentersParams ();
    std::vector<Vector<double> > centers;
    for (unsigned int i=0; i<centers_params.size(); i++)
    {
	    Vector<double> newcenter;
	    newcenter[0] = centers_params[i][0];
	    newcenter.push_back (centers_params[i][1]);
	    centers.push_back (newcenter);
    }
    Voronoi v;
    if (!v.readNodes(filename_neurons.c_str()))
    {
	    std::cerr << "unable to read file nodes.dat" << std::endl;
	    return 1;
    }
    SeqNeurons* nodes = v.getNeurons ();

    double mserror = getMSError (centers, *nodes) ;

    std::cout << "nodes: " << std::endl;
    v.showNodes ();
    std::cout << "mserror: " << mserror << std::endl;
    // v.getMaxMinValue ();
    // v.calcVoronoiGnuplot();
    // v.saveVoronoiGnuplot(filename_voronoi, "data.txt", "nodes.txt");
    // v.setSizefromData(500);
    // v.discretizeData ();
    // v.discretizeNeurons ();
    // v.calcVoronoiImage();
    // v.saveVoronoiImage(output.c_str());
    std::cout <<"run done"<<std::endl;
    return 0;
}
