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

#include <cstdlib>
#include <iostream>
#include "TGraph.h"
#include "DGraph.h"
#include "UGraph.h"
#include "Base_Graph.h"

using namespace std;
using namespace neuralgas;

int main(int argc, char *argv[])
{

/*    Base_Graph<float,int> bgraph(10);
    bgraph.initRandomGraph(3);
    bgraph.addNode();
    bgraph.addNode();
    bgraph.showGraph();
    bgraph.rmNode(2);
    bgraph.showGraph();
 */
 /*   UGraph<float,int> ugraph(10);
    ugraph.initRandomGraph(3);
    ugraph.addNode();
    ugraph.addNode();
    ugraph.addEdge(0,1);
    ugraph.addEdge(0,2);
    ugraph.isConnected(1);
    ugraph.showGraph();
    ugraph.rmEdge(0,2);
    ugraph.showGraph();
*/
    DGraph<float,int> dgraph(10);
    dgraph.initRandomGraph(2, 0, 2);
    dgraph.addNode();
    dgraph.getID();
    dgraph.addEdge(0,1);
    dgraph.addNode();
    dgraph.addEdge(0,2);
    dgraph.showGraph();

 /*   TGraph<float,int> Tgraph(10);
    Tgraph.initRandomGraph(2);
    
    Tgraph.addEdge(0,1);
    Tgraph.addNode();
    Tgraph.addEdge(0,2);
//    Tgraph.getID(0,2);
//Tgraph.getID(0,1);
    Tgraph.incAge(0,1);
    Tgraph.incAge(0,1);
    Tgraph.incAge(0,1);
    std::cout << Tgraph.getAge(0,1)<<std::endl;
    Tgraph.decAge(0,1);
    std::cout << Tgraph.getAge(0,1)<<std::endl;
    Tgraph.setAge(0,1,2.4);
    std::cout << Tgraph.getAge(0,1)<<std::endl;

   // std::cout << dgraph[0].num_in_edges<<std::endl;*/
    return EXIT_SUCCESS;
}
