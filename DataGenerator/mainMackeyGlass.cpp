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
#include "MackeyGlass.h"

using namespace std;
using namespace neuralgas;

int main(int argc, char *argv[])
{
    MackeyGlass mg;
    mg.setPastTimeSteps(17);
    mg.setBoundary(0.4);
    mg.setPower(10);
    int size = 100;
    mg.generate (size);
    showData (mg.getData());
    
    return EXIT_SUCCESS;
}
