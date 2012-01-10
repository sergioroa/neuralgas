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
* \file GaussianNoise.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef GAUSSIANNOISE_H
#define GAUSSIANNOISE_H

#include "DataGenerator.h"
#include <tools/math_helpers.h>
#include <string>

namespace neuralgas {

//! \class GaussianNoise
/*! \brief Generation of datasets by using Gaussian distributions with the possibility
     of adding white noise
 */
class GaussianNoise : public DataGenerator<double>
{
public:
	GaussianNoise();
	void generate(const int& number);
	virtual void reset();
	void setCanonicalDataset ();
	void setCustomizedDataset ();
	bool readCustomizedDataset (const char*);
	void setWhiteNoiseProb(const double&);
	void setLimits (double, double, double, double);
private:
	virtual  Vector<double>* generate();
	virtual  Vector<double>* next(){return generate();}
	std::vector<Vector<double> > centers_params;
	double min_x, min_y, max_x, max_y;
	double whitenoise_prob;

};

GaussianNoise::GaussianNoise() : DataGenerator<double>(2)
{
	srand( (unsigned)time( NULL ) );
	whitenoise_prob = -1.0;
}

void GaussianNoise::reset()
{
	DataGenerator<double>::reset();
}

void GaussianNoise::setCanonicalDataset ()
{
	setLimits (-3,3,-3,3);
	centers_params.resize (5);
	for (int i=0; i<5; i++)
		centers_params[i].resize (4);

	centers_params[0][0] = 0;
	centers_params[0][1] = 2;
	centers_params[0][2] = 0.1;
	centers_params[0][3] = 0.1;

	centers_params[1][0] = 0;
	centers_params[1][1] = 1;
	centers_params[1][2] = 0.1;
	centers_params[1][3] = 0.1;
	
	centers_params[2][0] = 2;
	centers_params[2][1] = 0;
	centers_params[2][2] = 0.3;
	centers_params[2][3] = 0.3;

	centers_params[3][0] = -1;
	centers_params[3][1] = 0;
	centers_params[3][2] = 0.2;
	centers_params[3][3] = 0.1;

	centers_params[4][0] = 0;
	centers_params[4][1] = -1;
	centers_params[4][2] = 0.1;
	centers_params[4][3] = 0.2;
}

void GaussianNoise::setCustomizedDataset ()
{
	int nr_centers;
	std::cout << "Type nr. of centers: ";
	std::cin >> nr_centers;

	centers_params.resize (nr_centers);
	for (int i=0; i<nr_centers; i++)
	{
		centers_params[i].resize (4);
		std::cout << "x = ";
		std::cin >> centers_params[i][0];
		std::cout << "y = ";
		std::cin >> centers_params[i][1];
		std::cout << "var x = ";
		std::cin >> centers_params[i][2];
		std::cout << "var y = ";
		std::cin >> centers_params[i][3];
	}

	double min_x, max_x, min_y, max_y;
	std::cout << "Min x: ";
	std::cin >> min_x;
	std::cout << "Max x: ";
	std::cin >> max_x;
	std::cout << "Min y: ";
	std::cin >> min_y;
	std::cout << "Max y: ";
	std::cin >> max_y;
	setLimits (min_x, max_x, min_y, max_y);
}

bool GaussianNoise::readCustomizedDataset (const char* filename)
{
	std::ifstream file (filename);
	if (!file)
		return false;
	int nr_centers;
	file >> nr_centers;
	centers_params.resize (nr_centers);
	for (int i=0; i< nr_centers; i++)
	{
		centers_params[i].resize (4);

		for (int j=0; j<4; j++)
		{
			file >> centers_params[i][j];
			std::cout << centers_params[i][j] << std::endl;
		}
		
	}

	double min_x, max_x, min_y, max_y;
	file >> min_x;
	file >> max_x;
	file >> min_y;
	file >> max_y;
	std::cout << min_x << ", " << max_x << ", " << min_y << ", " << max_y << std::endl;
	setLimits (min_x, max_x, min_y, max_y);

	file.close ();
	return true;

}


void GaussianNoise::setLimits (double min_x, double max_x, double min_y, double max_y)
{
	this->min_x = min_x;
	this->max_x = max_x;
	this->min_y = min_y;
	this->max_y = max_y;
	
}


void GaussianNoise::setWhiteNoiseProb(const double& prob)
{
	whitenoise_prob=prob;
}


void GaussianNoise::generate(const int& number)
{
	std::cout << "minx, maxx, miny, maxy: " << min_x << "," << max_x << "," << min_y << "," << max_y << std::endl;
	for (int i=0; i< number; i++)
	{
		double whitenoise_rand = double (rand() % 1000) / 1000;
		if (whitenoise_rand <= whitenoise_prob)
		{
			Vector<double>* whitenoise_point = new Vector<double>;
			whitenoise_point->at(0) = (double)((rand() / (static_cast<double>(RAND_MAX) + 1.0)) * (max_x - min_x) + min_x);
			whitenoise_point->push_back ((double)((rand() / (static_cast<double>(RAND_MAX) + 1.0)) * (max_y - min_y) + min_y));
			std::cout << "w: " << whitenoise_point->at(0) << ", " << whitenoise_point->at(1) << std::endl;
			_data->push_back (whitenoise_point);
		}
		else
			_data->push_back (generate());
	}
}

Vector<double>* GaussianNoise::generate()
{
	assert (centers_params.size() > 0);
	int center = rand() % centers_params.size();
	return normal_distribution (centers_params[center][0],
				    centers_params[center][1],
				    centers_params[center][2],
				    centers_params[center][3]);
}

} // namespace neuralgas

#endif
