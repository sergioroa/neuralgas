/** 
* \class GaussianNoise
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

namespace neuralgas {

class GaussianNoise : public DataGenerator<float>
{
public:
	GaussianNoise();
	void generate(const int& number);
	virtual void reset();
	void setCanonicalDataset ();
	void setCustomizedDataset () {}; // TO DO
	void setWhiteNoiseProb(const float&);
	void setLimits (float, float, float, float);
private:
	virtual  Vector<float>* generate();
	virtual  Vector<float>* next(){return generate();}
	std::vector<Vector<float> > centers_params;
	float min_x, min_y, max_x, max_y;
	float whitenoise_prob;

};

GaussianNoise::GaussianNoise() : DataGenerator<float>(2)
{
	srand( (unsigned)time( NULL ) );
	whitenoise_prob = -1.0;
}

void GaussianNoise::reset()
{
	DataGenerator<float>::reset();
}

void GaussianNoise::setCanonicalDataset ()
{
	setLimits (-5,5,-5,5);
	for (int i=0; i<5; i++)
	{
		Vector<float> params;
		centers_params.push_back (params);
		for (int j=0; j<4; j++)
			centers_params[i].push_back (0.0);
	}

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
	centers_params[2][2] = 0.1;
	centers_params[2][3] = 0.1;

	centers_params[3][0] = -1;
	centers_params[3][1] = 0;
	centers_params[3][2] = 0.2;
	centers_params[3][3] = 0.1;

	centers_params[4][0] = 0;
	centers_params[4][1] = -1;
	centers_params[4][2] = 0.1;
	centers_params[4][3] = 0.2;
}

void GaussianNoise::setLimits (float min_x, float max_x, float min_y, float max_y)
{
	this->min_x = min_x;
	this->max_x = max_x;
	this->min_y = min_y;
	this->max_y = max_y;
	
}


void GaussianNoise::setWhiteNoiseProb(const float& prob)
{
	whitenoise_prob=prob;
}


void GaussianNoise::generate(const int& number)
{
	std::cout << "minx, maxx, miny, maxy: " << min_x << "," << max_x << "," << min_y << "," << max_y << std::endl;
	for (int i=0; i< number; i++)
	{
		float whitenoise_rand = float (rand() % 1000) / 1000;
		if (whitenoise_rand <= whitenoise_prob)
		{
			Vector<float>* whitenoise_point = new Vector<float>;
			whitenoise_point->at(0) = (float)((rand() / (static_cast<float>(RAND_MAX) + 1.0)) * (max_x - min_x) + min_x);
			whitenoise_point->push_back ((float)((rand() / (static_cast<float>(RAND_MAX) + 1.0)) * (max_y - min_y) + min_y));
			std::cout << "w: " << whitenoise_point->at(0) << ", " << whitenoise_point->at(1) << std::endl;
			_data->push_back (whitenoise_point);
		}
		else
			_data->push_back (generate());
	}
}

Vector<float>* GaussianNoise::generate()
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