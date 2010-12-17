#include "Voronoi.h"

void Voronoi::addData(const std::string& line)
{
    int i=0;
    while(line[i]!=' ')
    {
        i++;
    }
    point p;
    p.x = float (atof( (line.substr(0,i-1)).c_str() ));
    p.y = float (atof( (line.substr(i+1, line.size()-i - 2 )).c_str() ));
    _data.push_back(p);

}

void Voronoi::addNeuron(const std::string& line)
{
    int i=0;
    while(line[i]!=' ')
    {
        i++;
    }
    point p;
    p.x = float (atof( (line.substr(0,i-1)).c_str() ));
    p.y = float (atof( (line.substr(i+1, line.size()-i - 2 )).c_str() ));
    _neurons.push_back(p);
}

void Voronoi::getData(const char* filename)
{
    std::string line;
    std::ifstream myfile (filename);
    if (myfile.is_open())
    {
        while ( myfile.good() )
        {
          getline (myfile,line);
          addData(line);
        }
        myfile.close();
        getMaxMinValue();
    }
}

void Voronoi::showData()
{
    for(int i=0; i < _data.size(); i++)
        std::cout << _data[i].x << " " <<_data[i].y << std::endl;
}

void Voronoi::getMaxMinValue()
{
    _maxX=_data[0].x;
    _minX=_data[0].x;
    _maxY=_data[0].y;
    _minY=_data[0].y;

    for (int i=1; i < _data.size(); i++)
    {
        if (_data[i].x > _maxX)
            _maxX = _data[i].x;
        if (_data[i].x < _minX)
            _minX = _data[i].x;
        if (_data[i].y > _maxY)
            _maxY = _data[i].y;
        if (_data[i].y < _minY)
            _minY = _data[i].y;

    }
}

void Voronoi::getNeurons(const char* filename)
{
    std::string line;
    std::ifstream myfile (filename);
    if (myfile.is_open())
    {
        while ( myfile.good() )
        {
          getline (myfile,line);
          addNeuron(line);
        }
        myfile.close();
    }

}


void Voronoi::discretize()
{
    float rangeX = _maxX - _minX;
    float rangeY = _maxY - _minY;
 

    float factorX = _width / rangeX;
    float factorY = _height / rangeY;
    float minX = factorX * _minX;
    float minY = factorY * _minY;
    std::cout << "minX "<<_minX <<" rescaled minX "<<minX<<" minY "<<_minY<<" rescaled minY "<<minY<<std::endl;

    for(int i=0; i < _data.size(); i++)
    {
        _data[i].x=( _data[i].x * factorX - minX );
        _data[i].y=( _data[i].y * factorY - minY);
    }

    for(int i=0; i < _neurons.size(); i++)
    {
        _neurons[i].x =( _neurons[i].x * factorX - minX);
        _neurons[i].y =( _neurons[i].y * factorY - minY);
    }
}

void Voronoi::setSize(const int& height, const int& width)
{
    _height = height;
    _width = width;
}

void Voronoi::setSizefromData(const int& somesideSize)
{
    float rangeX = _maxX - _minX;
    float rangeY = _maxY - _minY;

    if (rangeX > rangeY) {
	    _width = somesideSize;
	    float factorX = _width / rangeX;
	    _height = ceil(factorX * rangeY);
    }
    else {
	    _height = somesideSize;
	    float factorY = _height / rangeY;
	    _width = ceil(factorY * rangeX);
    }

    std::cout << "width: " << _width << std::endl;
    std::cout << "height: " << _height << std::endl;


}

void Voronoi::setNeurons()
{
    _xValues = new float[_neurons.size()];
    _yValues = new float[_neurons.size()];

    for(int i=0; i < _neurons.size(); i++)
    {
        _xValues[i] = _neurons[i].x;
        _yValues[i] = _neurons[i].y;
    }
}

void Voronoi::calcVoronoi()
{
    discretize();
    setNeurons();
    _vdg.generateVoronoi(_xValues,_yValues,_neurons.size(), 0,_width,0,_height,0);
    _vdg.resetIterator();
}

void Voronoi::drawLine(float& x1, float& y1, float& x2, float& y2, QImage& image)
{
    if (x1>x2)
    {
        float tmp=x2;
        x2=x1;
        x1=tmp;
        tmp=y2;
        y2=y1;
        y1=tmp;
    }

    float m= (y2-y1)/(x2-x1);
    QRgb value;

    value = qRgb(237, 187, 51); // 0xffedba31

    int range = int(x2-x1);
    for(int i =0; i <= range; i++)
    {
            int slope = int(abs(m))+1;
            for (int j=0; j < slope; j++)
                if ( m < 0 && y1+m*i-j > y2 )
                    image.setPixel(x1+i, y1+m*i-j,value);
                else if ( m >= 0 && y1+m*i+j < y2 )
                    image.setPixel(x1+i, y1+m*i+j,value);
    }

}

void Voronoi::save(const char* filename)
{
    QImage image(_width, _height, QImage::Format_RGB32);

    QRgb value;
    //value = qRgb(255, 255, 255);
    //image.fill(value);
    image.fill(0);

    // data
    value = qRgb(255, 255, 255); // 0xff7aa327
    //value = qRgb(0, 0, 255);
    for(int i=0; i < _data.size();i++)
    {
        image.setPixel(_data[i].x, _data[i].y, value);
    }
    // neurons
    value = qRgb(255, 0, 0); // 0xffbd9527
    for(int i=0; i < _neurons.size();i++)
    {
        image.setPixel(_neurons[i].x, _neurons[i].y, value);
        image.setPixel(_neurons[i].x+1, _neurons[i].y, value);
        image.setPixel(_neurons[i].x-1, _neurons[i].y, value);
        image.setPixel(_neurons[i].x, _neurons[i].y+1, value);
        image.setPixel(_neurons[i].x, _neurons[i].y-1, value);
        image.setPixel(_neurons[i].x+1, _neurons[i].y+1, value);
        image.setPixel(_neurons[i].x+1, _neurons[i].y-1, value);
        image.setPixel(_neurons[i].x-1, _neurons[i].y+1, value);
        image.setPixel(_neurons[i].x-1, _neurons[i].y-1, value);
    }
    // voronoi lines
    float x1,y1,x2,y2;
    while(_vdg.getNext(x1,y1,x2,y2))
    {
        drawLine(x1,y1,x2,y2,image);
    }
    image.save(filename);
}
