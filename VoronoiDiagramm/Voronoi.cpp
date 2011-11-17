#include "Voronoi.h"

namespace neuralgas {

Voronoi::Voronoi ()
{
    _data = NULL;
    _neurons = NULL;
    _xValues = NULL;
    _yValues = NULL;
    // data color
    dataColor = qRgb(255, 255, 255);
    backgroundColor = Qt::transparent;
}

Voronoi::~Voronoi ()
{
    if (_data != NULL)
    {
	for (unsigned int i=0; i < _data->size(); i++)
	    delete (*_data)[i];
	_data->clear();
	delete _data;
    }

    if (_neurons != NULL)
    {
	for (unsigned int i=0; i<_neurons->size(); i++)
	    for (unsigned int j = i+1; j < _neurons->at(i)->edges.size(); j++)
		if (_neurons->at(i)->edges[j] != NULL)
		    delete _neurons->at(i)->edges[j];

	for (unsigned int i=0; i<_neurons->size(); i++)
	    _neurons->at(i)->edges.clear();

	for (unsigned int i=0; i<_neurons->size(); i++)
	    delete _neurons->at(i);
	_neurons->clear();
	delete _neurons;
    }

    if (_xValues != NULL)
	delete _xValues;
    if (_yValues != NULL)
	delete _yValues;
    
}

void Voronoi::addData(const std::string& line)
{
    int i=0;
    while(line[i]!=' ')
    {
        i++;
    }
    Vector<double>* point = new Vector<double>(2);
    point->at(0) = atof( (line.substr(0,i)).c_str());
    point->at(1) = atof( (line.substr(i+1, line.size()-i - 2 )).c_str() );
    _data->push_back(point);
    std::cout << point->at(0) << " " << point->at(1) << std::endl;

}

void Voronoi::addNeuron(const std::string& line)
{
    int i=0;
    while(line[i]!=' ')
    {
        i++;
    }
    Base_Node<double, int>* neuron = new Base_Node<double, int>;
    neuron->weight.resize (2);
    neuron->weight[0] = atof( (line.substr(0,i)).c_str() );
    neuron->weight[1] = atof( (line.substr(i+1, line.size()-i - 2 )).c_str() );
    _neurons->push_back(neuron);
}

bool Voronoi::getData(const char* filename)
{
    _data = new std::vector < Vector<double>* >;
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
        // getMaxMinValue();
	return true;
    }
    return false;
}

void Voronoi::showData()
{
    for(unsigned int i=0; i < _data->size(); i++)
        std::cout << (*(*_data)[i])[0] << " " << (*(*_data)[i])[1] << std::endl;
}

void Voronoi::getMaxMinValue()
{
    _maxX=(*(*_data)[0])[0];
    _minX=(*(*_data)[0])[0];
    _maxY=(*(*_data)[0])[1];
    _minY=(*(*_data)[0])[1];
    
    for (unsigned int i=1; i < _data->size(); i++)
    {

        if ((*(*_data)[i])[0] > _maxX)
            _maxX = (*(*_data)[i])[0];
        if ((*(*_data)[i])[0] < _minX)
            _minX = (*(*_data)[i])[0];
        if ((*(*_data)[i])[1] > _maxY)
            _maxY = (*(*_data)[i])[1];
        if ((*(*_data)[i])[1] < _minY)
            _minY = (*(*_data)[i])[1];

	
    }
}

bool Voronoi::getNeurons(const char* filename)
{
    _neurons = new std::vector < Base_Node<double, int>* >;
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
	return true;
    }
    return false;

}

void Voronoi::discretizeData()
{

    for(unsigned int i=0; i < _data->size(); i++)
    {
        (*(*_data)[i])[0]=( (*(*_data)[i])[0] * factorX - minX );
	(*(*_data)[i])[1]=( (*(*_data)[i])[1] * factorY - minY);
    }

}

void Voronoi::discretizeNeurons()
{
    for(unsigned int i=0; i < _neurons->size(); i++)
    {
        (*_neurons)[i]->weight[0] =( (*_neurons)[i]->weight[0] * factorX - minX);
        (*_neurons)[i]->weight[1] =( (*_neurons)[i]->weight[1] * factorY - minY);
    }
}


void Voronoi::setSize(const int& height, const int& width)
{
    _height = height;
    _width = width;

    rangeX = _maxX - _minX;
    rangeY = _maxY - _minY;
 
    factorX = _width / rangeX;
    factorY = _height / rangeY;

    minX = factorX * _minX;
    minY = factorY * _minY;

    std::cout << "minX "<<_minX <<" rescaled minX "<<minX<<" minY "<<_minY<<" rescaled minY "<<minY<<std::endl;

}

void Voronoi::setSizefromData(const int& somesideSize)
{
    rangeX = _maxX - _minX;
    rangeY = _maxY - _minY;

    if (rangeX > rangeY) {
	    _width = somesideSize;
	    factorX = _width / rangeX;
	    _height = ceil(factorX * rangeY);
	    factorY = _height / rangeY;
    }
    else {
	    _height = somesideSize;
	    factorY = _height / rangeY;
	    _width = ceil(factorY * rangeX);
	    factorX = _width / rangeX;
    }

    std::cout << "width: " << _width << std::endl;
    std::cout << "height: " << _height << std::endl;

    minX = factorX * _minX;
    minY = factorY * _minY;

    std::cout << "minX "<<_minX <<" rescaled minX "<<minX<<" minY "<<_minY<<" rescaled minY "<<minY<<std::endl;

}

void Voronoi::setNeurons()
{
    _xValues = new double[_neurons->size()];
    _yValues = new double[_neurons->size()];

    for(unsigned int i=0; i < _neurons->size(); i++)
    {
        _xValues[i] = (*_neurons)[i]->weight[0];
        _yValues[i] = (*_neurons)[i]->weight[1];
    }
}

void Voronoi::calcVoronoiImage()
{
    // discretize();
    setNeurons();
    _vdg.generateVoronoi(_xValues,_yValues,_neurons->size(), 0,_width,0,_height,0);
    _vdg.resetIterator();
}

void Voronoi::calcVoronoiGnuplot()
{
    // discretize();
    setNeurons();
    _vdg.generateVoronoi(_xValues,_yValues,_neurons->size(), _minX,_maxX,_minY,_maxY,0);
    _vdg.resetIterator();
}


void Voronoi::drawDiagram (QPainter& painter)
{

    QMatrix m;
    m.translate( 0, _height );
    m.scale( 1, -1 );
    painter.setMatrix (m);
     
    // data
    painter.setPen (dataColor);
    painter.setBrush (QBrush(dataColor));

    for(unsigned int i=0; i < _data->size();i++)
	painter.drawEllipse (QPoint((*(*_data)[i])[0], (*(*_data)[i])[1]), 1, 1);
	    

    // neurons
    QRgb value = qRgb(255, 0, 0);
    painter.setPen (value);
    painter.setBrush (QBrush(value));
    for(unsigned int i=0; i < _neurons->size();i++)
	painter.drawEllipse (QPoint((*_neurons)[i]->weight[0], (*_neurons)[i]->weight[1]), 2, 2);

    // voronoi lines
    //value = qRgb(237, 187, 51); // 0xffedba31
    value = qRgb(0, 0, 255);
    painter.setPen (value);

    double x1,y1,x2,y2;
    while(_vdg.getNext(x1,y1,x2,y2))
	painter.drawLine (x1, y1, x2, y2);
    
}

void Voronoi::saveVoronoiImage(const char* filename)
{
    QImage image(_width, _height, QImage::Format_ARGB32);
    image.fill(backgroundColor);
    QPainter painter(&image);

    drawDiagram (painter);

    image.save(filename, "PNG");
}

void Voronoi::saveVoronoiGnuplot (std::string filename, std::string datafilename, std::string nodesfilename)
{
    gnuplotstream << "set terminal postscript eps enhanced color font \"Times-Roman,14\"" << std::endl;
    // gnuplotstream << "set output \"|epstopdf --filter > '" + fileName + ".pdf'" << endl;
    gnuplotstream << "set output \"" + filename + ".eps" << std::endl;
    gnuplotstream << "set autoscale" << std::endl;
    gnuplotstream << "unset key" << std::endl;
    gnuplotstream << "set style line 1 lt 1 lc rgb \"black\"" << std::endl;
    gnuplotstream << "set style line 2 lt 1 lw 2 lc rgb \"red\"" << std::endl;
    gnuplotstream << "set style line 3 lt 1 lc rgb \"blue\"" << std::endl;
    
    double x1,y1,x2,y2;
    while(_vdg.getNext(x1,y1,x2,y2))
	gnuplotstream << "set arrow from " << x1 << "," << y1 << " to " << x2 << "," << y2 << " nohead ls 3" << std::endl;

    gnuplotstream << "plot \"" << datafilename << "\" ls 1, \\" << std::endl; 
    gnuplotstream << "\"" << nodesfilename << "\" ls 2" << std::endl;

    std::string gnuplotfile = filename + ".gnu";
    std::ofstream StartPosgnuplotScript(gnuplotfile.c_str(), std::ios::out);
    StartPosgnuplotScript << gnuplotstream.str() << std::endl;
    StartPosgnuplotScript.close();
    
}

void Voronoi::setWhiteBackground ()
{
    dataColor = qRgb(0, 0, 0);
    backgroundColor = qRgb(255, 255, 255);
}

void Voronoi::setData (SeqData* d)
{
    _data = new SeqData ();
    for (unsigned int i=0; i<d->size(); i++)
    {
	Vector<double>* point = new Vector<double>(*d->at(i));
	_data->push_back (point);
    }
}

void Voronoi::setNeurons (SeqNeurons* n)
{
    _neurons = new SeqNeurons ();
    for (unsigned int i=0; i<n->size(); i++)
    {
	Base_Node<double, int>* neuron = new Base_Node<double, int>();
	neuron->weight = n->at(i)->weight;
	for(unsigned int k=0; k < n->size(); k++)
	    neuron->edges.push_back(NULL);

	_neurons->push_back (neuron);
    }

    for (unsigned int i=0; i<n->size(); i++)
    	for (unsigned int j=0; j<n->at(i)->edges.size(); j++)
    	    if ( n->at(i)->edges[j] != NULL )
		if (_neurons->at(i)->edges[j] == NULL)
		{
		    Base_Edge<int, double>* new_edge = new Base_Edge <int, double>;
		    new_edge->in = _neurons->at(i);
		    new_edge->out = _neurons->at(j);
		    _neurons->at(i)->edges[j] = new_edge;
		    _neurons->at(j)->edges[i] = new_edge;
		    _neurons->at(i)->num_connections++;
		    _neurons->at(j)->num_connections++;
		}
}


} // namespace neuralgas
