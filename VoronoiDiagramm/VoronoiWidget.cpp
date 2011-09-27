/*!
  \file   VoronoiWidget.cpp
  \author Sergio Roa
  \date   Fri Aug 19 19:56:42 2011
  
  \brief  
  
  
*/
#include "VoronoiWidget.h"
#include <cassert>

namespace neuralgas
{

VoronoiWidget::VoronoiWidget (QWidget *parent) : QWidget (parent)
{
	//setBackgroundColor(QColor (0,0,0));
	 QPalette palette;
	 palette.setColor(backgroundRole(), QColor(0,0,0));
	 setPalette(palette);

	 voronoi = new Voronoi;
}


void VoronoiWidget::paintEvent (QPaintEvent *event)
{
	// call setImageSize first
	assert (!image.isNull());
	
	QPainter painter (this);
	QRgb value;
	// value = qRgb(255, 255, 255);
	//image.fill(value);
	//image.fill(0);
	// image.fill (Qt::transparent);

	// data
	value = qRgb(255, 255, 255); // 0xff7aa327
	//value = qRgb(0, 0, 0);
	painter.setPen (value);
	painter.setBrush (QBrush(value));
	//value = qRgb(0, 0, 255);
	for(unsigned int i=0; i < voronoi->_data->size();i++)
		painter.drawEllipse (voronoi->_data->at(i)->at(0), voronoi->_data->at(i)->at(1), 2, 2);
	// neurons
	value = qRgb(255, 0, 0); // 0xffbd9527
	painter.setPen (value);
	painter.setBrush (QBrush(value));
	for(unsigned int i=0; i < voronoi->_neurons->size();i++)
		painter.drawEllipse (voronoi->_neurons->at(i)->weight[0], voronoi->_neurons->at(i)->weight[1], 5, 5);

	// voronoi lines
	//value = qRgb(237, 187, 51); // 0xffedba31
	value = qRgb(0, 0, 255); // 0xffedba31
	painter.setPen (value);
	
	double x1,y1,x2,y2;
	while(voronoi->_vdg.getNext(x1,y1,x2,y2))
		painter.drawLine (x1, y1, x2, y2);
	voronoi->_vdg.resetIterator();


}

void VoronoiWidget::setImageSize ()
{
	image = QPixmap (voronoi->_width, voronoi->_height);
	image.fill (Qt::transparent);
	resize(voronoi->_width, voronoi->_height);

}

VoronoiMainWindow::VoronoiMainWindow(QWidget *parent)
	: QMainWindow(parent)
{

	vw = new VoronoiWidget;
	QScrollArea* scroll(new QScrollArea);
	QHBoxLayout* layout(new QHBoxLayout(scroll)); 
	vw->setSizePolicy( QSizePolicy::Expanding, QSizePolicy::Expanding );
	layout->addWidget(vw);
	scroll->setWidget (vw);
	setCentralWidget( scroll );
   
}

VoronoiMainWindow::~VoronoiMainWindow()
{
}

void VoronoiMainWindow::customEvent(QEvent* e) {
	if (e->type() == 1001) {
		std::cout << "resizing..." << std::endl;
		ResizeEvent* rev = dynamic_cast<ResizeEvent*>(e);
		resize (rev->width, rev->height);
	}
	else if (e->type() == 1002) {
		std::cout << "showing..." << std::endl;
		show ();
	}
	
}

void VoronoiMainWindow::updateData ( SeqNeurons* neurons) {
	mutex->lock ();
	std::cout << "updating..." << std::endl;
	vw->voronoi->setNeurons (neurons);
	vw->voronoi->discretizeNeurons ();
	vw->voronoi->calcVoronoi ();
	vw->repaint ();
	condition->wakeAll ();
	mutex->unlock ();
	
}

void VoronoiMainWindow::initializeData ( SeqData* data, SeqNeurons* neurons, unsigned int sidesize)
{
	vw->voronoi->setData (data);
	vw->voronoi->getMaxMinValue();
	vw->voronoi->setSizefromData(sidesize);
	vw->voronoi->discretizeData ();
	vw->voronoi->setNeurons (neurons);
	vw->voronoi->discretizeNeurons ();
	vw->voronoi->calcVoronoi ();
	vw->setImageSize ();
	resize (vw->width()+10, vw->height()+10);

}


} // namespace neuralgas
