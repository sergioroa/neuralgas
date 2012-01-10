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

VoronoiWidget::~VoronoiWidget ()
{
	delete voronoi;
}

void VoronoiWidget::paintEvent (QPaintEvent *event)
{
	// call setImageSize first
	assert (!image.isNull());

	QPainter painter (this);
	//image.fill(value);
	//image.fill(0);
	// image.fill (Qt::transparent);

	voronoi->drawDiagram (painter);

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
	delete vw;
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
	vw->voronoi->calcVoronoiImage ();
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
	vw->voronoi->calcVoronoiImage ();
	vw->setImageSize ();
	resize (vw->width()+10, vw->height()+10);

}


} // namespace neuralgas
