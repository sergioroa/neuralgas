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
* \file VoronoiWidget.h
* \author Sergio Roa
* 
*  Copyright(c) 2011 Sergio Roa - All rights reserved
*  \version 1.0
*  \date    2011
*/

#ifndef VORONOIWIDGET_H
#define VORONOIWIDGET_H

#include <QWidget>
#include <QPixmap>
#include <QPainter>
#include <QEvent>
#include <QScrollArea>
#include <QVBoxLayout>
#include <QtGui/QMainWindow>
#include <QtGui/QApplication>
#include <VoronoiDiagramm/Voronoi.h>
#include <QMutex>
#include <QWaitCondition>

namespace neuralgas {


//! \class VoronoiWidget
/*! \brief Qt Widget for painting Voronoi Diagramms  
 */
class VoronoiWidget : public QWidget
{
	Q_OBJECT
public:
	//! \brief Constructor
	/*! 
	  
	  \param parent parent widget	  
	*/
	VoronoiWidget (QWidget *parent=0);
	//! \brief Destructor
	/*! 
	  
	  \param parent parent widget	  
	*/
	~VoronoiWidget ();
	/** 
	 * \brief paint the \p image
	 * 
	 * \param event paint event
	 */
	void paintEvent (QPaintEvent *event);

	/** 
	 * \brief initialize \p image size
	 * 
	 */
	void setImageSize ();

	//! data
	Voronoi *voronoi;

protected:
	//! image to paint
	QPixmap image;

};


//! \class ResizeEvent
/*! \brief Event class triggered when resizing the window
 */
class ResizeEvent : public QEvent
{
public:

	int width;
	int height;
	ResizeEvent(int w, int h) :
		QEvent ((QEvent::Type)1001),
		width(w),
		height(h)
	{
	}
	
};

//! \class ShowEvent
/*! \brief Event class triggered w showing the window
 */
class ShowEvent : public QEvent
{
public:
	ShowEvent() : QEvent ((QEvent::Type)1002) { }

};


//! \class VoronoiMainWindow
/*! \brief a Qt Main Window with a VoronoiWidget that handles events  
 */
class VoronoiMainWindow : public QMainWindow
{
	Q_OBJECT
	
public:
	//! \brief constructor
	/*! 
	  \param parent QWidget parent 
	*/
	VoronoiMainWindow(QWidget *parent = 0);
	
	//! \brief destructor
	~VoronoiMainWindow();
	
	//! \brief handle events like resize, update data, show
	/*! 
	  \param e Event
	*/
	virtual void customEvent(QEvent* e);

	//! The widget to print Voronoi diagrams
	VoronoiWidget* vw;
	//! \brief set Mutual exclusion object pointer
	/*! 
	  \param m mutex pointer
	*/
	void setMutex (QMutex* m) { mutex = m; }
	//! \brief set Wait condition object pointer
	/*! 
	  \param c waitcondition pointer
	*/
	void setWaitCondition (QWaitCondition *c) { condition = c; }

public slots:
	//! \brief update nodes
	/*! 
	  \param neurons or nodes
	*/
	void updateData (SeqNeurons* neurons);
	//! \brief initialize data, nodes and size of the painting area
	/*! 
	  \param data data set vectors
	  \param neurons nodes
	  \param sidesize size of painting area
	*/
	void initializeData ( SeqData* data, SeqNeurons* neurons, unsigned int sidesize);
	
protected:
	//! Mutual exclusion pointer for threads
	QMutex *mutex;
	//! Wait condition pointer for implementing semaphore
	QWaitCondition *condition;

};



} // namespace neuralgas

#endif
