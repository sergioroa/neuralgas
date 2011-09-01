/** 
* \class VoronoiWidget
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

namespace neuralgas {

class VoronoiWidget : public QWidget
{
	Q_OBJECT
public:
	//! \brief Constructor
	/*! 
	  
	  \param parent parent widget	  
	*/
	VoronoiWidget (QWidget *parent=0);
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

class ShowEvent : public QEvent
{
public:
	ShowEvent() : QEvent ((QEvent::Type)1002) { }

};

class UpdateDataEvent : public QEvent
{
public:
	int region;
	//const SeqData data;
	SeqNeurons* neurons;
	UpdateDataEvent(/*const SeqData& d, */SeqNeurons* n) :
		QEvent ((QEvent::Type)1003),
		// data (d),
		neurons (n)
	{
	}

};


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

	void updateData ( /*const SeqData& data, */SeqNeurons* neurons) {
		QApplication::postEvent(this, new UpdateDataEvent(/*data,*/ neurons) );

	}

	
	//! The widget to print Voronoi diagrams
	VoronoiWidget* vw;

private:
};



} // namespace neuralgas

#endif
