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
#include <QScrollArea>
#include <QVBoxLayout>
#include <QtGui/QMainWindow>
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
	//! The widget to print Voronoi diagrams
	VoronoiWidget* vw;

private:
};




} // namespace neuralgas

#endif
