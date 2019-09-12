#pragma once
#include <QWidget>
#include "Polygon.h"
#include "Rectangle.h"

class ChildWindow;
QT_BEGIN_NAMESPACE
class QImage;
class QPainter;
QT_END_NAMESPACE

enum DrawStatus
{
	kChoose,
	kPaste,
	kNone
};

enum ChooseStatus
{
	kDefault,
	kRect,
	kPoly,
	kFree,
	kTest
};

class ImageWidget :
	public QWidget
{
	Q_OBJECT

public:
	ImageWidget(ChildWindow *relatewindow);
	~ImageWidget(void);

	int ImageWidth();											// Width of image
	int ImageHeight();											// Height of image
	void set_draw_status_to_choose();
	void set_draw_status_to_paste();
	QImage* image();
	void set_source_window(ChildWindow* childwindow);

protected:
	void paintEvent(QPaintEvent *paintevent);
	void mousePressEvent(QMouseEvent *mouseevent);
	void mouseMoveEvent(QMouseEvent *mouseevent);
	void mouseReleaseEvent(QMouseEvent *mouseevent);

public slots:
	// File IO
	void Open(QString filename);								// Open an image file, support ".bmp, .png, .jpg" format
	void Save();												// Save image to current file
	void SaveAs();												// Save image to another file

	// Image processing
	void Invert();												// Invert pixel value in image
	void Mirror(bool horizontal = false, bool vertical = true);		// Mirror image vertically or horizontally
	void TurnGray();											// Turn image to gray-scale map
	void Restore();												// Restore image to origin

public:
	void SetChooseRectangle() { choose_status_ = kRect; rectangle_->reset(); }
	void SetChoosePolygon() { choose_status_ = kPoly; polygon_->reset(); }
	void SetChooseFree() { choose_status_ = kFree; polygon_->reset(); }
	void SetChooseforTest() { choose_status_ = kTest; }
	ChooseStatus GetChooseStatus() { return choose_status_; }
	void SetPaste(ChooseStatus cs) { draw_status_ = kPaste; choose_status_ = cs; }
	void SetDrawDefault() { draw_status_ = kNone; }
	void SetChooseDefault() { choose_status_ = kDefault; }


public:
	QPoint						point_start_;					// Left top point of rectangle region
	QPoint						point_end_;						// Right bottom point of rectangle region
	QPoint						point_transform_;				// Transform point

private:
	QImage						*image_;						// image 
	QImage						*image_backup_;

	// Pointer of child window
	ChildWindow					*source_window_;				// Source child window

	// Signs
	DrawStatus					draw_status_;					// Enum type of draw status
	bool						is_choosing_;
	bool						is_pasting_;

	ChooseStatus				choose_status_;

public:
	MyPolygon					*polygon_;
	MyRectangle					*rectangle_;
};

