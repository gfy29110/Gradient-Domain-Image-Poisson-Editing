#include "ImageWidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>
#include "ChildWindow.h"

using namespace std;

ImageWidget::ImageWidget(ChildWindow *relatewindow)
{
	image_ = new QImage();
	image_backup_ = new QImage();

	draw_status_ = kNone;
	is_choosing_ = false;
	is_pasting_ = false;
	choose_status_ = kDefault;

	point_start_ = QPoint(0, 0);
	point_end_ = QPoint(0, 0);
	point_transform_ = QPoint(0, 0);

	polygon_ = new MyPolygon;
	rectangle_ = new MyRectangle;

	source_window_ = NULL;
}

ImageWidget::~ImageWidget(void)
{
}

int ImageWidget::ImageWidth()
{
	return image_->width();
}

int ImageWidget::ImageHeight()
{
	return image_->height();
}

void ImageWidget::set_draw_status_to_choose()
{
	draw_status_ = kChoose;
}

void ImageWidget::set_draw_status_to_paste()
{
	draw_status_ = kPaste;
}

QImage* ImageWidget::image()
{
	return image_;
}

void ImageWidget::set_source_window(ChildWindow* childwindow)
{
	source_window_ = childwindow;
}

void ImageWidget::paintEvent(QPaintEvent *paintevent)
{
	QPainter painter;
	painter.begin(this);

	// Draw background
	painter.setBrush(Qt::lightGray);
	QRect back_rect(0, 0, width(), height());
	painter.drawRect(back_rect);

	// Draw image
	QRect rect = QRect(0, 0, image_->width(), image_->height());
	painter.drawImage(rect, *image_);

	// Draw choose region
	if (choose_status_ == kRect) {
		rectangle_->RectPaint(painter);
	}
	if (choose_status_ == kPoly || choose_status_ == kFree) {
		polygon_->PolyPaint(painter);
	}

	painter.end();
	update();
}

void ImageWidget::mousePressEvent(QMouseEvent *mouseevent)
{
	if (Qt::LeftButton == mouseevent->button())
	{
		switch (draw_status_)
		{
		case kChoose:
			point_start_ = point_end_ = mouseevent->pos();

			switch (choose_status_)
			{
			case kRect:
				rectangle_->reset();
				rectangle_->set_start(point_start_.x(), point_start_.y());
				rectangle_->set_end(point_end_.x(), point_end_.y());
				break;
			case kPoly:
				if (!is_choosing_) {
					polygon_->reset();
					polygon_->SetActivePoint(point_end_);
					polygon_->AddPointIn(point_end_);
				}
				else {
					polygon_->SetActivePoint(point_end_);
					polygon_->AddPointIn(point_end_);
				}

				break;
			case kFree:
				polygon_->reset();
				polygon_->SetActivePoint(point_end_);
				polygon_->AddPointIn(point_end_);
				break;
			default:
				break;
			}
			is_choosing_ = true;

			update();
			break;

		case kPaste:
		{
			is_pasting_ = true;

			// Start point in object image
			int source_x, source_y, source_width, source_height;
			switch (choose_status_) {
			case kRect:
				source_x = rectangle_->get_start_x();
				source_y = rectangle_->get_start_y();
				source_width = rectangle_->get_width();
				source_height = rectangle_->get_height();
				break;
			case kPoly:
				source_x = polygon_->get_min_x();
				source_y = polygon_->get_min_y();
				source_width = polygon_->get_width();
				source_height = polygon_->get_height();
				break;
			case kFree:
				source_x = polygon_->get_min_x();
				source_y = polygon_->get_min_y();
				source_width = polygon_->get_width();
				source_height = polygon_->get_height();
				break;
			default:
				source_x = 0;
				source_y = 0;
				source_width = 0;
				source_height = 0;
				break;
			}


			// Start point in source image
			int xpos = mouseevent->pos().rx();
			int ypos = mouseevent->pos().ry();

			// Width and Height of rectangle region
			int w = source_window_->imagewidget_->image_->width();
			int h = source_window_->imagewidget_->image_->height();

			// Paste
			if (w - xpos > source_width&&h - ypos > source_height) {
				switch (choose_status_) {
				case kRect:
					rectangle_->set_translation(xpos, ypos);
					point_transform_ = QPoint(xpos, ypos);
					break;
				case kPoly:
					polygon_->set_translation(xpos, ypos);
					point_transform_ = QPoint(xpos, ypos);
					break;
				case kFree:
					polygon_->set_translation(xpos, ypos);
					point_transform_ = QPoint(xpos, ypos);
					break;
				default:
					break;
				}

			}

		}

		update();
		break;

		default:
			break;
		}
	}

	if (Qt::RightButton == mouseevent->button())
	{
		if (choose_status_ == kPoly && draw_status_ == kChoose) {
			is_choosing_ = false;
			draw_status_ = kNone;
			polygon_->set_is_over();
		}
	}
}

void ImageWidget::mouseMoveEvent(QMouseEvent *mouseevent)
{
	switch (draw_status_)
	{
	case kChoose:
		// Store point position for rectangle region
		if (is_choosing_)
		{
			point_end_ = mouseevent->pos();

			switch (choose_status_) {
			case kRect:
				rectangle_->set_end(point_end_.x(), point_end_.y());
				break;
			case kPoly:
				polygon_->SetActivePoint(point_end_);
				break;
			case kFree:
				polygon_->SetActivePoint(point_end_);
				polygon_->AddPointIn(point_end_);
				break;
			default:
				break;
			}
		}
		update();
		break;

	case kPaste:
		// Paste rectangle region to object image
		if (is_pasting_)
		{
			// Start point in object image
			int source_x, source_y, source_width, source_height;
			switch (choose_status_) {
			case kRect:
				source_x = rectangle_->get_start_x();
				source_y = rectangle_->get_start_y();
				source_width = rectangle_->get_width();
				source_height = rectangle_->get_height();
				break;
			case kPoly:
				source_x = polygon_->get_min_x();
				source_y = polygon_->get_min_y();
				source_width = polygon_->get_width();
				source_height = polygon_->get_height();
				break;
			case kFree:
				source_x = polygon_->get_min_x();
				source_y = polygon_->get_min_y();
				source_width = polygon_->get_width();
				source_height = polygon_->get_height();
				break;
			default:
				source_x = 0;
				source_y = 0;
				source_width = 0;
				source_height = 0;
				break;
			}


			// Start point in source image
			int xpos = mouseevent->pos().rx();
			int ypos = mouseevent->pos().ry();

			// Width and Height of rectangle region
			int w = source_window_->imagewidget_->image_->width();
			int h = source_window_->imagewidget_->image_->height();

			// Paste
			if (w - xpos > source_width&&h - ypos > source_height) {
				switch (choose_status_) {
				case kRect:
					rectangle_->set_translation(xpos, ypos);
					point_transform_ = QPoint(xpos, ypos);
					break;
				case kPoly:
					polygon_->set_translation(xpos, ypos);
					point_transform_ = QPoint(xpos, ypos);
					break;
				case kFree:
					polygon_->set_translation(xpos, ypos);
					point_transform_ = QPoint(xpos, ypos);
					break;
				default:
					break;
				}

			}
		}
		update();
		break;

	default:
		break;
	}

	update();
}

void ImageWidget::mouseReleaseEvent(QMouseEvent *mouseevent)
{
	switch (draw_status_)
	{
	case kChoose:
		if (is_choosing_)
		{
			point_end_ = mouseevent->pos();

			switch (choose_status_) {
			case kRect:
				rectangle_->set_end(point_end_.x(), point_end_.y());
				is_choosing_ = false;
				draw_status_ = kNone;
				break;
			case kPoly:
				break;
			case kFree:
				polygon_->set_is_over();
				is_choosing_ = false;
				draw_status_ = kNone;
				break;
			default:
				break;
			}
		}

	case kPaste:
		if (is_pasting_)
		{
			is_pasting_ = false;
			//draw_status_ = kNone;
		}

	default:
		break;
	}

	update();
}

void ImageWidget::Open(QString filename)
{
	// Load file
	if (!filename.isEmpty())
	{
		image_->load(filename);
		*(image_backup_) = *(image_);
	}

	//	setFixedSize(image_->width(), image_->height());
	//	relate_window_->setWindowFlags(Qt::Dialog);
	//	relate_window_->setFixedSize(QSize(image_->width(), image_->height()));
	//	relate_window_->setWindowFlags(Qt::SubWindow);

		//image_->invertPixels(QImage::InvertRgb);
		//*(image_) = image_->mirrored(true, true);
		//*(image_) = image_->rgbSwapped();
	cout << "image size: " << image_->width() << ' ' << image_->height() << endl;
	update();
}

void ImageWidget::Save()
{
	SaveAs();
}

void ImageWidget::SaveAs()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), ".", tr("Images(*.bmp *.png *.jpg)"));
	if (filename.isNull())
	{
		return;
	}

	image_->save(filename);
}

void ImageWidget::Invert()
{
	for (int i = 0; i < image_->width(); i++)
	{
		for (int j = 0; j < image_->height(); j++)
		{
			QRgb color = image_->pixel(i, j);
			image_->setPixel(i, j, qRgb(255 - qRed(color), 255 - qGreen(color), 255 - qBlue(color)));
		}
	}

	// equivalent member function of class QImage
	// image_->invertPixels(QImage::InvertRgb);
	update();
}

void ImageWidget::Mirror(bool ishorizontal, bool isvertical)
{
	QImage image_tmp(*(image_));
	int width = image_->width();
	int height = image_->height();

	if (ishorizontal)
	{
		if (isvertical)
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, height - 1 - j));
				}
			}
		}
		else
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					image_->setPixel(i, j, image_tmp.pixel(i, height - 1 - j));
				}
			}
		}

	}
	else
	{
		if (isvertical)
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, j));
				}
			}
		}
	}

	// equivalent member function of class QImage
	//*(image_) = image_->mirrored(true, true);
	update();
}

void ImageWidget::TurnGray()
{
	for (int i = 0; i < image_->width(); i++)
	{
		for (int j = 0; j < image_->height(); j++)
		{
			QRgb color = image_->pixel(i, j);
			int gray_value = (qRed(color) + qGreen(color) + qBlue(color)) / 3;
			image_->setPixel(i, j, qRgb(gray_value, gray_value, gray_value));
		}
	}

	update();
}

void ImageWidget::Restore()
{
	*(image_) = *(image_backup_);
	point_start_ = point_end_ = QPoint(0, 0);
	update();
}