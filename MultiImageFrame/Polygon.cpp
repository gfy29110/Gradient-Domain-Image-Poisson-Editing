#include "Polygon.h"


MyPolygon::MyPolygon()
{
	vertexes = new vector<QPoint>;
	is_over = false;
	//is_start = false;
	x_max = 0;
	x_min = 0;
	y_max = 0;
	y_min = 0;
	translation_x = 0;
	translation_y = 0;
	active_point = QPoint(0, 0);

	is_bounday = false;
}


MyPolygon::~MyPolygon()
{
	vertexes->clear();
}

void MyPolygon::set_is_over() {
	is_over = true;
	is_bounday = true;
}

void MyPolygon::set_not_over() {
	is_over = false;
}

/*void MyPolygon::set_is_start() {
	is_start = true;
}

void MyPolygon::set_not_start() {
	is_start = false;
}*/

void MyPolygon::AddPointIn(QPoint p) {
	vertexes->push_back(p);
	if (vertexes->size() == 1) {
		x_max = p.rx();
		x_min = p.rx();
		y_max = p.ry();
		y_min = p.ry();
	}
	else {
		x_max = p.rx() > x_max ? p.rx() : x_max;
		x_min = x_min > p.rx() ? p.rx() : x_min;
		y_max = p.ry() > y_max ? p.ry() : y_max;
		y_min = y_min > p.ry() ? p.ry() : y_min;
	}
}

void MyPolygon::PolyPaint(QPainter &paint) {
	if (vertexes->size() > 1) {
		paint.setPen(poly_pen);
		for (int i = 0; i < vertexes->size() - 1; i++) {
			paint.drawLine(vertexes->at(i).x() + translation_x, vertexes->at(i).y() + translation_y, vertexes->at(i + 1).x() + translation_x, vertexes->at(i + 1).y() + translation_y);
		}
		paint.setPen(dot_pen);
		for (int i = 0; i < vertexes->size(); i++) {
			paint.drawPoint(vertexes->at(i).x() + translation_x, vertexes->at(i).y() + translation_y);
		}
		if (is_over) {
			paint.setPen(poly_pen);
			paint.drawLine(vertexes->at(0).x() + translation_x, vertexes->at(0).y() + translation_y, vertexes->at(vertexes->size() - 1).x() + translation_x, vertexes->at(vertexes->size() - 1).y() + translation_y);
			paint.setPen(dot_pen);
			if (is_bounday) {
				paint.setPen(rect_bound_pen);
				//paint.drawRect(x_min + translation_x, y_min + translation_y, x_max - x_min, y_max - y_min);
				paint.drawLine(x_min + translation_x, y_min + translation_y, x_min + translation_x, y_max + translation_y);
				paint.drawLine(x_min + translation_x, y_min + translation_y, x_max + translation_x, y_min + translation_y);
				paint.drawLine(x_max + translation_x, y_max + translation_y, x_min + translation_x, y_max + translation_y);
				paint.drawLine(x_max + translation_x, y_max + translation_y, x_max + translation_x, y_min + translation_y);
				paint.setPen(dot_pen);
				paint.drawPoint(x_min + translation_x, y_min + translation_y);
				paint.drawPoint(x_max + translation_x, y_max + translation_y);
				paint.drawPoint(x_min + translation_x, y_max + translation_y);
				paint.drawPoint(x_max + translation_x, y_min + translation_y);
			}
			else {
				paint.setPen(poly_pen);
				paint.drawLine(vertexes->at(vertexes->size() - 1).x() + translation_x, vertexes->at(vertexes->size() - 1).y() + translation_y, active_point.x() + translation_x, active_point.y() + translation_y);
			}
		}
	}
	else {
		if (vertexes->size() == 1&&!is_over) {
			paint.setPen(dot_pen);
			paint.drawPoint(vertexes->at(0).x() + translation_x, vertexes->at(0).y() + translation_y);
			paint.setPen(poly_pen);
			paint.drawLine(vertexes->at(vertexes->size() - 1).x() + translation_x, vertexes->at(vertexes->size() - 1).y() + translation_y, active_point.x() + translation_x, active_point.y() + translation_y);
		}
	}
}

void MyPolygon::set_translation(int x, int y)
{
	translation_x = x - x_min + 1;
	translation_y = y - y_min + 1;
}

void MyPolygon::reset() {
	is_over = false;
	//is_start = false;
	vertexes->clear();
	x_max = 0;
	x_min = 0;
	y_max = 0;
	y_min = 0;
}

void MyPolygon::paste(MyPolygon * poly)
{
	reset();

	for (int i = 0; i < poly->vertexes->size(); i++) {
		vertexes->push_back(poly->vertexes->at(i));
	}

	is_over = poly->is_over;
	x_max = poly->get_max_x();
	x_min = poly->get_min_x();
	y_max = poly->get_max_y();
	y_min = poly->get_min_y();
	translation_x = translation_y = 0;

	is_bounday = true;
}
