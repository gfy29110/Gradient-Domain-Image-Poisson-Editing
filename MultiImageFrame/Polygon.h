#pragma once
#include <QPainter>
#include <qpen.h>
#include <vector>

using namespace std;

class MyPolygon
{
	const QPen rect_bound_pen = QPen(Qt::black, 1, Qt::DashLine, Qt::SquareCap, Qt::RoundJoin);
	const QPen dot_pen = QPen(Qt::black, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	const QPen poly_pen = QPen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);

public:
	MyPolygon();
	~MyPolygon();
	void AddPointIn(QPoint p);
	void PolyPaint(QPainter &paint);
	void set_translation(int x, int y);
	void set_is_over();
	void set_not_over();
	/*void set_is_start();
	void set_not_start();*/
	int get_max_x() { return x_max; }
	int get_min_x() { return x_min; }
	int get_max_y() { return y_max; }
	int get_min_y() { return y_min; }
	int get_width() { return x_max - x_min; }
	int get_height() { return y_max - y_min; }
	vector<QPoint> *vertexes;
	void reset();
	bool get_is_over() { return is_over; }
	//bool get_is_start() { return is_start; }

	void paste(MyPolygon *poly);
	void SetActivePoint(QPoint a_p) { active_point = a_p; }
private:
	bool is_over;
	bool is_bounday;
	int x_max, x_min, y_max, y_min;
	int translation_x, translation_y;
	QPoint active_point;
};

