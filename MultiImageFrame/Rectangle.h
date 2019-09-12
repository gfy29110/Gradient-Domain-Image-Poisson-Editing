#pragma once
#include <QPainter>
#include <qpen.h>
#include <math.h>

class MyRectangle
{
	const QPen rect_pen = QPen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);

public:
	MyRectangle();
	~MyRectangle();
	void RectPaint(QPainter &paint);

	void set_start(int x, int y);
	void set_end(int x, int  y);
	void set_translation(int x, int y);

	int get_start_x() { return start_x; }
	int get_start_y() { return start_y; }
	int get_end_x() { return end_x; }
	int get_end_y() { return end_y; }

	int get_width() { return abs(end_x - start_x + 1); }
	int get_height() { return abs(end_y - start_y + 1); }

	void reset();

	void paste(MyRectangle *rect);

private:
	int start_x, start_y;
	int end_x, end_y;
	int translation_x, translation_y;
};

