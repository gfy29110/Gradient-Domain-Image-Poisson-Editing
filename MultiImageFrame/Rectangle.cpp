#include "Rectangle.h"



MyRectangle::MyRectangle()
{
	start_x = 0;
	start_y = 0;
	end_x = 0;
	end_y = 0;
	translation_x = 0;
	translation_y = 0;
}


MyRectangle::~MyRectangle()
{
}

void MyRectangle::RectPaint(QPainter & paint)
{
	paint.setPen(rect_pen);
	paint.drawLine(start_x + translation_x, start_y + translation_y, end_x + translation_x, start_y + translation_y);
	paint.drawLine(start_x + translation_x, start_y + translation_y, start_x + translation_x, end_y + translation_y);
	paint.drawLine(end_x + translation_x, end_y + translation_y, end_x + translation_x, start_y + translation_y);
	paint.drawLine(end_x + translation_x, end_y + translation_y, start_x + translation_x, end_y + translation_y);
}

void MyRectangle::set_start(int x, int y)
{
	start_x = x;
	start_y = y;
	end_x = x;
	end_y = y;
}

void MyRectangle::set_end(int x, int y)
{
	if (x > start_x) {
		end_x = x;
	}
	else {
		end_x = start_x;
		start_x = x;
	}
	if (y > start_y) {
		end_y = y;
	}
	else {
		end_y = start_y;
		start_y = y;
	}
}

void MyRectangle::set_translation(int x, int y)
{
	translation_x = x - start_x + 1;
	translation_y = y - start_y + 1;
}

void MyRectangle::reset()
{
	start_x = start_y = 0;
	end_x = end_y = 0;
	translation_x = translation_y = 0;
}

void MyRectangle::paste(MyRectangle * rect)
{
	start_x = rect->get_start_x();
	start_y = rect->get_start_y();
	end_x = rect->get_end_x();
	end_y = rect->get_end_y();
	translation_x = 0;
	translation_y = 0;
}
