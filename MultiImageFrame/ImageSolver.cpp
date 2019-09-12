#include "ImageSolver.h"



ImageSolver::ImageSolver()
{
}


ImageSolver::~ImageSolver()
{
}

void ImageSolver::ImageReigstInnerVetrtice(QImage *I_, MyPolygon mypolygon_)
{
}

void ImageSolver::ImageReigstInnerVetrtice(QImage * I_, MyRectangle myrectangle_)
{
}

void ImageSolver::ImageRegistforTest(QImage *I_, QImage *K_)
{
}

void ImageSolver::ImageRegistBoundary(QImage *J_, int Jx, int Jy)
{
}

void ImageSolver::ParameterSolve()
{
}

void ImageSolver::Reconstruction()
{
}

void ImageSolver::SaveData()
{
}

void ImageSolver::GetEnergy()
{
}

int ImageSolver::change_number(float f)
{
	int k = (int)f;
	if (k < 0) {
		k = 0;
	}
	if (k > 255) {
		k = 255;
	}
	return k;
}

bool ImageSolver::boundary_check_4_test(QRgb pixel)
{
	return((qRed(pixel) == 0 && qGreen(pixel) == 0 && qBlue(pixel) == 0));
}

vector<vector<MARK>> ImageSolver::boundary_expansion_4_test(vector<vector<MARK>> scaninglines)
{
	vector<vector<MARK>> ex_scaninglines;
	vector<MARK> front_line;
	ex_scaninglines.push_back(front_line);

	for (int y = 0; y < scaninglines.size(); y++) {
		if (scaninglines[y].size() != 0) {
			vector<MARK> inner_line;
			if (scaninglines[y][0].x == 0) {
				scaninglines[y][0].is_boundary = true;
			}
			else {
				MARK first_mark;
				first_mark.y = scaninglines[y][0].y;
				first_mark.x = scaninglines[y][0].x - 1;
				first_mark.is_boundary = true;
				inner_line.push_back(first_mark);
			}
			inner_line.insert(inner_line.end(), scaninglines[y].begin(), scaninglines[y].end());
			//this part need make diffenece if we passing the maximum in scaninling line
			MARK last_mark;
			last_mark.y = scaninglines[y][scaninglines[y].size() - 1].y;
			last_mark.x = scaninglines[y][scaninglines[y].size() - 1].x + 1;
			last_mark.is_boundary = true;
			inner_line.push_back(last_mark);
			ex_scaninglines.push_back(inner_line);
		}
	}

	vector<MARK> behind_line;
	ex_scaninglines.push_back(behind_line);

	//expansion
	for (int y = 1; y < ex_scaninglines.size() - 1; y++) {
		for (int x = 1; x < ex_scaninglines[y].size() - 1; x++) {
			MARK mark = ex_scaninglines[y][x];
			//the first and last line don't need cheek
			if (!mark.is_boundary) {
				//upper expansion
				//v_0
				vector<MARK>::iterator it = ex_scaninglines[y - 1].begin();
				if (ex_scaninglines[y - 1].size() != 0) {
					while (it != ex_scaninglines[y - 1].end()) {
						if (it->x >= mark.x - 1) {
							break;
						}
						it++;
					}
					if (it == ex_scaninglines[y - 1].end()) {
						MARK mark_0;
						mark_0.y = mark.y - 1;
						mark_0.x = mark.x - 1;
						mark_0.is_boundary = true;
						it = ex_scaninglines[y - 1].insert(it, mark_0);
					}
					else {
						if (it->x > mark.x - 1) {
							MARK mark_0;
							mark_0.y = mark.y - 1;
							mark_0.x = mark.x - 1;
							mark_0.is_boundary = true;
							it = ex_scaninglines[y - 1].insert(it, mark_0);
						}
					}
				}
				else {
					MARK mark_0;
					mark_0.y = mark.y - 1;
					mark_0.x = mark.x - 1;
					mark_0.is_boundary = true;
					it = ex_scaninglines[y - 1].insert(it, mark_0);
				}
				//v_1
				it++;
				if (it == ex_scaninglines[y - 1].end()) {
					MARK mark_1;
					mark_1.y = mark.y - 1;
					mark_1.x = mark.x;
					mark_1.is_boundary = true;
					it = ex_scaninglines[y - 1].insert(it, mark_1);
				}
				else {
					if (it->x > mark.x) {
						MARK mark_1;
						mark_1.y = mark.y - 1;
						mark_1.x = mark.x;
						mark_1.is_boundary = true;
						it = ex_scaninglines[y - 1].insert(it, mark_1);
					}
				}
				//v_2
				it++;
				if (it == ex_scaninglines[y - 1].end()) {
					MARK mark_2;
					mark_2.y = mark.y - 1;
					mark_2.x = mark.x + 1;
					mark_2.is_boundary = true;
					it = ex_scaninglines[y - 1].insert(it, mark_2);
				}
				else {
					if (it->x > mark.x + 1) {
						MARK mark_2;
						mark_2.y = mark.y - 1;
						mark_2.x = mark.x + 1;
						mark_2.is_boundary = true;
						it = ex_scaninglines[y - 1].insert(it, mark_2);
					}
				}

				//lower expansion
				//v_5
				it = ex_scaninglines[y + 1].begin();
				if (ex_scaninglines[y + 1].size() != 0) {
					while (it != ex_scaninglines[y + 1].end()) {
						if (it->x >= mark.x - 1) {
							break;
						}
						it++;
					}
					if (it == ex_scaninglines[y + 1].end()) {
						MARK mark_5;
						mark_5.y = mark.y + 1;
						mark_5.x = mark.x - 1;
						mark_5.is_boundary = true;
						it = ex_scaninglines[y + 1].insert(it, mark_5);
					}
					else {
						if (it->x > mark.x - 1) {
							MARK mark_5;
							mark_5.y = mark.y + 1;
							mark_5.x = mark.x - 1;
							mark_5.is_boundary = true;
							it = ex_scaninglines[y + 1].insert(it, mark_5);
						}
					}
				}
				else {
					MARK mark_5;
					mark_5.y = mark.y + 1;
					mark_5.x = mark.x - 1;
					mark_5.is_boundary = true;
					it = ex_scaninglines[y + 1].insert(it, mark_5);
				}
				//v_6
				it++;
				if (it == ex_scaninglines[y + 1].end()) {
					MARK mark_6;
					mark_6.y = mark.y + 1;
					mark_6.x = mark.x;
					mark_6.is_boundary = true;
					it = ex_scaninglines[y + 1].insert(it, mark_6);
				}
				else {
					if (it->x > mark.x) {
						MARK mark_6;
						mark_6.y = mark.y + 1;
						mark_6.x = mark.x;
						mark_6.is_boundary = true;
						it = ex_scaninglines[y + 1].insert(it, mark_6);
					}
				}
				//v_7
				it++;
				if (it == ex_scaninglines[y + 1].end()) {
					MARK mark_7;
					mark_7.y = mark.y + 1;
					mark_7.x = mark.x + 1;
					mark_7.is_boundary = true;
					it = ex_scaninglines[y + 1].insert(it, mark_7);
				}
				else {
					if (it->x > mark.x + 1) {
						MARK mark_7;
						mark_7.y = mark.y + 1;
						mark_7.x = mark.x + 1;
						mark_7.is_boundary = true;
						it = ex_scaninglines[y + 1].insert(it, mark_7);
					}
				}
			}
		}
	}

	for (int y = 1; y < ex_scaninglines.size() - 1; y++) {
		vector<MARK>::iterator it = ex_scaninglines[y].begin();
		while (it != ex_scaninglines[y].end()) {
			if (it->is_boundary) {
				it++;
			}
			else {
				//v_3
				vector<MARK>::iterator fit = it;
				fit--;
				if (fit->x + 1 < it->x) {
					MARK mark_3;
					mark_3.x = it->x - 1;
					mark_3.y = it->y;
					mark_3.is_boundary = true;
					it = ex_scaninglines[y].insert(it, mark_3);
					it++;
				}
				//v_4
				vector<MARK>::iterator bit = it;
				bit++;
				if (it->x + 1 < bit->x) {
					MARK mark_4;
					mark_4.x = it->x - 1;
					mark_4.y = it->y;
					mark_4.is_boundary = true;
					it = ex_scaninglines[y].insert(bit, mark_4);
				}
				it++;
			}
		}
	}
	return(ex_scaninglines);
}
