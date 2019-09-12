#include "PoissonSolver.h"



PoissonSolver::PoissonSolver()
{
	for (int i = 0; i < inner_mark.size(); i++) {
		inner_mark[i].clear();
	}
	inner_mark.clear();
	I_x0 = I_y0 = 0;
	J_x0 = J_y0 = 0;

	I = NULL;
	J = NULL;
	regist_I_time = regist_J_time = solve_time = reconstruct_time = 0.0;
}


PoissonSolver::~PoissonSolver()
{
	for (int i = 0; i < inner_mark.size(); i++) {
		inner_mark[i].clear();
	}
	inner_mark.clear();
	I_x0 = I_y0 = 0;
	J_x0 = J_y0 = 0;

	I = NULL;
	J = NULL;
	regist_I_time = regist_J_time = solve_time = reconstruct_time = 0.0;
}

void PoissonSolver::ImageReigstInnerVetrtice(QImage * I_, MyPolygon mypolygon_)
{
	clock_t start_t, end_t;
	start_t = clock();
	I = I_;
	ScanningLine *scanningline = new ScanningLine;
	scanningline->UpdateScanningLine(mypolygon_);
	Regist(scanningline->inner_mark);
	I_x0 = mypolygon_.get_min_x();
	I_y0 = mypolygon_.get_min_y();
	make_matrix();
	end_t = clock();
	regist_I_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
}

void PoissonSolver::ImageReigstInnerVetrtice(QImage * I_, MyRectangle myrectangle_)
{
	clock_t start_t, end_t;
	start_t = clock();
	I = I_;
	ScanningLine *scanningline = new ScanningLine;
	scanningline->UpdateScanningLine(myrectangle_);
	Regist(scanningline->inner_mark);
	I_x0 = myrectangle_.get_start_x();
	I_y0 = myrectangle_.get_start_y();
	make_matrix();
	end_t = clock();
	regist_I_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
}

void PoissonSolver::ImageRegistforTest(QImage *I_, QImage *K_)
{
	clock_t start_t, end_t;
	start_t = clock();
	I = I_;
	inner_mark.clear();
	I_x0 = K_->width();
	I_y0 = K_->height();
	
	for (int j = 0; j < K_->height(); j++) {
		vector<MARK> lines;
		lines.clear();
		for (int i = 0; i < K_->width(); i++) {
			if (qRed(K_->pixel(i, j)) == 0 && qGreen(K_->pixel(i, j)) == 0 && qBlue(K_->pixel(i, j)) == 0) {
				MARK p;
				p.x = i;
				p.y = j;
				if (i < I_x0) {
					I_x0 = i;
				}
				if (j < I_y0) {
					I_y0 = j;
				}
				lines.push_back(p);
			}
		}
		if (lines.size() > 0) {
			inner_mark.push_back(lines);
		}
	}
	make_matrix();
	end_t = clock();
	regist_I_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
}

void PoissonSolver::ImageRegistBoundary(QImage * J_, int Jx, int Jy)
{
	clock_t start_t, end_t;
	start_t = clock();
	J = J_;
	J_x0 = Jx;
	J_y0 = Jy;
	suit_color();
	end_t = clock();
	regist_J_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
}

void PoissonSolver::ParameterSolve()
{
	clock_t start_t, end_t;
	start_t = clock();
	llt.compute(matrix);
	new_red = Eigen::VectorXf(n);
	new_green = Eigen::VectorXf(n);
	new_blue = Eigen::VectorXf(n);;
	new_red = llt.solve(now_red);
	new_green = llt.solve(now_green);
	new_blue = llt.solve(now_blue);
	end_t = clock();
	solve_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
}

void PoissonSolver::Reconstruction()
{
	clock_t start_t, end_t;
	start_t = clock();
	for (int i = 0; i < inner_mark.size(); i++) {
		for (int j = 0; j < inner_mark[i].size(); j++) {
			int v_ID = line_num(i, j);
			QRgb new_color = qRgb(change_number(new_red(v_ID)), change_number(new_green(v_ID)), change_number(new_blue(v_ID)));
			J->setPixel(inner_mark[i][j].x + J_x0 - I_x0, inner_mark[i][j].y + J_y0 - I_y0, new_color);
		}
	}
	end_t = clock();
	reconstruct_time = (double)(end_t - start_t) / CLOCKS_PER_SEC;
}

void PoissonSolver::SaveData()
{
	ofstream data("PoissonSolvingData.txt");

	GetEnergy();
	data << "PoissonSolvingData:" << endl;
	data << "Time of registing interior vertices: " << regist_I_time << "s." << endl;
	data << "Time of registing boundary vertices: " << regist_J_time << "s." << endl;
	data << "Time of solving parameter: " << solve_time << "s." << endl;
	data << "Time of reconstruction: " << reconstruct_time << "s." << endl;

	data << "Final energy is: " << image_energy << endl;
	data.close();
}

float square_float(float a) 
{
	return(a*a);
}

void PoissonSolver::GetEnergy()
{
	image_energy = 0.0;
	for (int i = 0; i < inner_mark.size(); i++) {
		for (int j = 0; j < inner_mark[i].size(); j++) {
			QRgb color;
			MARKID k;
			int v_ID = line_num(i, j);
			QRgb av = I->pixel(inner_mark[i][j].x, inner_mark[i][j].y);
			QRgb au;
			k = search_ID(inner_mark[i][j].x + 1, inner_mark[i][j].y);
			au = I->pixel(inner_mark[i][j].x + 1, inner_mark[i][j].y);
			if (k.i > -1 && k.j > -1) {
				color = J->pixel(inner_mark[i][j].x + 1 + J_x0 - I_x0, inner_mark[i][j].y + J_y0 - I_y0);
				float dr = new_red(v_ID) - (float)qRed(color);
				float dg = new_green(v_ID) - (float)qGreen(color);
				float db = new_blue(v_ID) - (float)qBlue(color);
				float ar = (float)qRed(av) - (float)qRed(au);
				float ag = (float)qGreen(av) - (float)qGreen(au);
				float ab = (float)qBlue(av) - (float)qBlue(au);
				image_energy += square_float(dr - ar);
				image_energy += square_float(dg - ag);
				image_energy += square_float(db - ab);
			}
			else {
				int u_ID = line_num(k.i, k.j);
				float dr = new_red(v_ID) - new_red(u_ID);
				float dg = new_green(v_ID) - new_green(u_ID);
				float db = new_blue(v_ID) - new_blue(u_ID);
				float ar = (float)qRed(av) - (float)qRed(au);
				float ag = (float)qGreen(av) - (float)qGreen(au);
				float ab = (float)qBlue(av) - (float)qBlue(au);
				image_energy += square_float(dr - ar);
				image_energy += square_float(dg - ag);
				image_energy += square_float(db - ab);
			}
			k = search_ID(inner_mark[i][j].x, inner_mark[i][j].y + 1);
			au = I->pixel(inner_mark[i][j].x, inner_mark[i][j].y + 1);
			if (k.i > -1 && k.j > -1) {
				color = J->pixel(inner_mark[i][j].x + J_x0 - I_x0, inner_mark[i][j].y + 1 + J_y0 - I_y0);
				float dr = new_red(v_ID) - (float)qRed(color);
				float dg = new_green(v_ID) - (float)qGreen(color);
				float db = new_blue(v_ID) - (float)qBlue(color);
				float ar = (float)qRed(av) - (float)qRed(au);
				float ag = (float)qGreen(av) - (float)qGreen(au);
				float ab = (float)qBlue(av) - (float)qBlue(au);
				image_energy += square_float(dr - ar);
				image_energy += square_float(dg - ag);
				image_energy += square_float(db - ab);
			}
			else {
				int u_ID = line_num(k.i, k.j);
				float dr = new_red(v_ID) - new_red(u_ID);
				float dg = new_green(v_ID) - new_green(u_ID);
				float db = new_blue(v_ID) - new_blue(u_ID);
				float ar = (float)qRed(av) - (float)qRed(au);
				float ag = (float)qGreen(av) - (float)qGreen(au);
				float ab = (float)qBlue(av) - (float)qBlue(au);
				image_energy += square_float(dr - ar);
				image_energy += square_float(dg - ag);
				image_energy += square_float(db - ab);
			}

			k = search_ID(inner_mark[i][j].x - 1, inner_mark[i][j].y);
			au = I->pixel(inner_mark[i][j].x - 1, inner_mark[i][j].y);
			if (k.i > -1 && k.j > -1) {
				color = J->pixel(inner_mark[i][j].x + J_x0 - I_x0 - 1, inner_mark[i][j].y + J_y0 - I_y0);
				float dr = new_red(v_ID) - (float)qRed(color);
				float dg = new_green(v_ID) - (float)qGreen(color);
				float db = new_blue(v_ID) - (float)qBlue(color);
				float ar = (float)qRed(av) - (float)qRed(au);
				float ag = (float)qGreen(av) - (float)qGreen(au);
				float ab = (float)qBlue(av) - (float)qBlue(au);
				image_energy += square_float(dr - ar);
				image_energy += square_float(dg - ag);
				image_energy += square_float(db - ab);
			}
			k = search_ID(inner_mark[i][j].x, inner_mark[i][j].y - 1);
			au = I->pixel(inner_mark[i][j].x, inner_mark[i][j].y - 1);
			if (k.i > -1 && k.j > -1) {
				color = J->pixel(inner_mark[i][j].x + J_x0 - I_x0, inner_mark[i][j].y + J_y0 - I_y0 - 1);
				float dr = new_red(v_ID) - (float)qRed(color);
				float dg = new_green(v_ID) - (float)qGreen(color);
				float db = new_blue(v_ID) - (float)qBlue(color);
				float ar = (float)qRed(av) - (float)qRed(au);
				float ag = (float)qGreen(av) - (float)qGreen(au);
				float ab = (float)qBlue(av) - (float)qBlue(au);
				image_energy += square_float(dr - ar);
				image_energy += square_float(dg - ag);
				image_energy += square_float(db - ab);
			}
		}
	}
}

void PoissonSolver::Regist(vector<vector<MARK>>* IM)
{
	inner_mark.clear();
	for (auto c_i = IM->begin(); c_i != IM->end(); c_i++) {
		vector<MARK> lines;
		lines.clear();
		for (auto c_ii = c_i->begin(); c_ii != c_i->end(); c_ii++) {
			lines.push_back(*c_ii);
		}
		inner_mark.push_back(lines);
	}
}

void PoissonSolver::make_matrix()
{
	using namespace Eigen;

	n = 0;
	for (int i = 0; i < inner_mark.size(); i++) {
		n += inner_mark[i].size();
	}
	matrix = SparseMatrix <float>(n, n);
	std::vector <Eigen::Triplet<float>> triplet;
	orign_red = Eigen::VectorXf::Constant(n, 0);
	orign_green = Eigen::VectorXf::Constant(n, 0);
	orign_blue = Eigen::VectorXf::Constant(n, 0);
	for (int i = 0; i < inner_mark.size(); i++) {
		for (int j = 0; j < inner_mark[i].size(); j++) {
			MARKID k;
			int v_ID = line_num(i, j);
			int u_ID;

			Eigen::Triplet<float> temp(v_ID, v_ID, -4);
			triplet.push_back(temp);
			k = search_ID(inner_mark[i][j].x + 1, inner_mark[i][j].y);
			u_ID = line_num(k.i, k.j);
			if (u_ID != -1) {
				Eigen::Triplet<float> temp(v_ID, u_ID, 1);
				triplet.push_back(temp);
			}
			k = search_ID(inner_mark[i][j].x - 1, inner_mark[i][j].y);
			u_ID = line_num(k.i, k.j);
			if (u_ID != -1) {
				Eigen::Triplet<float> temp(v_ID, u_ID, 1);
				triplet.push_back(temp);

			}
			k = search_ID(inner_mark[i][j].x, inner_mark[i][j].y + 1);
			u_ID = line_num(k.i, k.j);
			if (u_ID != -1) {
				Eigen::Triplet<float> temp(v_ID, u_ID, 1);
				triplet.push_back(temp);
			}
			k = search_ID(inner_mark[i][j].x, inner_mark[i][j].y - 1);
			u_ID = line_num(k.i, k.j);
			if (u_ID != -1) {
				Eigen::Triplet<float> temp(v_ID, u_ID, 1);
				triplet.push_back(temp);
			}
			QRgb color = I->pixel(inner_mark[i][j].x, inner_mark[i][j].y);
			orign_red(v_ID) += -4 * qRed(color);
			orign_green(v_ID) += -4 * qGreen(color);
			orign_blue(v_ID) += -4 * qBlue(color);
			color = I->pixel(inner_mark[i][j].x + 1, inner_mark[i][j].y);
			orign_red(v_ID) += qRed(color);
			orign_green(v_ID) += qGreen(color);
			orign_blue(v_ID) += qBlue(color);
			color = I->pixel(inner_mark[i][j].x - 1, inner_mark[i][j].y);
			orign_red(v_ID) += qRed(color);
			orign_green(v_ID) += qGreen(color);
			orign_blue(v_ID) += qBlue(color);
			color = I->pixel(inner_mark[i][j].x, inner_mark[i][j].y + 1);
			orign_red(v_ID) += qRed(color);
			orign_green(v_ID) += qGreen(color);
			orign_blue(v_ID) += qBlue(color);
			color = I->pixel(inner_mark[i][j].x, inner_mark[i][j].y - 1);
			orign_red(v_ID) += qRed(color);
			orign_green(v_ID) += qGreen(color);
			orign_blue(v_ID) += qBlue(color);
		}
	}
	matrix.setFromTriplets(triplet.begin(), triplet.end());
}

void PoissonSolver::suit_color()
{
	if (matrix.size() != 0) {
		now_red = Eigen::VectorXf::Constant(n, 0);
		now_green = Eigen::VectorXf::Constant(n, 0);
		now_blue = Eigen::VectorXf::Constant(n, 0);
		for (int i = 0; i < inner_mark.size(); i++) {
			for (int j = 0; j < inner_mark[i].size(); j++) {
				int v_ID = line_num(i, j);
				now_red(v_ID) = orign_red(v_ID);
				now_green(v_ID) = orign_green(v_ID);
				now_blue(v_ID) = orign_blue(v_ID);
				QRgb color;
				MARKID k;
				k = search_ID(inner_mark[i][j].x + 1, inner_mark[i][j].y);
				if (k.i == -1 || k.j == -1) {
					color = J->pixel(inner_mark[i][j].x + 1 + J_x0 - I_x0, inner_mark[i][j].y + J_y0 - I_y0);//this part have some problem
					now_red(v_ID) -= qRed(color);
					now_green(v_ID) -= qGreen(color);
					now_blue(v_ID) -= qBlue(color);
				}
				k = search_ID(inner_mark[i][j].x - 1, inner_mark[i][j].y);
				if (k.i == -1 || k.j == -1) {
					color = J->pixel(inner_mark[i][j].x - 1 + J_x0 - I_x0, inner_mark[i][j].y + J_y0 - I_y0);
					now_red(v_ID) -= qRed(color);
					now_green(v_ID) -= qGreen(color);
					now_blue(v_ID) -= qBlue(color);
				}
				k = search_ID(inner_mark[i][j].x, inner_mark[i][j].y + 1);
				if (k.i == -1 || k.j == -1) {
					color = J->pixel(inner_mark[i][j].x + J_x0 - I_x0, inner_mark[i][j].y + 1 + J_y0 - I_y0);
					now_red(v_ID) -= qRed(color);
					now_green(v_ID) -= qGreen(color);
					now_blue(v_ID) -= qBlue(color);
				}
				k = search_ID(inner_mark[i][j].x, inner_mark[i][j].y - 1);
				if (k.i == -1 || k.j == -1) {
					color = J->pixel(inner_mark[i][j].x + J_x0 - I_x0, inner_mark[i][j].y - 1 + J_y0 - I_y0);
					now_red(v_ID) -= qRed(color);
					now_green(v_ID) -= qGreen(color);
					now_blue(v_ID) -= qBlue(color);
				}
			}
		}
	}
}

MARKID PoissonSolver::search_ID(int x, int y)
{
	MARKID token;
	token.i = -1;
	token.j = -1;

	for (int i = 0; i < inner_mark.size(); i++) {
		for (int j = 0; j < inner_mark[i].size(); j++) {
			if (x == inner_mark[i][j].x&&y == inner_mark[i][j].y) {
				token.i = i;
				token.j = j;
			}
		}
	}
	return token;
}

int PoissonSolver::line_num(int i, int j)
{
	if (i >= inner_mark.size() || i < 0) {
		return -1;
	}
	
	if (j >= inner_mark[i].size() || j < 0) {
		return -1;
	}
	else {
		int k = 0;
		for (int m = 0; m < i; m++) {
			k += inner_mark[m].size();
		}
		k += j;
		return k;
	}
}
