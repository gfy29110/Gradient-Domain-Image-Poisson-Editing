#pragma once
#include "ImageSolver.h"
#include <Eigen\Dense>
#include <Eigen\Sparse>
#include <math.h>
#include <fstream>
#include <iostream>

typedef struct MID {
	int i;
	int j;
}MARKID;

class PoissonSolver :
	public ImageSolver
{
public:
	PoissonSolver();
	~PoissonSolver();
	void ImageReigstInnerVetrtice(QImage *I_, MyPolygon mypolygon_);
	void ImageReigstInnerVetrtice(QImage *I_, MyRectangle myrectangle_);
	void ImageRegistforTest(QImage *I_, QImage *K_);
	void ImageRegistBoundary(QImage *J_, int Jx, int Jy);
	void ParameterSolve();
	void Reconstruction();
	void SaveData();
	void GetEnergy();

private:
	Eigen::SparseMatrix<float> matrix;
	Eigen::VectorXf orign_red, orign_green, orign_blue, now_red, now_green, now_blue, new_red, new_green, new_blue;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> llt;
	int n;
	//ScanningLine *scanningline;

	void Regist(vector<vector<MARK>> *IM);
	void make_matrix();
	void suit_color();
	MARKID search_ID(int x, int y);
	int line_num(int i, int j);
};

