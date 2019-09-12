#pragma once
#include "ImageSolver.h"
#include <queue>
#include <ctime>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>


class SquarePieces2 :
	public ImageSolver
{
public:
	SquarePieces2();
	~SquarePieces2();

	void ImageReigstInnerVetrtice(QImage *I_, MyPolygon mypolygon_);
	void ImageReigstInnerVetrtice(QImage *I_, MyRectangle myrectangle_);
	void ImageRegistforTest(QImage *I_, QImage *K_);
	void ImageRegistBoundary(QImage *J_, int Jx, int Jy);
	void ParameterSolve();
	void Reconstruction();
	void SaveData();
	void GetEnergy();

private:
	float rho;
	QRgb save_rgb;
	MyMesh::VertexHandle first_handle;
	ReconstructType RT;
	int iter;
	OpenMesh::HPropHandleT<Vec3f> Ux;

public:
	MyMesh *mesh;

private:
	void RegistSquarePiece(MyMesh::FaceHandle& f_handle, MARK p00, MARK p01, MARK p11, MARK p10, bool token);
	void RegistSquarePiece_J(MyMesh::FaceIter& f_it, bool token);
	void UpdateXOnePiece(MyMesh::FaceIter& f_it);
	void UpdateZOnePiece(MyMesh::FaceIter& f_it);
	float UpdateUOnePiece(MyMesh::FaceIter& f_it);
	void PassingXOnePiece(MyMesh::FaceIter& f_it);
	void PassingZOnePiece(MyMesh::FaceIter& f_it);
	void UpdateX();
	void UpdateZ();
	float UpdateU();
	float UpdateOneStep();
	void UpdateParameter();
	void ReconstructSeed(MyMesh::VertexHandle& v_it);
	void ReconstructBack();
	float EnergyFuntion();//This energy function countains boundary
	void ReconstructDFS();//DFS
	void ReconstructBFS();//BFS
	void ReconstructFFT();
	float BoundaryResidual();

private:
	void RegistMesh(bool token);
	Vec3f CalculateDUt(MyMesh::HalfedgeHandle hf);
	Vec3f color2vector(QRgb color);
	void CalculateCurl(Vec3f & curl, MyMesh::HalfedgeHandle hf);
	float EnergyofOnePiece(MyMesh::FaceIter f_it);
	void RecnstructOneEdge(MyMesh::VertexHandle& v_it, MyMesh::VertexHandle& u_it);
	void AdjustofHalfedge();
	void AdjustBoundaryData();
	gMARK information4face(const MyMesh::FaceIter& f_it);

	//This part for debug
	void BoundaryCheck();
};

