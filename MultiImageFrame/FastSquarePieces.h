#pragma once
#include "ImageSolver.h"
#include <queue>
#include <ctime>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>

struct FastTraits :public DefaultTraits
{
	typedef Vec2i Point;	// we only need int-value points
};

typedef PolyMesh_ArrayKernelT<FastTraits> PolyMesh;

class FastSquarePieces :
	public ImageSolver
{
public:
	FastSquarePieces();
	~FastSquarePieces();

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
	bool init_token;
	OpenMesh::EPropHandleT<Vec3f> Ux;
	OpenMesh::EPropHandleT<Vec3f> T;
	OpenMesh::EPropHandleT<Vec3f> X;
	OpenMesh::EPropHandleT<Vec3f> Z;
	OpenMesh::EPropHandleT<Vec3f> U;
	OpenMesh::FPropHandleT<Vec3f> Ut;
	OpenMesh::VPropHandleT<Vec3f> V;
	OpenMesh::VPropHandleT<bool> R;

public:
	PolyMesh *mesh;

private:
	void RegistMesh(bool token);
	Vec3f color2vector(QRgb color);
	void UpdateX();
	void UpdateZ();
	float UpdateU();
	void ReconstructSeed(MyMesh::VertexHandle& v_it);
	void ReconstructBack();
	float EnergyFuntion();//This energy function countains no boundary
	void ReconstructDFS();//DFS
	void ReconstructBFS();//BFS
	void ReconstructFFT();
	float BoundaryResidual();
	void RecnstructOneEdge(MyMesh::VertexHandle& v_it, MyMesh::VertexHandle& u_it);
	float UpdateOneStep();
	void UpdateXOnePiece(PolyMesh::FaceIter& f_it);
	void UpdateZOnePiece(PolyMesh::FaceIter& f_it);
	int Location4Edge(PolyMesh::EdgeIter& e_it);
	gMARK information4face(const MyMesh::FaceIter& f_it);

	void BoundaryCheck();
};

