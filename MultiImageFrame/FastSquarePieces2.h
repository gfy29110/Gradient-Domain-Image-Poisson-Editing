#pragma once
#include "ImageSolver.h"
#include <queue>
#include <ctime>
#include <fstream>
#include <omp.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>


struct FastTraits :public DefaultTraits
{
	typedef Vec2i Point;	// we only need int-value points
};

typedef PolyMesh_ArrayKernelT<FastTraits> PolyMesh;

class FastSquarePieces2 :
	public ImageSolver
{
public:
	FastSquarePieces2();
	~FastSquarePieces2();

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
	PolyMesh::VertexHandle first_handle;
	ReconstructType RT;
	int iter;
	OpenMesh::HPropHandleT<Vec3f> Ux;
	OpenMesh::HPropHandleT<Vec3f> T;
	OpenMesh::HPropHandleT<Vec3f> X;
	OpenMesh::HPropHandleT<Vec3f> Z;
	OpenMesh::HPropHandleT<Vec3f> U;
	OpenMesh::HPropHandleT<bool> W;
	OpenMesh::FPropHandleT<Vec3f> Ut;
	OpenMesh::VPropHandleT<Vec3f> V;
	OpenMesh::VPropHandleT<bool> R;
	OpenMesh::HPropHandleT<LocationType> L;

	int num_v, num_f, num_h, num_e;

public:
	PolyMesh *mesh;

private:
	void RegistSquarePiece(PolyMesh::FaceHandle& f_handle, MARK p00, MARK p01, MARK p11, MARK p10, bool token);
	void RegistSquarePiece_J(PolyMesh::FaceHandle f_it, bool token);
	void UpdateXOnePiece(PolyMesh::FaceHandle f_it);
	void UpdateZOnePiece(PolyMesh::FaceHandle f_it);
	float UpdateUOnePiece(PolyMesh::FaceHandle f_it);
	void PassingXOnePiece(PolyMesh::FaceHandle f_it);
	void PassingZOnePiece(PolyMesh::FaceHandle f_it);
	void UpdateX();
	void UpdateZ();
	float UpdateU();
	float UpdateOneStep();
	void UpdateParameter();
	void ReconstructSeed(PolyMesh::VertexHandle& v_it);
	void ReconstructBack();
	float EnergyFuntion();//This energy function countains boundary
	void ReconstructDFS();//DFS
	void ReconstructBFS();//BFS
	void ReconstructFFT();
	float BoundaryResidual();

private:
	void RegistMesh(bool token);
	Vec3f CalculateDUt(PolyMesh::HalfedgeHandle hf);
	Vec3f color2vector(QRgb color);
	void CalculateCurl(Vec3f & curl, PolyMesh::HalfedgeHandle hf);
	float EnergyofOnePiece(PolyMesh::FaceIter f_it);
	void RecnstructOneEdge(PolyMesh::VertexHandle& v_it, PolyMesh::VertexHandle& u_it);
	void AdjustofHalfedge();
	void AdjustBoundaryData();
	gMARK information4face(const PolyMesh::FaceIter& f_it);

	//This part for debug
	void BoundaryCheck();
};

