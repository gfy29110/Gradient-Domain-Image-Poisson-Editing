#pragma once
#ifndef _MY_IMAGESOLVER_
#define _MY_IMAGESOLVER_



#include "ScanningLine.h"
#include <ctime>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

using namespace OpenMesh;

enum LocationType {
	L_Location,
	U_Location,
	R_Location,
	D_Location
};

typedef struct gradientMARK
{
	int x;
	int y;
	bool boundary;
	Vec3f g_x;
	Vec3f g_y;

}gMARK;

enum ReconstructType {
	FFT,
	BFS,
	DFS
};

struct MyTraits :public DefaultTraits
{
	typedef Vec2i Point;	// we only need int-value points
	HalfedgeTraits{
	private:
		LocationType h_location;
		Vec3f X;
		Vec3f Z;
		Vec3f U;
		Vec3f T;
		bool weight;
	public:
		void setLocation(LocationType L) { h_location = L; }
		void setX(const Vec3f& value) { X = value; }
		void setZ(const Vec3f& value) { Z = value; }
		void setU(const Vec3f& value) { U = value; }
		void setT(const Vec3f& value) { T = value; }
		void setW(bool w) { weight = w; }
		LocationType getLocation() { return h_location; }
		const Vec3f& getX() const { return X; }
		const Vec3f& getZ() const { return Z; }
		const Vec3f& getU() const { return U; }
		const Vec3f& getT() const { return T; }
		bool getW() { return weight; }
		void init() {
			setLocation(L_Location);
			setX(Vec3f(0.0, 0.0, 0.0));
			setZ(Vec3f(0.0, 0.0, 0.0));
			setU(Vec3f(0.0, 0.0, 0.0));
			setT(Vec3f(0.0, 0.0, 0.0));
			setW(false);
		}
	};
	VertexTraits{
	private:
		Vec3f V;
		bool is_reconstructed;
		bool is_boundary;
	public:
		void setV(const Vec3f& value) { V = value; }
		void setR(bool token) { is_reconstructed = token; }
		void setB(bool token) { is_boundary = token; }
		const Vec3f& getV() const { return V; }
		bool getR() { return is_reconstructed; }
		bool getB() { return is_boundary; }
		void init() {
			setV(Vec3f(0.0, 0.0, 0.0));
			setR(false);
			setB(false);
		}
	};
	FaceTraits{
	private:
		Vec3f U_t;
	public:
		void setU(const Vec3f& value) { U_t = value; }
		const Vec3f& getU() const { return U_t; }
		void init() {
			setU(Vec3f(0.0, 0.0, 0.0));
		}
	};
};

typedef PolyMesh_ArrayKernelT<MyTraits> MyMesh;

class ImageSolver
{
public:
	ImageSolver();
	virtual ~ImageSolver();

	virtual void ImageReigstInnerVetrtice(QImage *I_, MyPolygon mypolygon_);
	virtual void ImageReigstInnerVetrtice(QImage *I_, MyRectangle myrectangle_);
	virtual void ImageRegistforTest(QImage *I_, QImage *K_);
	virtual void ImageRegistBoundary(QImage *J_, int Jx, int Jy);
	virtual void ParameterSolve();
	virtual void Reconstruction();
	virtual void SaveData();
	virtual void GetEnergy();

protected:
	vector<vector<MARK>> inner_mark;
	int I_x0, I_y0;
	int J_x0, J_y0;

	QImage* I;
	QImage* J;

	double regist_I_time, regist_J_time, solve_time, reconstruct_time;
	float image_energy;

	int change_number(float f);

	bool boundary_check_4_test(QRgb pixel);
	vector<vector<MARK>> boundary_expansion_4_test(vector<vector<MARK>> scaninglines);
};





#endif // !_MY_IMAGESOLVER_