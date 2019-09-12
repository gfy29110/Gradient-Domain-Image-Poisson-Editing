#include "FastSquarePieces2.h"



FastSquarePieces2::FastSquarePieces2()
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
	RT = BFS;
	save_rgb = QRgb(0);
	rho = 100.0;

	mesh = new PolyMesh;
	mesh->add_property(T);
	mesh->add_property(X);
	mesh->add_property(Z);
	mesh->add_property(U);
	mesh->add_property(Ux);
	mesh->add_property(W);
	mesh->add_property(Ut);
	mesh->add_property(V);
	mesh->add_property(R);
	mesh->add_property(L);

	int coreNum = omp_get_num_procs();
	cout << "The number of core is: " << coreNum << endl;
}


FastSquarePieces2::~FastSquarePieces2()
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
	RT = BFS;
	save_rgb = QRgb(0);
	rho = 100.0;

	mesh = new PolyMesh;
	mesh->add_property(T);
	mesh->add_property(X);
	mesh->add_property(Z);
	mesh->add_property(U);
	mesh->add_property(Ux);
	mesh->add_property(W);
	mesh->add_property(Ut);
	mesh->add_property(V);
	mesh->add_property(R);
	mesh->add_property(L);
}

void FastSquarePieces2::ImageReigstInnerVetrtice(QImage * I_, MyPolygon mypolygon_)
{
	clock_t start_t, end_t;
	start_t = clock();
	I = I_;
	ScanningLine* scanningline = new ScanningLine;
	scanningline->UpdateScanningLine(mypolygon_);
	I_x0 = mypolygon_.get_min_x();
	I_y0 = mypolygon_.get_min_y();
	inner_mark = *(scanningline->inner_mark);
	RegistMesh(true);
	end_t = clock();
	regist_I_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}

void FastSquarePieces2::ImageReigstInnerVetrtice(QImage * I_, MyRectangle myrectangle_)
{
	clock_t start_t, end_t;
	start_t = clock();
	I = I_;
	ScanningLine* scanningline = new ScanningLine;
	scanningline->UpdateScanningLine(myrectangle_);
	I_x0 = myrectangle_.get_start_x();
	I_y0 = myrectangle_.get_start_y();
	inner_mark = *(scanningline->inner_mark);
	RegistMesh(true);
	end_t = clock();
	regist_I_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}

void FastSquarePieces2::ImageRegistforTest(QImage *I_, QImage *K_)
{
	clock_t start_t, end_t;
	start_t = clock();
	I = I_;
	//vector<vector<MARK>> inner_mark;
	inner_mark.clear();
	I_x0 = K_->width();
	I_y0 = K_->height();
	for (int j = 0; j < K_->height(); j++) {
		vector<MARK> scanning_line;
		for (int i = 0; i < K_->width(); i++) {
			if (boundary_check_4_test(K_->pixel(i, j))) {
				MARK p;
				p.x = i;
				p.y = j;
				p.is_boundary = false;
				if (i < I_x0) {
					I_x0 = i;
				}
				if (j < I_y0) {
					I_y0 = j;
				}
				scanning_line.push_back(p);
			}
		}
		if (scanning_line.size() != 0) {
			inner_mark.push_back(scanning_line);
		}
	}

	inner_mark = boundary_expansion_4_test(inner_mark);

	I_x0--;
	I_y0--;



	for (int i = 0; i < inner_mark.size(); i++) {
		auto it = inner_mark[i].begin();
		auto fit = it + 1;
		while (fit != inner_mark[i].end()) {
			if (it->x == fit->x) {
				it = inner_mark[i].erase(it);
				fit = it + 1;
			}
			else {
				fit++;
				it++;
			}
		}
	}

	RegistMesh(true);
	end_t = clock();
	regist_I_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}

void FastSquarePieces2::ImageRegistBoundary(QImage * J_, int Jx, int Jy)
{
	clock_t start_t, end_t;
	start_t = clock();
	J = J_;
	J_x0 = Jx;
	J_y0 = Jy;
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_f; i++) {
		auto f_it = mesh->face_handle(i);
		RegistSquarePiece_J(f_it, true);
	}
	save_rgb = J->pixel(mesh->point(first_handle)[0] + J_x0 - I_x0, mesh->point(first_handle)[1] + J_y0 - I_y0);
#ifndef NDEBUG
	cout << "Size of Image I: " << I->width() << "*" << I->height() << endl;
	cout << "Size of Image J: " << J->width() << "*" << J->height() << endl;
#endif // !NDEBUG
	AdjustofHalfedge();
	//cout << "Number of Square is:" << square_num << endl;
#ifndef NDEBUG
	BoundaryCheck();
#endif // !NDEBUG
	end_t = clock();
	regist_J_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}

void FastSquarePieces2::ParameterSolve()
{
	UpdateParameter();
}

void FastSquarePieces2::Reconstruction()
{
	clock_t start_t, end_t;
	start_t = clock();
	switch (RT)
	{
	case FFT:
		ReconstructFFT();
		break;
	case BFS:
		ReconstructBFS();
		break;
	case DFS:
		ReconstructDFS();
		break;
	default:
		break;
	}
	end_t = clock();
	reconstruct_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}

void FastSquarePieces2::SaveData()
{
	ofstream data("FastSquarePieces2Data.txt");

	GetEnergy();
	data << "FastSquarePieces2Data:" << endl;
	data << "Iteration: " << iter << " step." << endl;
	data << "Reconstruction Type: ";
	switch (RT) {
	case FFT:
		data << "FFT";
		break;
	case DFS:
		data << "DFS";
		break;
	case BFS:
		data << "BFS";
		break;
	default:
		break;
	}
	data << endl;
	data << "Time of registing interior vertices: " << regist_I_time << "s." << endl;
	data << "Time of registing boundary vertices: " << regist_J_time << "s." << endl;
	data << "Time of solving parameter: " << solve_time << "s." << endl;
	data << "Time of reconstruction: " << reconstruct_time << "s." << endl;

	data << "Final energy is: " << image_energy << endl;

	data << "Residual of Boundary is: " << BoundaryResidual() << endl;

	data.close();
}

void FastSquarePieces2::GetEnergy()
{
	image_energy = 0.0;
	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
		for (auto voh_it = mesh->voh_begin(*v_it); voh_it != mesh->voh_end(*v_it); voh_it++) {
			if ((mesh->property(L, *voh_it) == L_Location || mesh->property(L, *voh_it)== D_Location) && !mesh->is_boundary(mesh->edge_handle(*voh_it))) {
				image_energy += (mesh->property(T, *voh_it) - mesh->property(X, *voh_it)).sqrnorm();
			}
		}
	}
}

//This part has replaced U and D
void FastSquarePieces2::RegistSquarePiece(PolyMesh::FaceHandle & f_handle, MARK p00, MARK p01, MARK p11, MARK p10, bool token)
{
	Vec3f P, Q;
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_handle);
	//01->11
	mesh->property(L, *fh_it) = U_Location;
	P = color2vector(I->pixel(p01.x, p01.y));
	Q = color2vector(I->pixel(p11.x, p11.y));
	mesh->property(T, *fh_it) = Q - P;
	if (token) {
		mesh->property(X, *fh_it) = Q - P;
		mesh->property(Z, *fh_it) = Q - P;
		mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
	}
	//00->01
	++fh_it;
	mesh->property(L, *fh_it) = L_Location;
	P = color2vector(I->pixel(p00.x, p00.y));
	Q = color2vector(I->pixel(p01.x, p01.y));
	mesh->property(T, *fh_it) = Q - P;
	if (token) {
		mesh->property(X, *fh_it) = Q - P;
		mesh->property(Z, *fh_it) = Q - P;
		mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
	}
	//10->00
	++fh_it;
	mesh->property(L, *fh_it) = D_Location;
	P = color2vector(I->pixel(p00.x, p00.y));
	Q = color2vector(I->pixel(p10.x, p10.y));
	mesh->property(T, *fh_it) = Q - P;
	if (token) {
		mesh->property(X, *fh_it) = Q - P;
		mesh->property(Z, *fh_it) = Q - P;
		mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
	}
	//11->10
	++fh_it;
	mesh->property(L, *fh_it) = R_Location;
	P = color2vector(I->pixel(p10.x, p10.y));
	Q = color2vector(I->pixel(p11.x, p11.y));
	mesh->property(T, *fh_it) = Q - P;
	if (token) {
		mesh->property(X, *fh_it) = Q - P;
		mesh->property(Z, *fh_it) = Q - P;
		mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
	}
}

void FastSquarePieces2::RegistSquarePiece_J(PolyMesh::FaceHandle f_it, bool token)
{
	Vec3f P, Q, U_t;
	U_t = Vec3f(0, 0, 0);
	PolyMesh::FaceVertexIter fv_it = mesh->fv_begin(f_it);
	//01->00->10->11
	PolyMesh::Point p00, p01, p11, p10;
	p01 = mesh->point(*fv_it);
	fv_it++;
	p00 = mesh->point(*fv_it);
	fv_it++;
	p10 = mesh->point(*fv_it);
	fv_it++;
	p11 = mesh->point(*fv_it);
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	//01->11
	if (mesh->is_boundary(mesh->edge_handle(*fh_it))) {
		P = color2vector(J->pixel(J_x0 - I_x0 + p01[0] + 1, J_y0 - I_y0 + p01[1] + 1));
		Q = color2vector(J->pixel(J_x0 - I_x0 + p11[0] + 1, J_y0 - I_y0 + p11[1] + 1));
		mesh->property(T, *fh_it) = Q - P;
		mesh->property(W, *fh_it) = true;
		if (token) {
			mesh->property(X, *fh_it) = Q - P;
			mesh->property(Z, *fh_it) = Q - P;
			mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
		}
		U_t = U_t - mesh->property(X, *fh_it);
	}
	else {
		U_t = U_t - mesh->property(X, *fh_it);
	}
	++fh_it;
	//00->01
	if (mesh->is_boundary(mesh->edge_handle(*fh_it))) {
		P = color2vector(J->pixel(J_x0 - I_x0 + p00[0] + 1, J_y0 - I_y0 + p00[1] + 1));
		Q = color2vector(J->pixel(J_x0 - I_x0 + p01[0] + 1, J_y0 - I_y0 + p01[1] + 1));
		mesh->property(T, *fh_it) = Q - P;
		mesh->property(W, *fh_it) = true;
		if (token) {
			mesh->property(X, *fh_it) = Q - P;
			mesh->property(Z, *fh_it) = Q - P;
			mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
		}
		U_t = U_t - mesh->property(X, *fh_it);
	}
	else {
		U_t = U_t - mesh->property(X, *fh_it);
	}
	++fh_it;
	//10->00
	if (mesh->is_boundary(mesh->edge_handle(*fh_it))) {
		P = color2vector(J->pixel(J_x0 - I_x0 + p00[0] + 1, J_y0 - I_y0 + p00[1] + 1));
		Q = color2vector(J->pixel(J_x0 - I_x0 + p10[0] + 1, J_y0 - I_y0 + p10[1] + 1));
		mesh->property(T, *fh_it) = Q - P;
		mesh->property(W, *fh_it) = true;
		if (token) {
			mesh->property(X, *fh_it) = Q - P;
			mesh->property(Z, *fh_it) = Q - P;
			mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
		}
		U_t = U_t + mesh->property(X, *fh_it);
	}
	else {
		U_t = U_t + mesh->property(X, *fh_it);
	}
	++fh_it;
	//11->10
	if (mesh->is_boundary(mesh->edge_handle(*fh_it))) {
		P = color2vector(J->pixel(J_x0 - I_x0 + p10[0] + 1, J_y0 - I_y0 + p10[1] + 1));
		Q = color2vector(J->pixel(J_x0 - I_x0 + p11[0] + 1, J_y0 - I_y0 + p11[1] + 1));
		mesh->property(T, *fh_it) = Q - P;
		mesh->property(W, *fh_it) = true;
		if (token) {
			mesh->property(X, *fh_it) = Q - P;
			mesh->property(Z, *fh_it) = Q - P;
			mesh->property(U, *fh_it) = Vec3f(0.0, 0.0, 0.0);
		}
		U_t = U_t + mesh->property(X, *fh_it);
	}
	else {
		U_t = U_t + mesh->property(X, *fh_it);
	}

	mesh->property(Ut, f_it) = U_t;
}

void FastSquarePieces2::UpdateXOnePiece(PolyMesh::FaceHandle f_it)
{
	float x_weight, y_weight;
	OpenMesh::Vec3f WZt = -mesh->property(Ut, f_it);
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	PolyMesh::FaceHalfedgeIter l_fh_it;
	PolyMesh::FaceHalfedgeIter d_fh_it;
	OpenMesh::Vec3f T_l, WZl, T_d, WZd;
	while (fh_it != mesh->fh_end(f_it)) {
		if (mesh->property(L, *fh_it) == L_Location) {
			if (mesh->property(W, *fh_it)) {
				T_l = mesh->property(T, *fh_it) - mesh->property(Ux, *fh_it);
			}
			else {
				T_l = mesh->property(T, *fh_it);
			}
			WZl = mesh->property(Z, *fh_it) - mesh->property(U, *fh_it);
			y_weight = mesh->property(W, *fh_it) ? rho : 1.0;
			l_fh_it = fh_it;
		}
		if (mesh->property(L, *fh_it) == D_Location) {
			if (mesh->property(W, *fh_it)) {
				T_d = mesh->property(T, *fh_it) - mesh->property(Ux, *fh_it);
			}
			else {
				T_d = mesh->property(T, *fh_it);
			}
			WZd = mesh->property(Z, *fh_it) - mesh->property(U, *fh_it);
			x_weight = mesh->property(W, *fh_it) ? rho : 1.0;
			d_fh_it = fh_it;
		}
		if (mesh->property(L, *fh_it) == U_Location) {
			WZt += mesh->property(Z, *fh_it);
		}
		if (mesh->property(L, *fh_it) == R_Location) {
			WZt -= mesh->property(Z, *fh_it);
		}
		fh_it++;
	}
	OpenMesh::Vec3f X_l, X_d;
	float denominator = x_weight * y_weight + 2 * (x_weight + y_weight)*rho + 3 * rho*rho;
	X_d = ((2 * rho + y_weight)* (x_weight*T_d + rho * WZd) + (y_weight*rho + rho * rho)*WZt + rho * (y_weight*T_l + rho * WZl)) / denominator;
	X_l = (rho * (x_weight*T_d + rho * WZd) - (x_weight*rho + rho * rho)*WZt + (2 * rho + x_weight)* (y_weight*T_l + rho * WZl)) / denominator;

	mesh->property(X, l_fh_it) = X_l;
	mesh->property(X, d_fh_it) = X_d;
}

void FastSquarePieces2::UpdateZOnePiece(PolyMesh::FaceHandle f_it)
{
	float x_weight, y_weight;
	OpenMesh::Vec3f WXt = mesh->property(Ut, f_it);
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	PolyMesh::FaceHalfedgeIter u_fh_it;
	PolyMesh::FaceHalfedgeIter r_fh_it;
	OpenMesh::Vec3f WXr, WXu;

	while (fh_it != mesh->fh_end(f_it)) {
		if (mesh->property(L, *fh_it) == L_Location) {
			WXt -= mesh->property(X, *fh_it);
		}
		if (mesh->property(L, *fh_it) == D_Location) {
			WXt += mesh->property(X, *fh_it);
		}
		if (mesh->property(L, *fh_it) == U_Location) {
			WXu = mesh->property(X, *fh_it) + mesh->property(U, *fh_it);
			u_fh_it = fh_it;
		}
		if (mesh->property(L, *fh_it) == R_Location) {
			WXr = mesh->property(X, *fh_it) + mesh->property(U, *fh_it);
			r_fh_it = fh_it;
		}
		fh_it++;
	}
	OpenMesh::Vec3f Z_r, Z_u;
	Z_u = (2 * WXu + WXr + WXt) / 3;
	Z_r = (WXu + 2 * WXr - WXt) / 3;

	mesh->property(Z, r_fh_it) = Z_r;
	mesh->property(Z, u_fh_it) = Z_u;
}

float FastSquarePieces2::UpdateUOnePiece(PolyMesh::FaceHandle f_it)
{
	float residual = 0.0;
	Vec3f U_t = Vec3f(0.0, 0.0, 0.0);
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	while (fh_it != mesh->fh_end(f_it)) {
		mesh->property(U, *fh_it) = mesh->property(U, *fh_it) + mesh->property(X, *fh_it) - mesh->property(Z, *fh_it);
		float part_residual = (mesh->property(X, *fh_it) - mesh->property(Z, *fh_it)).max_abs();
		if (residual < part_residual) {
			residual = part_residual;
		}
		U_t += CalculateDUt(*fh_it);
		++fh_it;
	}

	if (residual < U_t.max_abs()) {
		residual = U_t.max_abs();
	}
	U_t += mesh->property(Ut, f_it);
	mesh->property(Ut, f_it) = U_t;

	return(residual);
}

void FastSquarePieces2::PassingXOnePiece(PolyMesh::FaceHandle f_it)
{
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	while (fh_it != mesh->fh_end(f_it)) {
		if (mesh->property(L, *fh_it) == R_Location || mesh->property(L, *fh_it) == U_Location) {
			PolyMesh::HalfedgeHandle opposite_heh = mesh->opposite_halfedge_handle(*fh_it);
			if (!mesh->is_boundary(opposite_heh)) {
				mesh->property(X, *fh_it) = mesh->property(X, opposite_heh);
			}
			else {
				OpenMesh::Vec3f nX = ((mesh->property(Z, *fh_it) - mesh->property(U, *fh_it)) + (mesh->property(T, *fh_it) - mesh->property(Ux, *fh_it))) / 2;
				mesh->property(X, *fh_it) = nX;
			}
		}
		++fh_it;
	}
}

void FastSquarePieces2::PassingZOnePiece(PolyMesh::FaceHandle f_it)
{
	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	while (fh_it != mesh->fh_end(f_it)) {
		if (mesh->property(L, *fh_it) == L_Location || mesh->property(L, *fh_it) == D_Location) {
			PolyMesh::HalfedgeHandle opposite_heh = mesh->opposite_halfedge_handle(*fh_it);
			if (!mesh->is_boundary(opposite_heh)) {
				mesh->property(Z, *fh_it) = mesh->property(Z, opposite_heh);
			}
			else {
				mesh->property(Z, *fh_it) = mesh->property(X, *fh_it) + mesh->property(U, *fh_it);
			}
		}
		++fh_it;
	}
}

void FastSquarePieces2::UpdateX()
{
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_f; i++) {
		auto f_it = mesh->face_handle(i);
		UpdateXOnePiece(f_it);
	}
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_f; i++) {
		auto f_it = mesh->face_handle(i);
		PassingXOnePiece(f_it);
	}
}

void FastSquarePieces2::UpdateZ()
{
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_f; i++) {
		auto f_it = mesh->face_handle(i);
		UpdateZOnePiece(f_it);
	}

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_f; i++) {
		auto f_it = mesh->face_handle(i);
		PassingZOnePiece(f_it);
	}
}

float FastSquarePieces2::UpdateU()
{
	float Residual = 0.0;
#ifdef _USE_OPENMP
	int coreNum = omp_get_num_procs();
	float* maxArray = new float[coreNum];
	for (int i = 0; i < coreNum; i++) {
		maxArray[i] = 0.0;
	}
#pragma omp parallel for
	for (int i = 0; i < num_f; i++) {
		auto f_it = mesh->face_handle(i);
		int k = omp_get_thread_num();
		float part_enrgy = UpdateUOnePiece(f_it);
		if (maxArray[k] < part_enrgy) {
			maxArray[k] = part_enrgy;
		}
	}
	for (int i = 0; i < coreNum; i++) {
		if (Residual < maxArray[i]) {
			Residual = maxArray[i];
		}
	
	}

	for (int i = 0; i < coreNum; i++) {
		maxArray[i] = 0.0;
	}
#pragma omp parallel for
	for (int i = 0; i < num_h; ++i) {
		auto hf_it = mesh->halfedge_handle(i);
		int k = omp_get_thread_num();
		if (mesh->property(W, hf_it)) {
			Vec3f dU = mesh->property(X, hf_it) - mesh->property(T, hf_it);
			float part_residual = (dU).max_abs();
			mesh->property(Ux, hf_it) = mesh->property(Ux, hf_it) + dU;
			if (maxArray[k] < part_residual) {
				maxArray[k] = part_residual;
}
		}
	}
	for (int i = 0; i < coreNum; i++) {
		if (Residual < maxArray[i]) {
			Residual = maxArray[i];
		}

	}
#else
	for (auto f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it) {
		float part_enrgy = UpdateUOnePiece(f_it);
		if (Residual < part_enrgy) {
			Residual = part_enrgy;
		}
	}
	for (auto hf_it = mesh->halfedges_begin(); hf_it != mesh->halfedges_end(); ++hf_it) {
		if (mesh->property(W, *hf_it)) {
			Vec3f dU = mesh->property(X, *hf_it) - mesh->property(T, *hf_it);
			float part_residual = (dU).max_abs();
			mesh->property(Ux, *hf_it) = mesh->property(Ux, *hf_it) + dU;
			if (Residual < part_residual) {
				Residual = part_residual;
			}
		}
	}
#endif

	return(Residual);
}

float FastSquarePieces2::UpdateOneStep()
{
	float residual = 0.0;
	UpdateX();
	UpdateZ();
	residual = UpdateU();
	return(residual);
}

void FastSquarePieces2::UpdateParameter()
{
	float rho_max = 1e5;
	int i = 0;
	clock_t start_t, end_t;
	float total_time = 0.0;
	float R = 10.0;
	vector<float> energy_set;
	vector<float> time_set;
	energy_set.push_back(EnergyFuntion());
	while (i < 20000 && R > 1e-2) {
		start_t = clock();
		R = UpdateOneStep();
#ifndef NDEBUG
		cout << "Iteration " << i << ": Residual: " << R << endl;
#endif // !NDEBUG

		//rho = rho * 2;
		AdjustBoundaryData();
		end_t = clock();
		total_time += (float)(end_t - start_t) / CLOCKS_PER_SEC;
		time_set.push_back(total_time);
		i++;
		energy_set.push_back(EnergyFuntion());
	}

	iter = i;

	ofstream energy_file("energy.txt");
	energy_file << "Energy:" << endl;
	for (int i = 0; i < energy_set.size(); i++) {
		energy_file << energy_set[i] << endl;
	}
	energy_file.close();

	ofstream time_file("time.txt");
	time_file << "Time:" << endl;
	for (int i = 0; i < time_set.size(); i++) {
		time_file << time_set[i] << endl;
	}
	time_file.close();

	solve_time = total_time;
}

void FastSquarePieces2::ReconstructSeed(PolyMesh::VertexHandle & v_it)
{
	for (PolyMesh::VertexVertexCWIter vv_it = mesh->vv_cwiter(v_it); vv_it.is_valid(); ++vv_it) {
		if (!mesh->property(R, *vv_it)) {
			mesh->property(R, *vv_it) = true;
			RecnstructOneEdge(v_it, *vv_it);
			ReconstructSeed(*vv_it);
		}
	}
}



void FastSquarePieces2::ReconstructBack()
{

	for (int i = 0; i < num_v; i++) {
		PolyMesh::VertexHandle v_it = mesh->vertex_handle(i);
		int x = mesh->point(v_it)[0] + J_x0 - I_x0;
		int y = mesh->point(v_it)[1] + J_y0 - I_y0;
		int r = change_number(mesh->property(V, v_it)[0]);
		int g = change_number(mesh->property(V, v_it)[1]);
		int b = change_number(mesh->property(V, v_it)[2]);
		J->setPixel(x, y, qRgb(r, g, b));
	}
}

float FastSquarePieces2::EnergyFuntion()
{
	float energy = 0.0;
	/*for (PolyMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); f_it++) {
		energy += EnergyofOnePiece(f_it);
	}*/

	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++) {
		float weight = 1.0;
		if (mesh->is_boundary(*e_it)) {
			weight = 0;
		}
		auto hf = mesh->halfedge_handle(*e_it, 0);
		energy += weight * (mesh->property(T, hf) - mesh->property(X, hf)).sqrnorm();
	}

	return(energy);
}

void FastSquarePieces2::ReconstructDFS()
{
	clock_t start_t, end_t;
	start_t = clock();
	mesh->property(V, first_handle) = color2vector(save_rgb);
	mesh->property(R, first_handle) = true;
	ReconstructSeed(first_handle);
	ReconstructBack();
	end_t = clock();
	cout << "The cost of reconstruction is: " << (float)(end_t - start_t) / CLOCKS_PER_SEC << endl;
	reconstruct_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}

void FastSquarePieces2::ReconstructBFS()
{
	clock_t start_t, end_t;
	start_t = clock();
	queue<PolyMesh::VertexHandle> vertex_queue;
	mesh->property(V, first_handle) = color2vector(save_rgb);
	mesh->property(R, first_handle) = true;
	vertex_queue.push(first_handle);
	while (!vertex_queue.empty()) {
		for (auto vv_it = mesh->vv_iter(vertex_queue.front()); vv_it.is_valid(); ++vv_it) {
			if (!mesh->property(R, *vv_it)) {
				//mesh->data(*vv_it).setR(true);
				RecnstructOneEdge(vertex_queue.front(), *vv_it);
				vertex_queue.push(*vv_it);
			}
		}
		vertex_queue.pop();
	}
	//cout << "BFS is over." << endl;
	ReconstructBack();
	end_t = clock();
	cout << "The cost of reconstruction is: " << (float)(end_t - start_t) / CLOCKS_PER_SEC << endl;
	reconstruct_time = (float)(end_t - start_t) / CLOCKS_PER_SEC;
}



void FastSquarePieces2::ReconstructFFT()
{
	double lambda = 1e7;
	double mu = 1e-5;
	//This part may be inefficiency......

	int M = J->height();
	int N = J->width();
	cv::Mat I_ = cv::Mat::zeros(M, N, CV_32FC3);

	cv::Mat otfFx_real = cv::Mat::zeros(M, N, CV_32F);
	cv::Mat otfFy_real = cv::Mat::zeros(M, N, CV_32F);

	float* ptr = otfFx_real.ptr<float>(0);
	ptr[0] = -1.0;
	ptr[N - 1] = 1.0;
	ptr = otfFy_real.ptr<float>(0);
	ptr[0] = -1.0;
	ptr = otfFy_real.ptr<float>(M - 1);
	ptr[0] = 1.0;

	cv::Mat win_x = (cv::Mat_<float>(1, 2) << -1.0, 1.0);
	cv::Mat win_y = (cv::Mat_<float>(2, 1) << -1.0, 1.0);

	cv::Mat S_x(M, N, CV_32FC3);
	cv::Mat S_y(M, N, CV_32FC3);

	cv::Mat otfFx[] = { cv::Mat_<float>(otfFx_real), cv::Mat::zeros(otfFx_real.size(), CV_32F) };
	cv::Mat otfFy[] = { cv::Mat_<float>(otfFy_real), cv::Mat::zeros(otfFy_real.size(), CV_32F) };

	cv::Mat otfFx_complex, otfFy_complex;
	cv::merge(otfFx, 2, otfFx_complex);
	cv::merge(otfFy, 2, otfFy_complex);

	cv::Mat f_otfFx_complex, f_otfFy_complex;
	cv::dft(otfFx_complex, f_otfFx_complex);
	cv::dft(otfFy_complex, f_otfFy_complex);

	cv::split(f_otfFx_complex, otfFx);
	cv::split(f_otfFy_complex, otfFy);

	otfFx_complex.release();
	otfFy_complex.release();
	otfFx_real.release();
	otfFy_real.release();
	f_otfFx_complex.release();
	f_otfFy_complex.release();
	//66,67,68th line, prepare for otfFx & otfFy

	//this part should replace

	cv::Mat I_rgb[] = { cv::Mat::zeros(M, N, CV_32F), cv::Mat::zeros(M, N, CV_32F), cv::Mat::zeros(M, N, CV_32F) };

	float *pI_r, *pI_g, *pI_b;
	for (int i = 0; i < M; i++) {
		pI_r = I_rgb[0].ptr<float>(i);
		pI_g = I_rgb[1].ptr<float>(i);
		pI_b = I_rgb[2].ptr<float>(i);
		for (int j = 0; j < N; j++) {
			pI_r[j] = (float)qRed(J->pixel(j, i));
			pI_g[j] = (float)qGreen(J->pixel(j, i));
			pI_b[j] = (float)qBlue(J->pixel(j, i));
		}
	}

	//adjust image

	/*for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
		I_rgb[0].at<float>(mesh->point(*v_it)[1] + J_y0 - I_y0, mesh->point(*v_it)[0] + J_x0 - I_x0) = (float)qRed(I->pixel(mesh->point(*v_it)[0], mesh->point(*v_it)[1]));
		I_rgb[1].at<float>(mesh->point(*v_it)[1] + J_y0 - I_y0, mesh->point(*v_it)[0] + J_x0 - I_x0) = (float)qGreen(I->pixel(mesh->point(*v_it)[0], mesh->point(*v_it)[1]));
		I_rgb[2].at<float>(mesh->point(*v_it)[1] + J_y0 - I_y0, mesh->point(*v_it)[0] + J_x0 - I_x0) = (float)qBlue(I->pixel(mesh->point(*v_it)[0], mesh->point(*v_it)[1]));
	}*/

	merge(I_rgb, 3, I_);



	cv::Mat I_r[] = { cv::Mat_<float>(I_rgb[0]), cv::Mat::zeros(I_rgb[0].size(), CV_32F) };
	cv::Mat I_g[] = { cv::Mat_<float>(I_rgb[1]), cv::Mat::zeros(I_rgb[1].size(), CV_32F) };
	cv::Mat I_b[] = { cv::Mat_<float>(I_rgb[2]), cv::Mat::zeros(I_rgb[2].size(), CV_32F) };


	cv::Mat I_r_complex, I_g_complex, I_b_complex;
	merge(I_r, 2, I_r_complex);
	merge(I_g, 2, I_g_complex);
	merge(I_b, 2, I_b_complex);

	I_r[0].release();
	I_g[0].release();
	I_b[0].release();
	I_r[1].release();
	I_g[1].release();
	I_b[1].release();

	cv::Mat f_I_r_complex, f_I_g_complex, f_I_b_complex;
	cv::dft(I_r_complex, f_I_r_complex);
	cv::dft(I_g_complex, f_I_g_complex);
	cv::dft(I_b_complex, f_I_b_complex);

	I_r_complex.release();
	I_g_complex.release();
	I_b_complex.release();

	//69th line
	cv::Mat Denormin_real = otfFx[0].mul(otfFx[0]) + otfFx[1].mul(otfFx[1]) + otfFy[0].mul(otfFy[0]) + otfFy[1].mul(otfFy[1]);

	otfFx[0].release();
	otfFy[0].release();
	otfFx[1].release();
	otfFy[1].release();

	Denormin_real = /*1 + */lambda * Denormin_real;

	cv::Mat S_x_exp, S_y_exp;

	filter2D(I_, S_x, I_.depth(), win_x, cv::Point(0, 0));
	filter2D(I_, S_y, I_.depth(), win_y, cv::Point(0, 0));

	hconcat(S_x.colRange(0, N - 1), I_.col(0).clone() - I_.col(N - 1).clone(), S_x_exp);
	vconcat(S_y.rowRange(0, M - 1), I_.row(0).clone() - I_.row(M - 1).clone(), S_y_exp);

	for (auto f_it = mesh->faces_begin(); f_it != mesh->faces_end(); f_it++) {
		gMARK token = information4face(f_it);

		S_x_exp.at<cv::Vec3f>(token.y + J_y0 - I_y0, token.x + J_x0 - I_x0)[0] = token.g_x[0];
		S_x_exp.at<cv::Vec3f>(token.y + J_y0 - I_y0, token.x + J_x0 - I_x0)[1] = token.g_x[1];
		S_x_exp.at<cv::Vec3f>(token.y + J_y0 - I_y0, token.x + J_x0 - I_x0)[2] = token.g_x[2];

		S_y_exp.at<cv::Vec3f>(token.y + J_y0 - I_y0, token.x + J_x0 - I_x0)[0] = token.g_y[0];
		S_y_exp.at<cv::Vec3f>(token.y + J_y0 - I_y0, token.x + J_x0 - I_x0)[1] = token.g_y[1];
		S_y_exp.at<cv::Vec3f>(token.y + J_y0 - I_y0, token.x + J_x0 - I_x0)[2] = token.g_y[2];


	}

	//75,76th line

	cv::Mat Normin = cv::Mat::zeros(M, N, CV_32FC3);
	Normin.col(0) = Normin.col(0) + S_x_exp.col(N - 1) - S_x_exp.col(0);
	for (int i = 1; i < N; i++) {
		Normin.col(i) = Normin.col(i) + S_x_exp.col(i - 1) - S_x_exp.col(i);
	}
	Normin.row(0) = Normin.row(0) + S_y_exp.row(M - 1) - S_y_exp.row(0);
	for (int j = 1; j < M; j++) {
		Normin.row(j) = Normin.row(j) + S_y_exp.row(j - 1) - S_y_exp.row(j);
	}

	S_x_exp.release();
	S_y_exp.release();

	//77,78th line
	cv::Mat S_rgb[3];
	split(Normin, S_rgb);

	Normin.release();

	cv::Mat S_r[] = { cv::Mat_<float>(S_rgb[0]), cv::Mat::zeros(S_rgb[0].size(), CV_32F) };
	cv::Mat S_g[] = { cv::Mat_<float>(S_rgb[1]), cv::Mat::zeros(S_rgb[1].size(), CV_32F) };
	cv::Mat S_b[] = { cv::Mat_<float>(S_rgb[2]), cv::Mat::zeros(S_rgb[2].size(), CV_32F) };

	S_rgb[0].release();
	S_rgb[1].release();
	S_rgb[2].release();

	cv::Mat S_r_complex, S_g_complex, S_b_complex;
	merge(S_r, 2, S_r_complex);
	merge(S_g, 2, S_g_complex);
	merge(S_b, 2, S_b_complex);

	S_r[0].release();
	S_g[0].release();
	S_b[0].release();
	S_r[1].release();
	S_g[1].release();
	S_b[1].release();

	cv::Mat f_S_r_complex, f_S_g_complex, f_S_b_complex;
	dft(S_r_complex, f_S_r_complex);
	dft(S_g_complex, f_S_g_complex);
	dft(S_b_complex, f_S_b_complex);




	S_r_complex.release();
	S_g_complex.release();
	S_b_complex.release();

	//prepare for fft2(Normin2)
	cv::Mat FS_r_complex = (/*f_I_r_complex + */lambda * f_S_r_complex);
	cv::Mat FS_g_complex = (/*f_I_g_complex + */lambda * f_S_g_complex);
	cv::Mat FS_b_complex = (/*f_I_b_complex + */lambda * f_S_b_complex);

	f_S_r_complex.release();
	f_S_g_complex.release();
	f_S_b_complex.release();
	f_I_r_complex.release();
	f_I_g_complex.release();
	f_I_b_complex.release();

	cv::Mat FS_r[2], FS_g[2], FS_b[2];
	split(FS_r_complex, FS_r);
	split(FS_g_complex, FS_g);
	split(FS_b_complex, FS_b);

	cv::Mat new_FS_r_r = FS_r[0].mul(1 / Denormin_real);
	cv::Mat new_FS_g_r = FS_g[0].mul(1 / Denormin_real);
	cv::Mat new_FS_b_r = FS_b[0].mul(1 / Denormin_real);
	cv::Mat new_FS_r_i = FS_r[1].mul(1 / Denormin_real);
	cv::Mat new_FS_g_i = FS_g[1].mul(1 / Denormin_real);
	cv::Mat new_FS_b_i = FS_b[1].mul(1 / Denormin_real);

	new_FS_r_r.copyTo(FS_r[0]);
	new_FS_g_r.copyTo(FS_g[0]);
	new_FS_b_r.copyTo(FS_b[0]);
	new_FS_r_i.copyTo(FS_r[1]);
	new_FS_g_i.copyTo(FS_g[1]);
	new_FS_b_i.copyTo(FS_b[1]);

	merge(FS_r, 2, FS_r_complex);
	merge(FS_g, 2, FS_g_complex);
	merge(FS_b, 2, FS_b_complex);

	FS_r[0].release();
	FS_g[0].release();
	FS_b[0].release();
	FS_r[1].release();
	FS_g[1].release();
	FS_b[1].release();

	//79th line
	cv::Mat if_FS_r_complex, if_FS_g_complex, if_FS_b_complex;
	idft(FS_r_complex, if_FS_r_complex, cv::DFT_REAL_OUTPUT);
	idft(FS_g_complex, if_FS_g_complex, cv::DFT_REAL_OUTPUT);
	idft(FS_b_complex, if_FS_b_complex, cv::DFT_REAL_OUTPUT);

	FS_r_complex.release();
	FS_g_complex.release();
	FS_b_complex.release();

	cv::Mat S_r_r = if_FS_r_complex / (M*N);
	cv::Mat S_g_r = if_FS_g_complex / (M*N);
	cv::Mat S_b_r = if_FS_b_complex / (M*N);

	I_rgb[0] = cv::Mat_<float>(S_r_r.clone());
	I_rgb[1] = cv::Mat_<float>(S_g_r.clone());
	I_rgb[2] = cv::Mat_<float>(S_b_r.clone());

	S_r_r.release();
	S_g_r.release();
	S_b_r.release();

	merge(I_rgb, 3, I_);
	//80th line

	float dr = (float)qRed(J->pixel(0, 0)) - I_rgb[0].at<float>(0, 0);
	float dg = (float)qGreen(J->pixel(0, 0)) - I_rgb[1].at<float>(0, 0);
	float db = (float)qBlue(J->pixel(0, 0)) - I_rgb[2].at<float>(0, 0);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			int r = change_number(I_rgb[0].at<float>(i, j) + dr);
			int g = change_number(I_rgb[1].at<float>(i, j) + dg);
			int b = change_number(I_rgb[2].at<float>(i, j) + db);

			QRgb pixel_ptr = qRgb(r, g, b);
			J->setPixel(j, i, pixel_ptr);
		}
	}

	I_rgb[0].release();
	I_rgb[1].release();
	I_rgb[2].release();

	I_.release();
}

float FastSquarePieces2::BoundaryResidual()
{
	float residual = 0.0;
	for (auto hf_it = mesh->halfedges_begin(); hf_it != mesh->halfedges_end(); hf_it++) {
		if (mesh->is_boundary(*hf_it)) {
			residual += (mesh->property(X, *hf_it) - mesh->property(T, *hf_it)).sqrnorm();
		}
	}
	return residual;
}

Vec3f FastSquarePieces2::color2vector(QRgb color)
{
	float r = (float)qRed(color);
	float g = (float)qGreen(color);
	float b = (float)qBlue(color);
	return Vec3f(r, g, b);
}

void FastSquarePieces2::CalculateCurl(Vec3f & curl, PolyMesh::HalfedgeHandle hf)
{
	switch (mesh->property(L, hf)) {
	case L_Location:
		curl += mesh->property(X, hf);
		break;
	case U_Location:
		curl += mesh->property(X, hf);
		break;
	case R_Location:
		curl -= mesh->property(X, hf);
		break;
	case D_Location:
		curl -= mesh->property(X, hf);
		break;
	}
}

float FastSquarePieces2::EnergyofOnePiece(PolyMesh::FaceIter f_it)
{
	float energy_1 = 0.0;
	float energy_2 = 0.0;
	Vec3f curl = Vec3f(0.0, 0.0, 0.0);

	PolyMesh::FaceHalfedgeIter fh_it = mesh->fh_begin(f_it);
	energy_1 += (mesh->property(W, *fh_it) ? 0.0 : 1.0)*(mesh->property(T, *fh_it) - mesh->property(X, *fh_it)).sqrnorm();
	CalculateCurl(curl, *fh_it);
	fh_it++;
	energy_1 += (mesh->property(W, *fh_it) ? 0.0 : 1.0)*(mesh->property(T, *fh_it) - mesh->property(X, *fh_it)).sqrnorm();
	CalculateCurl(curl, *fh_it);
	fh_it++;
	energy_1 += (mesh->property(W, *fh_it) ? 0.0 : 1.0)*(mesh->property(T, *fh_it) - mesh->property(X, *fh_it)).sqrnorm();
	CalculateCurl(curl, *fh_it);
	fh_it++;
	energy_1 += (mesh->property(W, *fh_it) ? 0.0 : 1.0)*(mesh->property(T, *fh_it) - mesh->property(X, *fh_it)).sqrnorm();
	CalculateCurl(curl, *fh_it);

	energy_2 = curl.sqrnorm();
	float energy = energy_1 + energy_2;
	return(energy_1);
}

void FastSquarePieces2::RecnstructOneEdge(PolyMesh::VertexHandle & v_it, PolyMesh::VertexHandle & u_it)
{
	HalfedgeHandle heh = mesh->find_halfedge(v_it, u_it);
	mesh->property(R, u_it) = true;

	int d_x = mesh->point(u_it)[0] - mesh->point(v_it)[0];
	int d_y = mesh->point(u_it)[1] - mesh->point(v_it)[1];
	if (d_x*d_y != 0) {
#ifndef NDEBUG
		cout << "This is wrong" << endl;
#endif // !NDEBUG
		abort();
	}
	if (d_x == 1 || d_y == 1) {
		mesh->property(V, u_it) = mesh->property(V, v_it) + mesh->property(X, heh);
	}
	if (d_x == -1 || d_y == -1) {
		mesh->property(V, u_it) = mesh->property(V, v_it) - mesh->property(X, heh);
	}
	/*
	switch (mesh->data(heh).getLocation())
	{
	case L_Location:
		mesh->data(u_it).setV(mesh->data(v_it).getV() + mesh->data(heh).getX());
		break;
	case R_Location:
		mesh->data(u_it).setV(mesh->data(v_it).getV() - mesh->data(heh).getX());
		break;
	case U_Location:
		mesh->data(u_it).setV(mesh->data(v_it).getV() + mesh->data(heh).getX());
		break;
	case D_Location:
		mesh->data(u_it).setV(mesh->data(v_it).getV() - mesh->data(heh).getX());
	default:
		break;
	}
	*/
}

void FastSquarePieces2::AdjustofHalfedge()
{
	for (auto he_it = mesh->halfedges_begin(); he_it != mesh->halfedges_end(); he_it++) {
		//add property
		mesh->property(Ux, *he_it).vectorize(0.0f);
		if (mesh->is_boundary(*he_it)) {
			PolyMesh::HalfedgeHandle opposite_heh = mesh->opposite_halfedge_handle(*he_it);
			switch (mesh->property(L, opposite_heh))
			{
			case L_Location:
				mesh->property(L, *he_it) = R_Location;
				mesh->property(W, *he_it) = true;
				//mesh->data(opposite_heh).setW(true);
				break;
			case R_Location:
				mesh->property(L, *he_it) = L_Location;
				mesh->property(W, *he_it) = true;
				//mesh->data(opposite_heh).setW(true);
				break;
			case U_Location:
				mesh->property(L, *he_it) = D_Location;
				mesh->property(W, *he_it) = true;
				//mesh->data(opposite_heh).setW(true);
				break;
			case D_Location:
				mesh->property(L, *he_it) = U_Location;
				mesh->property(W, *he_it) = true;
				//mesh->data(opposite_heh).setW(true);
				break;
			default:
				break;
			}
		}
	}
}

void FastSquarePieces2::AdjustBoundaryData()
{
	for (auto he_it = mesh->halfedges_begin(); he_it != mesh->halfedges_end(); he_it++) {
		if (mesh->is_boundary(*he_it)) {
			PolyMesh::HalfedgeHandle opposite_heh = mesh->opposite_halfedge_handle(*he_it);
			mesh->property(T, *he_it) = mesh->property(T, opposite_heh);
			mesh->property(X, *he_it) = mesh->property(X, opposite_heh);
			mesh->property(U, *he_it) = mesh->property(U, opposite_heh);
			mesh->property(Z, *he_it) = mesh->property(Z, opposite_heh);
			mesh->property(Ux, *he_it) = mesh->property(Ux, opposite_heh);
		}
	}
}

gMARK FastSquarePieces2::information4face(const PolyMesh::FaceIter & f_it)
{
	gMARK token;

	PolyMesh::Point start_point = mesh->point(*mesh->fv_begin(*f_it));

	int x_min_ = start_point[0];
	int y_min_ = start_point[1];
	PolyMesh::FaceVertexIter v_it = mesh->fv_begin(*f_it);

	for (auto fv_it = mesh->fv_begin(*f_it); fv_it != mesh->fv_end(*f_it); fv_it++) {
		if (mesh->point(*fv_it)[0] < x_min_) {
			x_min_ = mesh->point(*fv_it)[0];
			v_it = fv_it;
		}
		if (mesh->point(*fv_it)[1] < y_min_) {
			y_min_ = mesh->point(*fv_it)[1];
			v_it = fv_it;
		}
	}

	token.x = x_min_;
	token.y = y_min_;
	token.boundary = mesh->is_boundary(*v_it);

	for (auto u_it = mesh->vv_begin(*v_it); u_it != mesh->vv_end(*v_it); u_it++) {
		if (mesh->point(*u_it)[0] > x_min_) {
			token.g_x = mesh->property(X, mesh->find_halfedge(v_it, u_it));
		}

		if (mesh->point(*u_it)[1] > y_min_) {
			token.g_y = mesh->property(X, mesh->find_halfedge(v_it, u_it));
		}
	}

	return token;
}

void FastSquarePieces2::BoundaryCheck()
{
	bool token = true;

	for (auto e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++) {
		if (mesh->is_boundary(*e_it)) {
			if (mesh->property(W, mesh->halfedge_handle(*e_it, 0)) == false || mesh->property(W, mesh->halfedge_handle(*e_it, 1)) == false) {
				token = false;
				if (mesh->property(W, mesh->halfedge_handle(*e_it, 0)) == false && mesh->property(W, mesh->halfedge_handle(*e_it, 1)) == false) {
					cout << "both of boundary wrong" << endl;
				}
			}
		}
		else {
			if (mesh->property(W, mesh->halfedge_handle(*e_it, 0)) == true || mesh->property(W, mesh->halfedge_handle(*e_it, 1)) == true) {
				token = false;
				cout << "inner wrong" << endl;
			}
		}
	}

	if (token) {
		cout << "Correct in boundary edge rigist!" << endl;
	}
	else {
		cout << "Wrong in boundary edge rigist!" << endl;
	}
}


void FastSquarePieces2::RegistMesh(bool token)
{
	clock_t start_t, end_t;
	//add vertex in mesh
	start_t = clock();
	vector<vector<PolyMesh::VertexHandle>> handle_table;
	for (int i = 0; i < inner_mark.size(); i++) {
		vector<PolyMesh::VertexHandle> inner_line;
		for (int j = 0; j < inner_mark[i].size(); j++) {
			PolyMesh::VertexHandle vh = mesh->add_vertex(PolyMesh::Point(inner_mark[i][j].x, inner_mark[i][j].y));
			//mesh->data(vh).init();
			//mesh->data(vh).setB(inner_mark[i][j].is_boundary);
			inner_line.push_back(vh);
		}
		handle_table.push_back(inner_line);
	}

	//seed handle
	first_handle = handle_table[0][0];
	//regist faces
	int square_num = 0;
	for (int i = 1; i < inner_mark.size() - 1; i++) {
		for (int j = 1; j < inner_mark[i].size() - 1; j++) {
			if (!inner_mark[i][j].is_boundary) {
				//order is 00->01->11->10
				//basic face
				vector<PolyMesh::VertexHandle> face_vhandles;
				PolyMesh::FaceHandle face_handle;
				//searching for upper and lower vertex
				int k = 0;
				for (int k_ = 0; k_ < inner_mark[i - 1].size(); k_++) {
					if (inner_mark[i][j].x == inner_mark[i - 1][k_].x) {
						k = k_;
					}
				}
				int l = 0;
				for (int l_ = 0; l_ < inner_mark[i + 1].size(); l_++) {
					if (inner_mark[i][j].x == inner_mark[i + 1][l_].x) {
						l = l_;
					}
				}
				//Basic Square Pieces
				face_vhandles.clear();
				face_vhandles.push_back(handle_table[i][j]);
				face_vhandles.push_back(handle_table[i - 1][k]);
				face_vhandles.push_back(handle_table[i - 1][k + 1]);
				face_vhandles.push_back(handle_table[i][j + 1]);
				face_handle = mesh->add_face(face_vhandles);
				RegistSquarePiece(face_handle, inner_mark[i - 1][k], inner_mark[i][j], inner_mark[i][j + 1], inner_mark[i - 1][k + 1], token);
				square_num++;

				//Left Square Pieces
				if (inner_mark[i][j - 1].is_boundary) {
					face_vhandles.clear();
					face_vhandles.push_back(handle_table[i][j - 1]);
					face_vhandles.push_back(handle_table[i - 1][k - 1]);
					face_vhandles.push_back(handle_table[i - 1][k]);
					face_vhandles.push_back(handle_table[i][j]);
					face_handle = mesh->add_face(face_vhandles);
					RegistSquarePiece(face_handle, inner_mark[i - 1][k - 1], inner_mark[i][j - 1], inner_mark[i][j], inner_mark[i - 1][k], token);
					square_num++;
				}

				//Down Square Pieces
				if (inner_mark[i + 1][l].is_boundary&&inner_mark[i + 1][l + 1].is_boundary) {
					face_vhandles.clear();
					face_vhandles.push_back(handle_table[i + 1][l]);
					face_vhandles.push_back(handle_table[i][j]);
					face_vhandles.push_back(handle_table[i][j + 1]);
					face_vhandles.push_back(handle_table[i + 1][l + 1]);
					face_handle = mesh->add_face(face_vhandles);
					RegistSquarePiece(face_handle, inner_mark[i][j], inner_mark[i + 1][l], inner_mark[i + 1][l + 1], inner_mark[i][j + 1], token);
					square_num++;
				}

				//Left-Down Square Pieces
				if (inner_mark[i + 1][l].is_boundary&&inner_mark[i + 1][l - 1].is_boundary&&inner_mark[i][j - 1].is_boundary) {
					face_vhandles.clear();
					face_vhandles.push_back(handle_table[i + 1][l - 1]);
					face_vhandles.push_back(handle_table[i][j - 1]);
					face_vhandles.push_back(handle_table[i][j]);
					face_vhandles.push_back(handle_table[i + 1][l]);
					face_handle = mesh->add_face(face_vhandles);
					RegistSquarePiece(face_handle, inner_mark[i][j - 1], inner_mark[i + 1][l - 1], inner_mark[i + 1][l], inner_mark[i][j], token);
					square_num++;
				}
			}
		}
	}
	auto v_it_ = mesh->vertices_end();
	v_it_--;
	num_v = v_it_->idx() + 1;
	auto e_it_ = mesh->edges_end();
	e_it_--;
	num_e = e_it_->idx() + 1;
	auto f_it_ = mesh->faces_end();
	f_it_--;
	num_f = f_it_->idx() + 1;
	auto h_it_ = mesh->halfedges_end();
	h_it_--;
	num_h = h_it_->idx() + 1;


#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_h; i++) {
		mesh->property(W, mesh->halfedge_handle(i)) = false;
	}


	
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < num_v; i++) {
		mesh->property(R, mesh->vertex_handle(i)) = false;
	}
	end_t = clock();
	//cout << "The time of regist mesh is: " << (double)(end_t - start_t) / CLOCKS_PER_SEC << "s." << endl;
}

Vec3f FastSquarePieces2::CalculateDUt(PolyMesh::HalfedgeHandle hf)
{
	Vec3f curl;
	switch (mesh->property(L, hf)) {
	case L_Location:
		curl = -mesh->property(X, hf);
		break;
	case U_Location:
		curl = -mesh->property(Z, hf);
		break;
	case R_Location:
		curl = mesh->property(Z, hf);
		break;
	case D_Location:
		curl = mesh->property(X, hf);
		break;
	}
	return(curl);
}