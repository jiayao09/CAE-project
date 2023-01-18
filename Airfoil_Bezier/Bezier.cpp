#include "Header.h"
#include <cmath>
#include <iostream>

void Basisfun(const float ub, const float* knot, int span, int k, float* basis_matrix) {
	float* left = new float[k];
	float* right = new float[k];
	float temp = 0;

	for (int j = 0; j < k; j++) {
		left[j] = 0;
		right[j] = 0;
	}
	basis_matrix[0] = 1.0;
	for (int j = 1; j < k; j++) {
		left[j] = ub - knot[span - j + 1];
		right[j] = knot[span + j] - ub;
		float saved = 0.0;
		for (int r = 0; r < j; r++) {
			temp = basis_matrix[r] / (right[r + 1] + left[j - r]);
			basis_matrix[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		basis_matrix[j] = saved;

	}
}
void InterForSurfU(const float* Data_point, int k, int n, const float* ub, const float* knot, float* control_point, int order) {
	float* Amatrix = new float[n * n];
	float* basis_matrix = new float[n];
	float* D_matrix = new float[n];
	float* c_matrix = new float[n];
	for (int i = 0; i < n; i++) {
		basis_matrix[i] = 0;
		D_matrix[i] = Data_point[order * n + i];
		c_matrix[i] = 0;
	}
	for (int i = 1; i < n - 1; i++) {
		int span = 0;
		while (ub[i] > knot[span])
		{
			span += 1;
		}
		span = span - 1;

		Basisfun(ub[i], knot, span, k, basis_matrix);

		for (int j = 0; j < n; j++) {
			if (j > span - k && j < span + 1) {
				Amatrix[i * n + j] = basis_matrix[j - span + k - 1];
			}
			else {
				Amatrix[i * n + j] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		Amatrix[i] = 0;
		Amatrix[n * (n - 1) + i] = 0;
	}
	Amatrix[0] = 1;
	Amatrix[n * n - 1] = 1;

	linearslove(Amatrix, D_matrix, c_matrix, n, n, order);
	for (int i = 0; i < n; i++) {
		control_point[order * n + i] = c_matrix[i];
	}
	delete[] D_matrix;
	delete[] c_matrix;
	delete[] basis_matrix;
	delete[] Amatrix;
}

void InterForSurfV(const float* Data_point, int k, int n, const float* ub, const float* knot, float* control_point, int order) {
	float* Amatrix = new float[n * n];
	float* basis_matrix = new float[n];
	float* D_matrix = new float[n];
	float* c_matrix = new float[n];
	for (int i = 0; i < n; i++) {
		D_matrix[i] = Data_point[i * (SECTIONNO + EXTRAP) + order];
		c_matrix[i] = 0;
		basis_matrix[i] = 0;
	}
	for (int i = 1; i < n - 1; i++) {
		int span = 0;
		while (ub[i] > knot[span])
		{
			span += 1;
		}
		span = span - 1;

		Basisfun(ub[i], knot, span, k, basis_matrix);

		for (int j = 0; j < n; j++) {
			if (j > span - k && j < span + 1) {
				Amatrix[i * n + j] = basis_matrix[j - span + k - 1];
			}
			else {
				Amatrix[i * n + j] = 0;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		Amatrix[i] = 0;
		Amatrix[n * (n - 1) + i] = 0;
	}
	Amatrix[0] = 1;
	Amatrix[n * n - 1] = 1;

	linearslove(Amatrix, D_matrix, c_matrix, n, n, order);

		for (int i = 0; i < n; i++) {
			control_point[i * (SECTIONNO + EXTRAP) + order] = c_matrix[i];
		}	

	delete[] D_matrix;
	delete[] c_matrix;
	delete[] basis_matrix;
	delete[] Amatrix;
}





void  BsplineSurfInter(const int* input, const float* Dx, const float* Dy, const float* Dz, float* Control_pointx, float* Control_pointy, float* Control_pointz, float* knotU, float* knotV) {
	// ku DEGREEU, kv DEGREEV, nu SECTIONNO, nv section
	// Set up input Data point 
	int nu = SECTIONNO + EXTRAP;
	int nv = input[0];
	float* Data_pointx = new float[nu * nv];
	float* Data_pointy = new float[nu * nv];
	float* Data_pointz = new float[nu * nv];

	for (int i = 0; i < input[0]; i++) {
		for (int j = 0; j < EXTRAP; j++) {
			Data_pointx[i * nu + j] = Dx[(i + 1) * SECTIONNO - 1] - Dx[(i + 1) * SECTIONNO - 1] * j / EXTRAP;
			Data_pointy[i * nu + j] = 0;
			Data_pointz[i * nu + j] = Dz[i] / 2 * 1000;

		}
		for (int j = 0; j < SECTIONNO; j++) {
			Data_pointx[i * nu + j + EXTRAP] = Dx[i * SECTIONNO + j];
			Data_pointy[i * nu + j + EXTRAP] = Dy[i * SECTIONNO + j];
			Data_pointz[i * nu + j + EXTRAP] = Dz[i] / 2 * 1000;
		}
	}

	// Step 1: Computing parameters: ubu(), ubv()
	float* ubu = new float[nu];
	float* ubv = new float[nv];
	float* cdu = new float[nu];
	float* cdv = new float[nv];
	for (int i = 0; i < nu; i++) {
		ubu[i] = 0;
		cdu[i] = 0;
	}
	for (int i = 0; i < nv; i++) {
		ubv[i] = 0;
		cdv[i] = 0;
	}

	for (int i = 0; i < nv; i++) {
		float totalChord = 0;
		for (int j = 0; j < nu - 1; j++) {
			float deltax = (Data_pointx[i * nu + j + 1] - Data_pointx[i * nu + j]) * (Data_pointx[i * nu + j + 1] - Data_pointx[i * nu + j]);
			float deltay = (Data_pointy[i * nu + j + 1] - Data_pointy[i * nu + j]) * (Data_pointy[i * nu + j + 1] - Data_pointy[i * nu + j]);
			float deltaz = (Data_pointz[i * nu + j + 1] - Data_pointz[i * nu + j]) * (Data_pointz[i * nu + j + 1] - Data_pointz[i * nu + j]);
			cdu[j] = pow(deltax + deltay + deltaz, 0.5);
			totalChord += cdu[j];
		}
		float dt = 0;
		for (int j = 0; j < nu - 1; j++) {
			dt += cdu[j];
			ubu[j + 1] += dt / totalChord;
		}

	}

	for (int i = 1; i < nu; i++) {
		ubu[i] = ubu[i] / nv;
	}
	ubu[nu - 1] = 1.0;

	for (int i = 0; i < nu; i++) {
		float totalChordv = 0;
		for (int j = 0; j < nv - 1; j++) {
			float deltax = (Data_pointx[i + nu * (j + 1)] - Data_pointx[i + nu * j]) * (Data_pointx[i + nu * (j + 1)] - Data_pointx[i + nu * j]);
			float deltay = (Data_pointy[i + nu * (j + 1)] - Data_pointy[i + nu * j]) * (Data_pointy[i + nu * (j + 1)] - Data_pointy[i + nu * j]);
			float deltaz = (Data_pointz[i + nu * (j + 1)] - Data_pointz[i + nu * j]) * (Data_pointz[i + nu * (j + 1)] - Data_pointz[i + nu * j]);
			cdv[j] = pow(deltax + deltay + deltaz, 0.5);
			totalChordv += cdv[j];
		}
		float dt = 0;
		for (int j = 0; j < nv - 1; j++) {
			dt += cdv[j];
			ubv[j + 1] += dt / totalChordv;
		}
	}

	for (int i = 1; i < nv; i++) {
		ubv[i] = ubv[i] / nu;
	}
	ubv[nv - 1] = 1.0;


	// Step 2: Computing knot vector : U, V
	for (int i = 0; i < nu - DEGREEU; i++) {
		knotU[i + DEGREEU] = 0.;
		for (int j = i; j < i + DEGREEU - 1; j++) {
			knotU[i + DEGREEU] += ubu[j + 1];
		}
		knotU[i + DEGREEU] = knotU[i + DEGREEU] / (DEGREEU - 1);
	}
	for (int i = 0; i < DEGREEU; i++) {
		knotU[i] = 0;
	}
	for (int i = nu; i < DEGREEU + nu; i++) {
		knotU[i] = 1.0;
	}

	for (int i = 0; i < nv - DEGREEV; i++) {
		knotV[i + DEGREEV] = 0.;
		for (int j = i; j < i + DEGREEV - 1; j++) {
			knotV[i + DEGREEV] += ubv[j + 1];
		}
		knotV[i + DEGREEV] = knotV[i + DEGREEV] / (DEGREEV - 1);
	}
	for (int i = 0; i < DEGREEV; i++) {
		knotV[i] = 0;
	}
	for (int i = nv; i < DEGREEV + nv; i++) {
		knotV[i] = 1.0;
	}

	float* Qx = new float[nv * nu];
	float* Qy = new float[nv * nu];
	float* Qz = new float[nv * nu];
	for (int i = 0; i < nv * nu; i++) {
		Qx[i] = 0.0;
		Qy[i] = 0.0;
		Qz[i] = 0.0;
	}
	// Step 3: Get the intermediate control points Q
	for (int i = 0; i < nv; i++) {
		int k = DEGREEU;
		InterForSurfU(Data_pointx, k, nu, ubu, knotU, Qx, i);
		InterForSurfU(Data_pointy, k, nu, ubu, knotU, Qy, i);
		InterForSurfU(Data_pointz, k, nu, ubu, knotU, Qz, i);
	}

	// Step 4: Get the final Control points P

	for (int i = 0; i < nu; i++) {
		int k = DEGREEV;
		InterForSurfV(Qx, k, nv, ubv, knotV, Control_pointx, i);
		InterForSurfV(Qy, k, nv, ubv, knotV, Control_pointy, i);
		InterForSurfV(Qz, k, nv, ubv, knotV, Control_pointz, i);
	}
	
	delete[] Data_pointx;
	delete[] Data_pointy;
	delete[] Data_pointz;
	delete[] Qx;
	delete[] Qy;
	delete[] Qz;
	delete[] ubu;
	delete[] ubv;
	delete[] cdu;
	delete[] cdv;
}
