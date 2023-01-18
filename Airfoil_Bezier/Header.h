#pragma once
#include <cmath>
#include <string>
#include <iostream>
#include <omp.h>
#include <fstream>
#include "C:\Users\pangj\OneDrive\Desktop\project\FinalProject\eigen-3.4.0\Eigen\Dense"
using namespace Eigen;
using namespace std;

#define PI 3.14159265
#define FACTORA 0.1
#define FACTORB 0.01
#define PROPCOEFFA 0.5
#define PROPCOEFFB 0.5
#define NUMBERP 40 // # of points on the B-spline curve
#define SECTIONNO 14
#define DEGREEU 4 //order of the curve u
#define DEGREEV 3 // order of the curve v
#define EXTRAP 11

void cal_sectionblade(const int n_s, int* output_sec);
void cal_coeffab(const int* input, const float capacity, const float head, const int n, const int n_s, int l_axial, float* coeff);
void bladeangle(const int* input, const float* coeff, const int n, const float capacity, const float head, float* output_angle);
void cal_chordlength(const int* input, const float* coeff, const float capacity, float* l_chord);
void cal_impellerDia(const int* input, const float* coeff, const float* l_chord, const float capacity, const int n_s, float* impeller_dia);
void airfoil_791(const int* input, float* impeller_dia, float head, const float* l_chord, float* airfoil_thicknessmax, float* airfoil_thickness, float* airfoil_length);


void Basisfun(const float ub, const float* knot, int span, int k, float* basis_matrix);
void InterForSurfV(const float* Data_point, int k, int n, const float* ub, const float* knot, float* control_point, int order);
void InterForSurfU(const float* Data_point, int k, int n, const float* ub, const float* knot, float* control_point, int order);
void BsplineSurfInter(const int* input, const float* Dx, const float* Dy, const float* Dz, float* Control_pointx, float* Control_pointy, float* Control_pointz, float* knotU, float* knotV);

void linearslove(const float* A, const float* B, float* C, int x, int y, int order);
void disp(const float* array, int size, int y);
void disp(const float* array, int size);

void writefile(const int* input, float* Control_pointx, float* Control_pointy, float* Control_pointz, float* knotU, float* knotV, int dvn);