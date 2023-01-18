#include "Header.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>

using namespace std;

int main(int argc, char** argv) {
	
	// Get input basic information
	// 1. check IGS file, and input number
	int dvn = 10; // resolution
	float capacity = 0.35;      // m^3/s
	float head = 6.72;	// m
	int rotary_speed = 1450;	// r/min
	int l_axial = 26;	//mm
	int n_s = (3.65 * rotary_speed * sqrt(capacity)) / (pow(head, 0.75));

	int* sec_blade = new int[2];
	cal_sectionblade(n_s, sec_blade);	// Determine the number of setion and number of blades
	if(sec_blade[0] == 0 || sec_blade[1] == 0) { 
		std::cout << "Error in input number" << std::endl; 
		exit(0);
	}

	float* coeff = new float[2];
	float* l_chord = new float[sec_blade[0]];
	float* blade_angle = new float[sec_blade[0]];
	cal_coeffab(sec_blade, capacity, head, rotary_speed, n_s, l_axial, coeff);
	cal_chordlength(sec_blade, coeff, capacity, l_chord);	// caculate the airfoil chord length 
	bladeangle(sec_blade, coeff, rotary_speed, capacity, head, blade_angle); // caculate the blade angle
	
	float* impeller_dia= new float[sec_blade[0]];
	cal_impellerDia(sec_blade, coeff, l_chord, capacity, n_s, impeller_dia);
	float hub_ratio = -0.00029 * n_s + 0.71587;
	float hub_dia = hub_ratio * impeller_dia[sec_blade[0] - 1];

	float* airfoil_thicknessmax = new float[sec_blade[0]];
	float* airfoil_thickness = new float[SECTIONNO * sec_blade[0]];
	float* airfoil_length = new float[SECTIONNO * sec_blade[0]];
	airfoil_791(sec_blade, impeller_dia, head, l_chord, airfoil_thicknessmax, airfoil_thickness, airfoil_length);

	float* Control_pointx = new float[(SECTIONNO + EXTRAP) * sec_blade[0]];
	float* Control_pointy = new float[(SECTIONNO + EXTRAP) * sec_blade[0]];
	float* Control_pointz = new float[(SECTIONNO + EXTRAP) * sec_blade[0]];
	float* knotU = new float[SECTIONNO + EXTRAP + DEGREEU];
	float* knotV = new float[sec_blade[0] + DEGREEV];
	for (int i = 0; i < (SECTIONNO + EXTRAP) * sec_blade[0]; i++) {
		Control_pointx[i] = 0.0;
		Control_pointy[i] = 0.0;
		Control_pointz[i] = 0.0;
	}
	BsplineSurfInter(sec_blade, airfoil_length, airfoil_thickness, impeller_dia, Control_pointx, Control_pointy, Control_pointz, knotU, knotV);
	
	writefile(sec_blade, Control_pointx, Control_pointy, Control_pointz,  knotU, knotV, dvn);

	delete[] coeff;
	delete[] l_chord;
	delete[] blade_angle;
	delete[] impeller_dia;
	delete[] airfoil_thickness;
	delete[] airfoil_thicknessmax;
	return 0;
}
