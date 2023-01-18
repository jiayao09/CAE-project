#include "Header.h"
#include <cmath>
#include <iostream>


void cal_sectionblade(const int n_s, int* output_sec) {
	int section = 0;
	int blade = 0;
	// Determine number of the section
	if (n_s < 450) section = 4;
	else if (n_s >= 450 && n_s <= 800) section = 5;
	else if (n_s > 800) section = 6;
	else section = 0;
	// Determine number of blade
	if (n_s >= 0 && n_s <= 600) blade = 5;
	else if (n_s > 600 && n_s <= 850) blade = 4;
	else if (n_s > 850 && n_s <= 1500) section = 3;
	else section = 0;

	output_sec[0] = section;
	output_sec[1] = blade;
	// std::cout << "section = " << section << "blade = " << blade << std::endl;
}


void cal_coeffab(const int* input, const float capacity, const float head, const int n, const int n_s, int l_axial, float* output_coffab)
{
	float u = 0.00843 * pow(capacity, 0.5) * PI * n;
	float v_rotate = 12.9 * head / u;
	float axial_length;
	float length;
	float error;
	float v_m;
	float beta_1;
	float beta_2;
	float angle;
	float coefficientamax[3][3] = { {5.5, 5.8, 6.3},{7.3, 7.6, 8.9},{8.9, 9.4, 8.9} };
	float coefficientamin[3][3] = { {2.8, 3.4, 4.2},{5.5, 5.8, 6.3},{7.3, 7.6, 10.2} };
	float coefficientbmax[4][3] = { {0.22, 0.24, 0.28},{0.16, 0.19, 0.21},{0.13, 0.13, 0.16},{0.07,0.08,0.12} };
	float coefficientbmin[4][3] = { {0.16, 0.19, 0.21},{0.13, 0.13, 0.16},{0.07, 0.08, 0.12},{0.03,0.03,0.05} };
	int section_b = 0;
	if (n_s < 380) section_b = 0;
	else if (n_s >= 380 && n_s < 610) section_b = 1;
	else if (n_s >= 610 && n_s < 930) section_b = 2;
	else section_b = 3;
	float range_a = coefficientamax[input[1] - 3][input[0] - 4] - coefficientamin[input[1] - 3][input[0] - 4];
	float range_b = coefficientbmax[section_b][input[0] - 4] - coefficientbmin[section_b][input[0] - 4];
	float factor_a = 0;
	float factor_b = 0;
	float coff_a = 0;
	float coff_b = 0;
	int count_a = range_a / FACTORA + 2;
	int count_b = range_b / FACTORB + 2;
	float* output_error = new float[count_a * count_b];
	float* output_j = new float[count_a * count_b];
	for (int i = 0; i < count_a; i++) {
		factor_a = (FACTORA * i) / range_a;
		coff_a = (range_a)*factor_a + coefficientamin[input[1] - 3][input[0] - 4];
		length = (89.3 * pow(capacity, 0.5) * coff_a) / input[1];
		for (int j = 0; j < count_b; j++) {
			factor_b = (FACTORB * j) / range_b;
			coff_b = (range_b)*factor_b + coefficientbmin[section_b][input[0] - 4];
			v_m = 4 / (PI * (0.276 - coff_b));
			beta_1 = atan(v_m / u);
			beta_2 = atan(v_m / (u - v_rotate));
			angle = (beta_1 + beta_2) / 2;
			axial_length = length * sin(angle);
			error = (l_axial / axial_length) * 100;
			output_error[i * count_b + j] = 100 - error;
			output_j[i * count_b + j] = j;
			if (error >= 100) {
				output_error[i * count_b + j] = 100;
			}
		}
	}

	int index = 0;
	float min = output_error[0];
#pragma omp parallel
	{
		int size = count_a * count_b;
		int index_local = index;
		float min_local = min;
#pragma omp parallel for nowait
		for (int i = 1; i < size; i++) {
			if (output_error[i] < min_local) {
				min_local = output_error[i];
				index_local = i;
			}
		}
#pragma omp critical 
		{
			if (min_local < min) {
				min = min_local;
				index = index_local;
			}
		}
	}

	int indexi = (index - output_j[index]) / count_b;
	int indexj = output_j[index];
	factor_a = (FACTORA * indexi) / range_a;
	coff_a = (range_a)*factor_a + coefficientamin[input[1] - 3][input[0] - 4];
	factor_b = (FACTORB * indexj) / range_b;
	coff_b = (range_b)*factor_b + coefficientbmin[section_b][input[0] - 4];
	output_coffab[0] = coff_a;
	output_coffab[1] = coff_b;

	delete[] output_error;
	delete[] output_j;
}


void bladeangle(const int* input, const float* coeff, const int n, const float capacity, const float head, float* output_angle)
{
	float coefficientmax[3][6] = { {2.24, 1.73, 1.56, 1.0, 0.0, 0.0},{2.18, 1.68, 1.34, 1.15, 1.0, 0.0},{2.06, 1.53, 1.42, 1.21, 0.92, 1.0} };
	float coefficientmin[3][6] = { {1.92, 1.52, 1.36, 1.0, 0.0, 0.0},{1.84, 1.43, 1.22, 1.06, 1.0, 0.0},{1.72, 1.21, 1.17, 0.97, 0.93, 1.0} };
	int section = input[0] - 4;
	float u = 0.00843 * pow(capacity, 0.5) * PI * n;
	float v_rotate = 12.9 * head / u;
	float v_m = 4 / (PI * (0.276 - coeff[1]));
	float beta_1 = atan(v_m / u);
	float beta_2 = atan(v_m / (u - v_rotate));
	float angle = (beta_1 + beta_2) / 2;
	float coff_b1 = 0;
	// std::cout << "angle = " << angle << std::endl;
	for (int i = 0; i < input[0]; i++) {
		float range = coefficientmax[section][i] - coefficientmin[section][i];
		coff_b1 = coefficientmin[section][i] + (range) * PROPCOEFFA;
		output_angle[i] = angle * coff_b1;
	}
}


void cal_chordlength(const int* input, const float* coeff, const float capacity, float* output_chord) {
	float coefficientmax[3][6] = { {0.728, 0.873,0.981,1.0, 0.0, 0.0},{0.685, 0.787, 0.894, 0.963, 1.0, 0.0},{0.553, 0.653, 0.781, 0.842, 0.925, 1.0} };
	float coefficientmin[3][6] = { {0.651, 0.793,0.894,1.0, 0.0, 0.0},{0.623, 0.712, 0.826, 0.931, 1.0, 0.0},{0.489, 0.586, 0.705, 0.793, 0.856, 1.0} };
	int section = input[0] - 4;
	float length = (89.3 * pow(capacity, 0.5) * coeff[0]) / input[1];
	float range = 0;
	float coff_a1 = 0;
	// std::cout << "length = " << length << std::endl;
	for (int i = 0; i < input[0]; i++) {
		range = coefficientmax[section][i] - coefficientmin[section][i];
		coff_a1 = coefficientmin[section][i] + range * PROPCOEFFA;
		output_chord[i] = length * coff_a1;

	}
}


void cal_impellerDia(const int* input, const float* coeff, const float* l_chord, const float capacity, const int n_s, float* impeller_dia) {
	float coefficientc[3][6] = { {0.5, 0.64, 0.76, 0.98, 0.0, 0.0},{0.5, 0.61, 0.73, 0.85, 0.97, 0.0},{0.5, 0.53, 0.57, 0.69, 0.82, 0.93} };
	float hub_ratio = -0.00029 * n_s + 0.71587;
	float lg = 0.02377;
	int coefficientK[3] = { 21, 19, 17 };  // section 4 (19.3-22.45), section 5 (17.8-20.14), section (15.8-19.6)
	int section = input[0] - 4;
	for (int i = 0; i < input[0]; i++) {
		float l = l_chord[i];
		impeller_dia[i] = coefficientc[section][i] * (10.5 + lg * l) / coefficientK[section] * pow(capacity, 0.5);
	}
}

void airfoil_791(const int* input, float* impeller_dia, const float head, const float* l_chord, float* airfoil_thicknessmax, float* airfoil_thickness, float* airfoil_length) {
	float thickness_ratio[SECTIONNO] = { 0.0, 0.296, 0.405, 0.489,0.778,0.92, 0.978, 1.0, 0.883, 0.756, 0.544, 0.356, 0.2, 0.0 };
	float length_ratio[SECTIONNO] = { 0.0, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 };
	for (int i = 0; i < input[0]; i++) {
		airfoil_thicknessmax[i] = impeller_dia[i] * pow(head, 0.5) * 0.014 * 1000;
		for (int j = 0; j < SECTIONNO; j++) {
			airfoil_length[i * SECTIONNO + j] = length_ratio[j] * l_chord[i];
			airfoil_thickness[i * SECTIONNO + j] = thickness_ratio[j] * airfoil_thicknessmax[i];
		}
	}
}
