#include "Header.h"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <fstream>

void linearslove(const float* A, const float* B, float* C, int x, int y, int order) {
	MatrixXd A1(x, y);
	VectorXd B1(y);
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			A1(i, j) = A[i * y + j];
			B1(j) = B[j];
		}
	
	}

	// std::cout << "Here is the matrix A:\n" << A1 << std::endl;
	// std::cout << "Here is the vector b:\n" << B1 << std::endl;
	Eigen::VectorXd sol = A1.colPivHouseholderQr().solve(B1);
	// std::cout << "The solution is:\n" << sol << std::endl;

	for (int i = 0; i < y; i++) {
		C[i] = sol(i);
	}
}


void disp(const float* array, int size, int x) {
	int y = size / x;

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			cout << array[i * x + j] << " , ";
		}
		cout << endl;
	}
}



void disp(const float* array, int size) {

	for (int i = 0; i < size; i++) {
		cout << array[i] << endl;
	}
	cout << endl;
}


void writefile(const int* input, float* Control_pointx, float* Control_pointy, float* Control_pointz, float* knotU, float* knotV, int dvn) {
	ofstream MyFile("ControlPoint_info.txt");
	if (MyFile.is_open())
	{
		// Write STL info
		MyFile << "N = " << dvn << endl;
		MyFile << "ku = " << DEGREEU << endl;
		MyFile << "nu = " << SECTIONNO + EXTRAP << endl;
		MyFile << "kv = " << DEGREEV << endl;
		MyFile << "nv = " << input[0] << endl;

		// Write control Point info
		MyFile << "Control_pointx" << endl;
		for (int i = 0; i < input[0]; i++) {
			for (int j = 0; j < (SECTIONNO + EXTRAP); j++) {
				MyFile << Control_pointx[i * (SECTIONNO + EXTRAP) + j] << ",";
			}
			MyFile << endl;
		}

		MyFile << "Control_pointy" << endl;
		for (int i = 0; i < input[0]; i++) {
			for (int j = 0; j < (SECTIONNO + EXTRAP); j++) {
				MyFile << Control_pointy[i * (SECTIONNO + EXTRAP) + j] << ",";
			}
			MyFile << endl;
		}

		MyFile << "Control_pointz" << endl;
		for (int i = 0; i < input[0]; i++) {
			for (int j = 0; j < (SECTIONNO + EXTRAP); j++) {
				MyFile << Control_pointz[i * (SECTIONNO + EXTRAP) + j] << ",";
			}
			MyFile << endl;
		}

		// Write U info
		MyFile << "knot U" << endl;
		for (int i = 0; i < SECTIONNO + EXTRAP + DEGREEU; i++) {
			MyFile << knotU[i] << ",";
		}
		MyFile << endl;
		// Write V info
		MyFile << "knot V" << endl;
		for (int i = 0; i < input[0] + DEGREEV; i++) {
			MyFile << knotV[i] << ",";
		}
		MyFile << endl;
		// Close the file
	}
	MyFile.close();
}
