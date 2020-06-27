#include <iostream>
#include <array>
#include "matrix.h"

using namespace std;

int main() {

	vector<vector<double>> vec= {
	{ 2, -1, 1 },
	{ 1, 3, -2 },
	{ 0, 1, -2 },
	};
	Matrix m1(vec);
	cout << "original matrix" << endl;
	m1.printMatrix();
	cout << endl;
	double det = m1.determinant();
	cout << "determinant: " << det << endl << endl;
	if (det != 0) {
		Matrix inv_m1 = m1.inverse();
		cout << "inverse:" << endl;
		inv_m1.printMatrix();
		cout << endl;
		Matrix identity = (inv_m1*m1);
		identity.trimMatrix();
		cout << "identity:" << endl;
		identity.printMatrix();
		cout << endl;
	}

	pair<Matrix, Matrix> p = m1.qr_decomposition();
	cout << "original matrix" << endl;
	m1.printMatrix();
	cout << endl;
	cout << "Q" << endl;
	p.first.printMatrix();
	cout << endl;
	cout << "R" << endl;
	p.second.printMatrix();
	cout << endl;

	Matrix s = p.first*p.second;
	s.trimMatrix();
	cout << "Q*R" << endl;
	s.printMatrix();
	if (s == m1) {
		cout << "QR passed the test." << endl;
	}
	else {
		cout << endl << "QR error: " << max((s - m1).apply(abs)) << endl;
	}
	cout << endl;

	array<Matrix, 3> arr = m1.singular_value_decomposition();
	Matrix & U = arr[0], & E = arr[1], & VT=arr[2];
	cout << "U" << endl;
	U.printMatrix();
	cout << endl;
	cout << "E" << endl;
	E.printMatrix();
	cout << endl;
	cout << "VT" << endl;
	VT.printMatrix();
	cout << endl;
	cout << "MATRIX" << endl;
	m1.printMatrix();
	cout << endl;
	U.trimMatrix();
	E.trimMatrix();
	VT.trimMatrix();
	Matrix res = U * E * VT;
	res.trimMatrix();
	cout << "U*E*VT" << endl;
	res.printMatrix();
	if (res == m1) {
		cout << "SVD passed the test." << endl;
	}
	else {
		cout << endl << "SVD error: " << max((res - m1).apply(abs)) << endl;
	}
	return 0;
}