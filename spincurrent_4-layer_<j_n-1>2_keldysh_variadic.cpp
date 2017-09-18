#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
/* #include "cunningham_points_adaptive.h" */
#include "cunningham.h"

using namespace Eigen;
using namespace std;

typedef complex<double> dcomp;
struct a_struct{
	const int N = 50;
	const double a = 1.;
};

VectorXd f(double x, double z, const double a, const int N) {
//F(theta)/NM(n)/F(0)/NM
	const double t = 0.5;
	const double v1 = -1.5;
	const double v2 = -2.9238795325;
	const double theta = M_PI/2.;
	const double nab = 0.4;
	dcomp i;
//	double e = 1.6e-19;
//	double h = 6.626e-34;
	double F = cos(x*a)+cos(z*a);
// v2=-2.9238795325 has integer periodicity in the spacer by dispersion relation
	const static Matrix2cd T = (Matrix2cd(2,2) << t,0,0,t).finished();
	const static Matrix2cd V2 = (Matrix2cd(2,2) << v2,0.,0.,v2).finished();
	const static Matrix2cd V3 = (Matrix2cd(2,2) << v1+nab,0.,0.,v1-nab).finished();
	const static Matrix2cd I = (Matrix2cd(2,2) << 1.,0.,0.,1.).finished();
	const static Matrix2cd S = (Matrix2cd(2,2) << cos(theta/2.),sin(theta/2.),-sin(theta/2.),cos(theta/2.)).finished();
	const static Matrix2cd V1 = S.inverse()*V3*S;
	i=-1;
	i=sqrt(i);
	const dcomp E=(1e-5)*i;
	Matrix2cd OMV2=E*I-V2-2.*T*F;

	Matrix4cd X2,O2;
	X2 << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV2(0,0)/t,OMV2(0,1)/t,
		0,	-t,	OMV2(1,0)/t,OMV2(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces2;
	ces2.compute(X2);
	O2 = ces2.eigenvectors();
	Matrix2cd b2 = O2.topRightCorner(2,2);
	Matrix2cd d2 = O2.bottomRightCorner(2,2);
	Matrix2cd GR = b2*d2.inverse();

	Matrix2cd OM = E*I - 2.*T*F;

	Matrix2cd OMV1=E*I-V1-2.*T*F;

	Matrix4cd X,O,Oinv;
	X << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV1(0,0)/t,OMV1(0,1)/t,
		0,	-t,	OMV1(1,0)/t,OMV1(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	Matrix2cd b = O.topRightCorner(2,2);
	Matrix2cd d = O.bottomRightCorner(2,2);
	Matrix2cd GL = b*d.inverse();

	Matrix2cd Pauli, xPau, yPau, zPau;
	xPau << 0.,1.,1.,0.;
	yPau << 0.,-i,i,0.;
	zPau << 1.,0.,0.,-1.;

	Pauli = yPau;

	VectorXd spincurrent(N);
	Matrix2cd A,B,TOT,GM;
	GM = (OM - V3 -T*GR*T).inverse();
//lim is thickness of layer 3
	const int lim = 9;
//build thickness of layer 3 to lim layers
	for (int it=0; it < lim; ++it){

		GM = (OM - V3 -T*GM*T).inverse();
	}
//adlayer layer 2 from layer 1 to spacer thickness, N
	for (int it=0; it < N; ++it){

		GL = (OM - V2 -T*GL*T).inverse();
		A = (I-GM*T.adjoint()*GL*T).inverse();
		B = (I-GM.adjoint()*T.adjoint()*GL.adjoint()*T).inverse();
		TOT = (GL*T*A*B*GM.adjoint()*T.adjoint() - A*B + 0.5*(A+B))*Pauli;
		spincurrent[it] = (1./(8.*M_PI*M_PI*M_PI))*real(TOT(0,0)+TOT(1,1));//remember conversion of summation to integral (a/2pi)^2
	}
//	spincurrent[]
//	jy = (1./(4.*M_PI))*real(TOT(0,0)+TOT(1,1))
	return spincurrent;
	
}

int main() {
	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );
	const int N = 50;
	const double a = 1.;
	VectorXd integral;
// only include the next term for gamma point only calculations
//	integral = integral + f(0.,0.,N,a);
	/* integral = kspace(params, &f, 4096, 1e-5); */
	integral = kspace(&f, 3, 1e-2, 768, a, N);
	Myfile<<"N , Gamma"<<endl;
	for (int J=0;J<abs(N);++J) {

		Myfile<<J+1<<" , "<<integral[J]<<endl;
	}

	cout<<"finished!"<<endl;


return 0;
}

