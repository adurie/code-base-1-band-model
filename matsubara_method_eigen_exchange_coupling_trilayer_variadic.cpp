#include <iostream>
#include <utility>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
/* #include "cunningham.h" */
#include "cunningham_spawn.h"

using namespace Eigen;
using namespace std;

typedef complex<double> dcomp;

Matrix2cd S(double theta){
	Matrix2cd rotate;
	rotate << cos(theta/2.),sin(theta/2.),-sin(theta/2.),cos(theta/2.);
	return rotate;
}

VectorXcd Rspace(double x, double z, const double a, dcomp E, const int N) {
//...F(0)|NM(n)|F(theta)...
	const double t = 0.5;
	const double v1 = -1.6295;
	const double v2 = -2.8;
	const double v3 = -1.6295;
	const double nab = -0.6795;
	double F = cos(x*a)+cos(z*a);

	static const Matrix2cd T((Matrix2cd() << t,0.,0.,t).finished());
	static const Matrix2cd V2((Matrix2cd() << v2,0.,0.,v2).finished());
	static const Matrix2cd V3((Matrix2cd() << v3+nab,0.,0.,v3-nab).finished());
	static const Matrix2cd I((Matrix2cd() << 1.,0.,0.,1.).finished());
	static const Matrix2cd V1((Matrix2cd() << v1+nab,0.,0.,v1-nab).finished());
//initialise surface Green's function using mobius transformation for rotated layer at theta = 0
	Matrix2cd V3_0;
	V3_0 = S(0).inverse()*V3*S(0);
	Matrix2cd OMV3 = E*I-V3_0-2.*T*F;

	Matrix4cd X,O;
	X << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV3(0,0)/t,OMV3(0,1)/t,
		0,	-t,	OMV3(1,0)/t,OMV3(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	Matrix2cd b = O.topRightCorner(2,2);
	Matrix2cd d = O.bottomRightCorner(2,2);
	Matrix2cd GR_0 = b*d.inverse();

//initialise surface Green's function using mobius transformation for rotated layer at theta = PI 
	Matrix2cd V3_PI;
	V3_PI = S(M_PI).inverse()*V3*S(M_PI);
	Matrix2cd OMV2=E*I-V3_PI-2.*T*F;

	X << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV2(0,0)/t,OMV2(0,1)/t,
		0,	-t,	OMV2(1,0)/t,OMV2(1,1)/t;
	ces.compute(X);
	O = ces.eigenvectors();
	b = O.topRightCorner(2,2);
	d = O.bottomRightCorner(2,2);
	Matrix2cd GR_PI = b*d.inverse();

//initialise surface Green's function using mobius transformation for fixed layer 
	Matrix2cd OMV1=E*I-V1-2.*T*F;

	X << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV1(0,0)/t,OMV1(0,1)/t,
		0,	-t,	OMV1(1,0)/t,OMV1(1,1)/t;
	ces.compute(X);
	O = ces.eigenvectors();
	b = O.topRightCorner(2,2);
	d = O.bottomRightCorner(2,2);
	Matrix2cd GL = b*d.inverse();

	Matrix2cd Rsigma_0, Rsigma_PI;

	Matrix2cd OM = E*I - 2.*T*F;
	dcomp Fsigma;
	VectorXcd result(N);
	result.fill(0.);
//adlayer layer 2 from layer 1 to spacer thickness, N
	for (int it=0; it != N; ++it){

		GL = (OM - V2 -T*GL*T).inverse();
		Rsigma_0 = (I-GR_0*T.adjoint()*GL*T);
		Rsigma_PI = (I-GR_PI*T.adjoint()*GL*T);
		Fsigma = (1./M_PI)*log((Rsigma_0*Rsigma_PI.inverse()).determinant());
		result[it] = Fsigma;
	}
	
	return  result;
}

VectorXd f(const double a, const int N) {
	dcomp i;
	i = -1.;
	i = sqrt(i);
	dcomp E = 0.;
	const double Ef = 0.0;
	const double kT = 8.617342857e-5*300/13.6058;
	VectorXcd result_complex(N);
	result_complex.fill(0.);
	for (int j=0; j!=10; j++){
		E = Ef + (2.*j + 1.)*kT*M_PI*i;
		result_complex = result_complex + kspace(&Rspace, 2, 5e-2, a, E, N);
	}
	VectorXd result_return = result_complex.real();

	return kT*result_return;
}

int main() 
{

	/* cout<<"Name the data file\n"; */
	string Mydata = "test";
	/* getline(cin, Mydata); */
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );
	const double a = 1.;
	const int N = 50;

	/* dcomp i; */
	/* i = -1; */
	/* i = sqrt(i); */
	/* dcomp Ef = 0. + 1e-5*i; */	

	VectorXd result(N);
	//next line is gamma point only. Bypasses integration
	/* result = f(0,0,params); */
	result = f(a, N);
	result /= 4.*M_PI*M_PI;
	Myfile<<"N , Gamma"<<endl;

	for (int i=0; i < N ; ++i)
		Myfile << i+1 <<" ,  "<< -2.*M_PI*result[i] << endl;

	cout<<"finished!"<<endl;


return 0;
}
