#define EIGEN_DONT_PARALLELIZE
/* #define EIGEN_USE_MKL_ALL */
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "cunningham_points_integration.h"
#include <gsl/gsl_integration.h>

using namespace Eigen;
using namespace std;
typedef complex<double> dcomp;

struct a_struct{
	const double a = 1.;
	const double t = 0.5;
	const double Ef = -0.1;
	const double kT = 8.617342857e-5*400/13.6058;
	int N;
};

Matrix2cd S(double theta){
	Matrix2cd rotate;
	rotate << cos(theta/2.),sin(theta/2.),-sin(theta/2.),cos(theta/2.);
	return rotate;
}

Matrix2cd greens(Matrix2cd &OM, a_struct &params)
{
	Matrix4cd X,O;
	X << 	0,	0,	1/params.t,	0,
		0,	0,	0,	1/params.t,
		-params.t,	0,	OM(0,0)/params.t,OM(0,1)/params.t,
		0,	-params.t,	OM(1,0)/params.t,OM(1,1)/params.t;
	ComplexEigenSolver<Matrix4cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	Matrix2cd b = O.topRightCorner(2,2);
	Matrix2cd d = O.bottomRightCorner(2,2);
	Matrix2cd GR = b*d.inverse();
	return GR;
}

double Rspace(double x, double z, a_struct &params, dcomp E) {
//...F(0)|NM(n)|F(theta)...
	const double v1 = -1.5;
	const double v2 = -3.0238795325;
	const double v3 = -1.5;
	const double nab = 0.75;
	double F = cos(x*params.a)+cos(z*params.a);

	static const Matrix2cd T((Matrix2cd() << params.t,0.,0.,params.t).finished());
	static const Matrix2cd V2((Matrix2cd() <<v2,0.,0.,v2).finished());
	static const Matrix2cd V3((Matrix2cd() << v3+nab,0.,0.,v3-nab).finished());
	static const Matrix2cd I((Matrix2cd() << 1.,0.,0.,1.).finished());
	static const Matrix2cd V1((Matrix2cd() << v1+nab,0.,0.,v1-nab).finished());
//initialise surface Green's function using mobius transformation for rotated layer at theta = 0
	Matrix2cd V3_0;
	V3_0 = S(0).inverse()*V3*S(0);
	Matrix2cd OMV3=E*I-V3_0-2.*T*F;
	Matrix2cd GR_0;
	GR_0 = greens(OMV3, params);

//initialise surface Green's function using mobius transformation for rotated layer at theta = PI 
	Matrix2cd V3_PI;
	V3_PI = S(M_PI).inverse()*V3*S(M_PI);
	Matrix2cd OMV2=E*I-V3_PI-2.*T*F;
	Matrix2cd GR_PI;
	GR_PI = greens(OMV2, params);

//initialise surface Green's function using mobius transformation for fixed layer 
	Matrix2cd OMV1=E*I-V1-2.*T*F;
	Matrix2cd GL;
	GL = greens(OMV1, params);

	Matrix2cd OM = E*I - V2 - 2.*T*F;

	/* Matrix2cd root_arg = OM*T.inverse()*OM*T.inverse() - 4.*I;   //Mobius Transformation Method */
	/* Matrix2cd root; */
	/* root <<	sqrt(root_arg(0,0)),	0, */
	/* 		0,	sqrt(root_arg(1,1)); */
	/* Matrix2cd X_p = 0.5*(OM*T.inverse() + root); */
	/* Matrix2cd X_m = 0.5*(OM*T.inverse() - root); */
	/* Matrix2cd delta_Nm1, delta_Np1, delta_N, delta_1; */
	/* delta_Nm1 = X_p.pow(N - 1) - X_m.pow(N - 1); */
	/* delta_Np1 = X_p.pow(N + 1) - X_m.pow(N + 1); */
	/* delta_N = X_p.pow(N) - X_m.pow(N); */
	/* delta_1 = X_p - X_m; */
	/* Matrix2cd GN = T.inverse()*delta_1.inverse()*(delta_N - delta_Nm1*T*GL)*(delta_Np1 - delta_N*T*GL).inverse()*delta_1; */

	Matrix2cd Rsigma_0, Rsigma_PI;

	for (int it=0; it != params.N; ++it)			//Adlayer Method
		GL = (OM - T*GL*T).inverse();

	double Fsigma;
	Rsigma_0 = (I-GR_0*T.adjoint()*GL*T);
	Rsigma_PI = (I-GR_PI*T.adjoint()*GL*T);
	Fsigma = imag((1./M_PI)*log((Rsigma_0*Rsigma_PI.inverse()).determinant()));
	
	double fermi = 1./(1.+exp((real(E)-params.Ef)/params.kT));
	double integrand = Fsigma*fermi;
	/* cout<<"E = "<<E<<endl; */
	/* cout<<"integrand = "<<integrand<<endl; */
	return integrand;
}

double temp(double E, void * alpha) {
	a_struct params = *(a_struct *) alpha;
	dcomp i, im, Ec;
	i = -1;
	i = sqrt(i);
	im = 1e-5*i;	
	Ec = E + im;
	double result = kspace(params, &Rspace, 0, 0, Ec);
	return result;
}

int main() 
{

	a_struct params;
	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	int max_N = 50;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;

	gsl_function F;
	F.function = &temp;
	F.params = &params;
	const double start = -6.;
	const double end = params.Ef + .1;

	/* dcomp i; */
	/* i = -1; */
	/* i = sqrt(i); */
	/* dcomp Ef = 0. + 1e-5*i; */	

	//next line is gamma point only. Bypasses integration
	/* result = f(0,0,params); */
	/* result *= 2.*kT*M_PI; */
	/* Myfile<<"N , Gamma"<<endl; */

	for (params.N = 0; params.N<=max_N; params.N++ ){
		gsl_integration_qags(&F, start, end, 0, 3e-2, 1000, w, &result, &error);
		cout<<error<<endl;
		Myfile.open( Mydata.c_str(),ios::app );

		Myfile<< params.N <<" , "<< result << endl;
		cout<<100*params.N/(max_N*1.)<<"\% completed"<<endl;
		Myfile.close();
	}

	cout<<"finished!"<<endl;


return 0;
}
