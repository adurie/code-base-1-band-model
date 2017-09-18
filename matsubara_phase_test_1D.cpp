#include <iostream>
#include <utility>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
/* #include "cunningham_points_adaptive.h" */

using namespace Eigen;
using namespace std;

typedef complex<double> dcomp;

/* Matrix2cd S(double theta){ */
/* 	Matrix2cd rotate; */
/* 	rotate << cos(theta/2.),sin(theta/2.),-sin(theta/2.),cos(theta/2.); */
/* 	return rotate; */
/* } */

dcomp greens(dcomp OM, double t)
{
	Matrix2cd X,O;
	X << 	0,	1/t,
		-t,	OM/t;

	ComplexEigenSolver<Matrix2cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	dcomp b = O(0,1);
	dcomp d = O(1,1);
	dcomp GR = b/d;
	return GR;
}

VectorXcd Rspace(const double a, dcomp E, const int N) {
//...F(0)|NM(n)|F(theta)...
	const double t_1 = 0.15;
	const double t_2 = 0.4;
	const double t_3 = 0.15;
	const double v1 = -0.6;
	const double v2 = -0.79;
	const double v3 = -0.6;
	const double nab = -0.1;
	double F = 0.;//cos(x*a)+cos(z*a);

	/* static const Matrix2cd T_1((Matrix2cd() << t_1,0.,0.,t_1).finished()); */
	/* static const Matrix2cd T_2((Matrix2cd() << t_2,0.,0.,t_2).finished()); */
	/* static const Matrix2cd T_3((Matrix2cd() << t_3,0.,0.,t_3).finished()); */
	/* static const Matrix2cd V2((Matrix2cd() << v2,0.,0.,v2).finished()); */
	/* static const Matrix2cd V3((Matrix2cd() << v3+nab,0.,0.,v3-nab).finished()); */
	/* static const Matrix2cd I((Matrix2cd() << 1.,0.,0.,1.).finished()); */
	/* static const Matrix2cd V1((Matrix2cd() << v1+nab,0.,0.,v1-nab).finished()); */
	
//initialise surface Green's function using mobius transformation for rotated layer at theta = 0
	double V3_0_u = v3+nab;
	dcomp OMV3_u= E-V3_0_u-2.*t_3*F;
	dcomp GR_0_u = greens(OMV3_u, t_3);

	double V3_0_d = v3-nab;
	dcomp OMV3_d= E-V3_0_d-2.*t_3*F;
	dcomp GR_0_d = greens(OMV3_d, t_3);

//initialise surface Green's function using mobius transformation for rotated layer at theta = PI 
	double V3_PI_u = v3-nab;
	dcomp OMV2_u= E-V3_PI_u-2.*t_3*F;
	dcomp GR_PI_u = greens(OMV2_u, t_3);

	double V3_PI_d = v3+nab;
	dcomp OMV2_d= E-V3_PI_d-2.*t_3*F;
	dcomp GR_PI_d = greens(OMV2_d, t_3);

//initialise surface Green's function using mobius transformation for fixed layer 
	double V1_u = v1+nab;
	dcomp OMV1_u= E-V1_u-2.*t_1*F;
	dcomp GL_u = greens(OMV1_u, t_1);

	double V1_d = v1-nab;
	dcomp OMV1_d= E-V1_d-2.*t_1*F;
	dcomp GL_d = greens(OMV1_d, t_1);

	dcomp Rsigma_0_u, Rsigma_0_d, Rsigma_PI_u, Rsigma_PI_d;
	dcomp Fsigma;

	dcomp OM = E - 2.*t_2*F;
	VectorXcd result(N);
	result.fill(0.);
//adlayer layer 2 from layer 1 to spacer thickness, N
	for (int it=0; it != N; ++it){

		GL_u = 1./(OM - v2 -t_2*GL_u*t_2);
		GL_d = 1./(OM - v2 -t_2*GL_d*t_2);
		Rsigma_0_u = (1.-GR_0_u*t_2*GL_u*t_2);
		Rsigma_0_d = (1.-GR_0_d*t_2*GL_d*t_2);
		Rsigma_PI_u = (1.-GR_PI_u*t_2*GL_u*t_2);
		Rsigma_PI_d = (1.-GR_PI_d*t_2*GL_d*t_2);
		Fsigma = (1./M_PI)*log(Rsigma_0_d);//*Rsigma_0_u/(Rsigma_PI_u*Rsigma_PI_d));
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
		result_complex = result_complex + Rspace(a, E, N);
	}
	VectorXd result_return = result_complex.real();

	return kT*result_return;
}

int main() 
{

	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
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
	/* result /= 4.*M_PI*M_PI; */
	Myfile<<"N , Gamma"<<endl;

	for (int i=0; i < N ; ++i)
		Myfile << i+1 <<" ,  "<< -2.*M_PI*result[i] << endl;

	cout<<"finished!"<<endl;


return 0;
}
