//#define FFTWPP_SINGLE_THREAD
#include "Array.h"
#include "fftw++-2.03/fftw++.h"
#include <iostream>
#include <utility>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

// Compile with:
//clang++-3.8 -I/usr/include/fftw++-2.03 -fopenmp=libiomp5 J_SPA_FFT.cpp /usr/include/fftw++-2.03/fftw++.cc -lfftw3 -lfftw3_omp -std=c++11
// FFT to solve J(N) vs N using SPA - requires inverse transform, taken complex conjugate of FFT, reversed sign of imag(FFT) and reverse sign of output. Also divide by total number of FFT points taken.
using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;
using namespace Eigen;
typedef complex<double> dcomp;	

struct a_struct{
	const double a = 1.;
	const double v2 = -2.8;
	const double kT = 8.617342857e-5*300/13.6058;
	dcomp E;
	const double x = 0.;
	const double z = 0.;
	double P;	//period of J against spacer thickness
	const double t = 0.5;
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

dcomp Rspace(double N, a_struct &params) {
//...F(0)|NM(n)|F(theta)...
	const double v1 = -1.6295;
	const double v3 = -1.6295;
	const double nab = -0.6795;
	double F = cos(params.x*params.a)+cos(params.z*params.a);

	static const Matrix2cd T((Matrix2cd() << params.t,0.,0.,params.t).finished());
	static const Matrix2cd V2((Matrix2cd() << params.v2,0.,0.,params.v2).finished());
	static const Matrix2cd V3((Matrix2cd() << v3+nab,0.,0.,v3-nab).finished());
	static const Matrix2cd I((Matrix2cd() << 1.,0.,0.,1.).finished());
	static const Matrix2cd V1((Matrix2cd() << v1+nab,0.,0.,v1-nab).finished());
//initialise surface Green's function using mobius transformation for rotated layer at theta = 0
	Matrix2cd V3_0;
	V3_0 = S(0).inverse()*V3*S(0);
	Matrix2cd OMV3=params.E*I-V3_0-2.*T*F;
	Matrix2cd GR_0;
	GR_0 = greens(OMV3, params);

//initialise surface Green's function using mobius transformation for rotated layer at theta = PI 
	Matrix2cd V3_PI;
	V3_PI = S(M_PI).inverse()*V3*S(M_PI);
	Matrix2cd OMV2=params.E*I-V3_PI-2.*T*F;
	Matrix2cd GR_PI;
	GR_PI = greens(OMV2, params);

//initialise surface Green's function using mobius transformation for fixed layer 
	Matrix2cd OMV1=params.E*I-V1-2.*T*F;
	Matrix2cd GL;
	GL = greens(OMV1, params);

	Matrix2cd OM = params.E*I - V2 - 2.*T*F;
	Matrix2cd root_arg = OM*T.inverse()*OM*T.inverse() - 4.*I;
	Matrix2cd root;
	root <<	sqrt(root_arg(0,0)),	0,
			0,	sqrt(root_arg(1,1));
	Matrix2cd X_p = 0.5*(OM*T.inverse() + root);
	Matrix2cd X_m = 0.5*(OM*T.inverse() - root);
	Matrix2cd delta_Nm1, delta_Np1, delta_N, delta_1;
	delta_Nm1 = X_p.pow(N - 1) - X_m.pow(N - 1);
	delta_Np1 = X_p.pow(N + 1) - X_m.pow(N + 1);
	delta_N = X_p.pow(N) - X_m.pow(N);
	delta_1 = X_p - X_m;

	Matrix2cd GN = T.inverse()*delta_1.inverse()*(delta_N - delta_Nm1*T*GL)*(delta_Np1 - delta_N*T*GL).inverse()*delta_1;

	Matrix2cd Rsigma_0, Rsigma_PI;
	dcomp Fsigma;
	Rsigma_0 = (I-GR_0*T.adjoint()*GN*T);
	Rsigma_PI = (I-GR_PI*T.adjoint()*GN*T);
	Fsigma = ((1./M_PI)*log((Rsigma_0*Rsigma_PI.inverse()).determinant()));
	
	return Fsigma;
}

Matrix2d hessian(a_struct &params)
{
	double Cz = cos(params.a*params.z);
	double Sz = sin(params.a*params.z);
	double Cx = cos(params.a*params.x);
	double Sx = sin(params.a*params.x);
	double disp = (params.v2-real(params.E))/(2.*params.t) + Cx + Cz;
	double el_11 = -params.a*Cx/sqrt(-disp*disp +1.) + params.a*Sx*Sx*disp/pow((-disp*disp + 1.),1.5);
	double el_22 = -params.a*Cz/sqrt(-disp*disp +1.) + params.a*Sz*Sz*disp/pow((-disp*disp + 1.),1.5);
	double el_12 = params.a*Sx*Sz*disp/pow((-disp*disp +1.),1.5);
	Matrix2d Hes;
	Hes << el_11,	el_12,
	    	el_12,	el_22;
	/* cout<<el_11<<'\t'<<el_12<<'\n'<<el_12<<'\t'<<el_22<<'\n'<<endl; */
	return Hes;
}

double Psi(dcomp func){
	dcomp Func = conj(func);
	/* double psi = atan2(-imag(Func),real(Func)); */
	double psi = abs(atan(-imag(Func)/real(Func)));
	if ((imag(Func) > 0) && (real(Func) < 0))
		psi = M_PI - psi;
	if ((imag(Func) < 0) && (real(Func) > 0))
		psi = -psi;
	if ((imag(Func) < 0) && (real(Func) < 0))
		psi = -M_PI + psi;
	return psi;
}
  
int main()
{
  //fftw::maxthreads=get_max_threads();
        a_struct params;
	dcomp i;
	i = -1.;
	i = sqrt(i);
	params.E = 0.0 + i*1e-10;

	double F = cos(params.x*params.a)+cos(params.z*params.a);
	double y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	double dy_dE = -1./(2.*params.t*sin(y));
	Matrix2d Hes = hessian(params);
	double det_hes = sqrt(abs(Hes(0,0)*Hes(1,1)-Hes(0,1)*Hes(1,0)));
	/* cout<<det_hes<<endl; */
	dcomp tau;
	/* cout<<Hes.trace()<<endl; */
	SelfAdjointEigenSolver<Matrix2d> es;
	es.compute(Hes);
	VectorXd eigenvals = es.eigenvalues();

	if ((eigenvals(0) < 0) && (eigenvals(1) < 0))
		tau = -i;
	if ((eigenvals(0) > 0) && (eigenvals(1) > 0))
		tau = i;
	if (((eigenvals(0) > 0) && (eigenvals(1) < 0)) || ((eigenvals(0) < 0) && (eigenvals(1) > 0)))
		tau = 1.;

	unsigned int n=128; 
	size_t align = sizeof(Complex);
	array1<Complex> f(n,align);
	array1<Complex> fu(n,align);
	array1<Complex> fd(n,align);
	array1<Complex> f2u(n,align);
	array1<Complex> f2d(n,align);

	fft1d Forward(n,-1);
	fft1d Backward(n,1);

	for(unsigned int k=0; k < n; k++){
		y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
		params.P = M_PI/y;
		f[k] = Rspace(k*(params.P)/(n-1), params);

		const double delta = 1e-2;
		params.E -= delta;
		y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
		params.P = M_PI/y;
		fd[k] = Rspace(k*(params.P)/(n-1), params);
		params.E += delta;

		params.E += delta;
		y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
		params.P = M_PI/y;
		fu[k] = Rspace(k*(params.P)/(n-1), params);
		params.E -= delta;	
	
		const double delta_check = 1.2;
		params.E -= delta_check*delta;
		y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
		params.P = M_PI/y;
		f2d[k] = Rspace(k*(params.P)/(n-1), params);
		params.E += delta_check*delta;

		params.E += delta_check*delta;
		y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
		params.P = M_PI/y;
		f2u[k] = Rspace(k*(params.P)/(n-1), params);
		params.E -= delta_check*delta;
	}

	Backward.fft(f);
	Backward.fft(fu);
	Backward.fft(fd);
	Backward.fft(f2u);
	Backward.fft(f2d);

	VectorXd dPsi_dE;
	dPsi_dE.resize(n);
	double psi, psi_u, psi_2u, psi_d, psi_2d, dsi_dE;

	for(unsigned int k=0; k < n; k++){
		const double delta = 1e-2;
		const double delta_check = 1.2;
		psi = Psi(f[k]*(static_cast<double>(n)));
		psi_u = Psi(fu[k]*(static_cast<double>(n)));
		psi_2u = Psi(f2u[k]*(static_cast<double>(n)));
		psi_d = Psi(fd[k]*(static_cast<double>(n)));
		psi_2d = Psi(f2d[k]*(static_cast<double>(n)));
		dsi_dE = (psi_u - psi_d)/(2*delta);
		double test = (psi_2u - psi_2d)/(2*delta_check*delta);
		double abs_err = abs(abs(dsi_dE/test)-1);
		dPsi_dE[k] = dsi_dE;
		/* if ( abs_err > 0.05) */
			/* cout<<"caution: for s = "<<k<<", dpsi/dE has absolute error of "<<abs_err<<endl; */
		/* if ( abs(abs(2.*psi/(psi_u + psi_d))-1) > 0.05) */
			/* cout<<"caution: for s = "<<k<<", psi is not on the gradient line"<<endl; */
	}

	string Mydata = "out.txt";
	ofstream Myfile;
	Myfile.open(Mydata.c_str(), ios::trunc);
	Myfile<<"N , J(N)"<<endl;
	double spincurrent; 
	y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	params.P = M_PI/y;

	for (double j = 1; j<=51; j = j + 0.02){
		spincurrent = 0.;

		for(unsigned int s=1; s <n; s++){

			spincurrent += real(tau*conj(f[s])*exp(i*static_cast<double>(j*s)*2.*y)
				/(2*s*det_hes*sinh(M_PI*params.kT*(2.*j*s*dy_dE 
					/* )))); */	
					+ dPsi_dE[s]))));

		}

	Myfile<< j <<" , "<< (1/static_cast<double>(n))*(params.kT/(2.*j))*spincurrent << endl;
	}
	return 0;
}
