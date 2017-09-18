#include <iostream>
#include <utility>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace std;

typedef complex<double> dcomp;

struct a_struct{
	const double a = 1.;
	const double v2 = -.79;
	const double kT = 8.617342857e-5*300/13.6058;
	dcomp E;
	/* const double x = 0.; */
	/* const double z = 0.; */
	int s;
	double P;	//period of J against spacer thickness
	const double t = 0.4;
};

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

dcomp Rspace(double N, a_struct &params) {
//...F(0)|NM(n)|F(theta)...
	const double v1 = -0.6;
	/* const double v1 = -2.5; */
	const double v3 = -0.6;
	/* const double v3 = -2.5; */
	const double nab = -0.1;
	double F = 0;//cos(params.x*params.a)+cos(params.z*params.a);

	const double t_1 = 0.15;
	/* static const Matrix2cd T_1((Matrix2cd() <<t_1,0.,0.,t_1).finished()); */
	/* static const Matrix2cd T_2((Matrix2cd() <<params.t,0.,0.,params.t).finished()); */
	const double t_3 = 0.15;
	/* static const Matrix2cd T_3((Matrix2cd() <<t_3,0.,0.,t_3).finished()); */
	/* static const Matrix2cd V2((Matrix2cd() << params.v2,0.,0.,params.v2).finished()); */
	/* static const Matrix2cd V3((Matrix2cd() << v3+nab,0.,0.,v3-nab).finished()); */
	/* static const Matrix2cd I((Matrix2cd() << 1.,0.,0.,1.).finished()); */
	/* static const Matrix2cd V1((Matrix2cd() << v1+nab,0.,0.,v1-nab).finished()); */
//initialise surface Green's function using mobius transformation for rotated layer at theta = 0
	double V3_0_u = v3+nab;
	dcomp OMV3_u=params.E-V3_0_u-2.*t_3*F;
	dcomp GR_0_u = greens(OMV3_u, t_3);

	double V3_0_d = v3-nab;
	dcomp OMV3_d=params.E-V3_0_d-2.*t_3*F;
	dcomp GR_0_d = greens(OMV3_d, t_3);

//initialise surface Green's function using mobius transformation for rotated layer at theta = PI 
	double V3_PI_u = v3-nab;
	dcomp OMV2_u=params.E-V3_PI_u-2.*t_3*F;
	dcomp GR_PI_u = greens(OMV2_u, t_3);

	double V3_PI_d = v3+nab;
	dcomp OMV2_d=params.E-V3_PI_d-2.*t_3*F;
	dcomp GR_PI_d = greens(OMV2_d, t_3);

//initialise surface Green's function using mobius transformation for fixed layer 
	double V1_u = v1+nab;
	dcomp OMV1_u=params.E-V1_u-2.*t_1*F;
	dcomp GL_u = greens(OMV1_u, t_1);

	double V1_d = v1-nab;
	dcomp OMV1_d=params.E-V1_d-2.*t_1*F;
	dcomp GL_d = greens(OMV1_d, t_1);

	dcomp OM = params.E - params.v2 - 2.*params.t*F;
	dcomp root = sqrt((OM/params.t)*(OM/params.t) - 4.);
	dcomp X_p = 0.5*(OM/params.t + root);
	dcomp X_m = 0.5*(OM/params.t - root);
	dcomp delta_Nm1, delta_Np1, delta_N, delta_1;
	delta_Nm1 = pow(X_p,(N - 1)) - pow(X_m,(N - 1));
	delta_Np1 = pow(X_p,(N + 1)) - pow(X_m,(N + 1));
	delta_N = pow(X_p,N) - pow(X_m,N);
	delta_1 = X_p - X_m;

	dcomp GN_u = (1./(params.t*delta_1))*((delta_N - delta_Nm1*params.t*GL_u)/(delta_Np1 - delta_N*params.t*GL_u))*delta_1;
	dcomp GN_d = (1./(params.t*delta_1))*((delta_N - delta_Nm1*params.t*GL_d)/(delta_Np1 - delta_N*params.t*GL_d))*delta_1;

	dcomp Rsigma_0_u, Rsigma_0_d, Rsigma_PI_u, Rsigma_PI_d;
	dcomp Fsigma;
	Rsigma_0_u = (1.-GR_0_u*params.t*GN_u*params.t);
	Rsigma_0_d = (1.-GR_0_d*params.t*GN_d*params.t);
	Rsigma_PI_u = (1.-GR_PI_u*params.t*GN_u*params.t);
	Rsigma_PI_d = (1.-GR_PI_d*params.t*GN_d*params.t);
	Fsigma = (1./M_PI)*log(Rsigma_0_d*Rsigma_0_u/(Rsigma_PI_u*Rsigma_PI_d));
	
	dcomp i;
	i = -1.;
	i = sqrt(i);
	dcomp factor = exp(-i*2.*M_PI*static_cast<dcomp>(params.s*N)/params.P);
	/* cout<<"factor = "<<factor<<"N, s = "<<N<<" "<<params.s<<endl; */
	Fsigma *= factor;

	return (1/params.P)*Fsigma;
}

dcomp integration(a_struct &params)
{
	double start = 18.;
	double end = start + params.P;
	int n = 150;
	double h = (end - start)/(2.*n);
	dcomp result = 0.;
	double N = start;
	result += Rspace(N, params);
	N = end;
	result += Rspace(N, params);
	#pragma omp parallel
	{
	#pragma omp for nowait reduction(+:result)
	for (int j = 1; j <= n; j++){
		N = start + (2*j - 1.)*h;	
		result += 4.*Rspace(N, params);
		if (j != n){
			N = start + 2*j*h;
			result += 2.*Rspace(N, params);
		}
	}}
	result *= h/3.;
	/* cout<<"s = "<<params.s<<", c_s = "<<result<<endl; */
	return result;
}

double dsi_dE(double psi, a_struct &params)
{
	
	const double delta = 5e-3;
	params.E -= delta;
	double F = 0.;//cos(params.x*params.a)+cos(params.z*params.a);
	double y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	params.P = M_PI/y;
	dcomp delta_E_n = integration(params);
	params.E += delta;

	params.E += delta;
	y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	params.P = M_PI/y;
	dcomp delta_E_p = integration(params);
	params.E -= delta;	

	params.E -= 1.2*delta;
	y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	params.P = M_PI/y;
	dcomp delta_E_2n = integration(params);
	params.E += 1.2*delta;

	params.E += 1.2*delta;
	y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	params.P = M_PI/y;
	dcomp delta_E_2p = integration(params);
	params.E -= 1.2*delta;

	double psi_n = atan2(imag(delta_E_n),real(delta_E_n));

	double psi_p = atan2(imag(delta_E_p),real(delta_E_p));

	double psi_2n = atan2(imag(delta_E_2n),real(delta_E_2n));

	double psi_2p = atan2(imag(delta_E_2p),real(delta_E_2p));
/* cout<<psi<<'\t'<<psi_n<<'\t'<<psi_p<<endl; */
	double dsi_dE = (psi_p - psi_n)/(2*delta);
	double test = (psi_2p - psi_2n)/(2.4*delta);
	double abs_err = abs(abs(dsi_dE/test)-1);
	/* cout<<1/sinh(dsi_dE)<<endl; */

	if ( abs_err > 0.05)
		cout<<"caution: for s = "<<params.s<<", dpsi/dE has absolute error of "<<abs_err<<endl;
	if ( abs(abs(2.*psi/(psi_p + psi_n))-1) > 0.05)
		cout<<"caution: for s = "<<params.s<<", psi is not on the gradient line"<<endl;

	/* cout<<params.s<<" "<<dsi_dE<<endl; */
	return dsi_dE;
}

/* Matrix2d hessian(a_struct &params) */
/* { */
/* 	double Cz = cos(params.a*params.z); */
/* 	double Sz = sin(params.a*params.z); */
/* 	double Cx = cos(params.a*params.x); */
/* 	double Sx = sin(params.a*params.x); */
/* 	double disp = (params.v2-real(params.E))/(2.*params.t) + Cx + Cz; */
/* 	double el_11 = -params.a*Cx/sqrt(-disp*disp +1.) + params.a*Sx*Sx*disp/pow((-disp*disp + 1.),1.5); */
/* 	double el_22 = -params.a*Cz/sqrt(-disp*disp +1.) + params.a*Sz*Sz*disp/pow((-disp*disp + 1.),1.5); */
/* 	double el_12 = params.a*Sx*Sz*disp/pow((-disp*disp +1.),1.5); */
/* 	Matrix2d Hes; */
/* 	Hes << el_11,	el_12, */
/* 	    	el_12,	el_22; */
/* 	/1* cout<<el_11<<'\t'<<el_12<<'\n'<<el_12<<'\t'<<el_22<<'\n'<<endl; *1/ */
/* 	return Hes; */
/* } */
		
int main() 
{
	//FOR C_S SIGN OF IMAGINARY COMPONENT IS DEPENDENT ON SIGN OF S. THIS SHOULD NOT BE THE CASE, NOT THE CASE FOR REAL COMPONENT OR ANDREY'S REAL OR IMAGINARY.
//why doesn't dsi_dE contribute to spincurrent?
	a_struct params;
	double N = 50;

	dcomp i;
	i = -1.;
	i = sqrt(i);
	params.E = 0.0 + i*1e-10;

	//dy/dE = -1/sqrt(-(E+v2-F)^2+1
	//or more obviously.. since E = u + 2t(cos(x)+cos(y)+cos(z), dE/dy = -2tsin(y)..
	double F = 0.;//cos(params.x*params.a)+cos(params.z*params.a);
	const double y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	double dy_dE = -1./(2.*params.t*sin(y));

	// see SPA check for J(N).mw... d derived from this worksheet by solving determinant 
	// of Hessian matrix, only valid for Ef and v2 as given in worksheet
	
	/* Matrix2d Hes = hessian(params); */
	/* double det_hes = sqrt(abs(Hes(0,0)*Hes(1,1)-Hes(0,1)*Hes(1,0))); */
	/* /1* cout<<det_hes<<endl; *1/ */
	/* dcomp tau; */
	/* /1* cout<<Hes.trace()<<endl; *1/ */
	/* SelfAdjointEigenSolver<Matrix2d> es; */
	/* es.compute(Hes); */
	/* VectorXd eigenvals = es.eigenvalues(); */

	/* if ((eigenvals(0) < 0) && (eigenvals(1) < 0)) */
	/* 	tau = -i; */
	/* if ((eigenvals(0) > 0) && (eigenvals(1) > 0)) */
	/* 	tau = i; */
	/* if (((eigenvals(0) > 0) && (eigenvals(1) < 0)) || ((eigenvals(0) < 0) && (eigenvals(1) > 0))) */
	/* 	tau = 1.; */

	int fourier_start = -5;
	int fourier_end = 0;
	int size = (fourier_end-fourier_start);
	VectorXcd C_s;
	VectorXd dPsi_dE;
	C_s.resize(size);
	dPsi_dE.resize(size);
	
	dcomp c_s;

	int index;
	for (params.s = fourier_start; params.s != fourier_end; params.s++){
		index = params.s - fourier_start;
		if (params.s == 0){
			C_s[index] = 0;
			dPsi_dE[index] = 0;
			continue;}
		params.P = M_PI/y;
		c_s = integration(params);
		double psi = atan2(imag(c_s),real(c_s));
		double dpsi_dE = dsi_dE(psi, params);
		C_s[index] = c_s;
		dPsi_dE[index] = dpsi_dE;
	}

	double spincurrent, result;

	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );
	Myfile<<"N , J(N)"<<endl;

	/* double j = 20.; */

	for (double j = 1; j<=N; j = j + 0.1){
		spincurrent = 0.;
		// s here is the summation of Fourier terms
		
		for (params.s = fourier_start; params.s != fourier_end; params.s++){
			index = params.s - fourier_start;
			result = real(C_s[index]*exp(i*static_cast<double>(j*params.s)*2.*y)
				/(sinh(M_PI*params.kT*(2.*j*params.s*dy_dE + dPsi_dE[index]))));
			spincurrent += result;
			/* cout<<"s: "<<params.s<<"\tc_s: "<<C_s[index]<<endl; */
		}
		Myfile<< j-1 <<" , "<< -(params.kT*M_PI)*spincurrent << endl;

	}

	cout<<"finished!"<<endl;


return 0;
}
