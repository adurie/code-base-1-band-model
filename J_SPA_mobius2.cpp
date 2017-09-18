#include <iostream>
#include <utility>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
// This is certainly not the most efficient form this could be in... see eqn 2 of Andreys phase paper
using namespace Eigen;
using namespace std;

typedef complex<double> dcomp;

struct a_struct{
	const double a = 1.;
	const double v2 = -2.8;
	const double kT = 8.617342857e-5*300/13.6058;
	dcomp E;
	const double x = 0.;
	const double z = 0.;
	int s;
	double P;	//period of J against spacer thickness
	const double t = 0.5;
};

dcomp greens(dcomp OM, a_struct &params)
{
	Matrix2cd X,O;
	X << 	0,	1/params.t,
		-params.t,	OM/params.t;
	ComplexEigenSolver<Matrix2cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	dcomp b = O(0,1);
	dcomp d = O(1,1);
	dcomp GR = b/d;
	/* cout<<GR(0,0)<<GR(0,1)<<GR(1,0)<<GR(1,1)<<endl; */

	return GR;
}

dcomp Fgreens(dcomp OM, a_struct &params, dcomp G)
{
	Matrix2cd O, fO;

	Matrix2cd X;
	X << 	0,	1/params.t,
		-params.t,	OM/params.t;
	ComplexEigenSolver<Matrix2cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();

	fO = O.inverse();
	dcomp f = (fO(0,0)*G + fO(0,1))/(fO(1,0)*G + fO(1,1));
	/* cout<<f(0,0)<<f(0,1)<<f(1,0)<<f(1,1)<<endl; */
	return f;
}

dcomp Rspace(a_struct &params) {
//...F(0)|NM(n)|F(theta)...
	const double v1 = -1.6295;
	/* const double v1 = -2.8; */
	const double v3 = -1.6295;
	/* const double v3 = -2.8; */
	const double nab = -0.6795;
	double F = cos(params.x*params.a)+cos(params.z*params.a);

	double T = params.t;
	double V2 = params.v2;

	const double v3_u = v3+nab;
	const double v3_d = v3-nab;
	const double v1_u = v1+nab;
	const double v1_d = v1-nab;

	dcomp OM = params.E - V2 - 2.*T*F;
//initialise surface Green's function using mobius transformation for rotated layer at theta = 0
	dcomp OMV3=params.E-v3_u-2.*T*F;
	dcomp GR_0;
	GR_0 = greens(OMV3, params);
	dcomp fR_u;
	fR_u = Fgreens(OM, params, GR_0);

//initialise surface Green's function using mobius transformation for rotated layer at theta = PI 
	dcomp OMV2=params.E-v3_d-2.*T*F;
	dcomp GR_PI;
	GR_PI = greens(OMV2, params);
	dcomp fR_d;
	fR_d = Fgreens(OM, params, GR_PI);

//initialise surface Green's function using mobius transformation for fixed layer 
	dcomp OMV1=params.E-v1_u-2.*T*F;
	dcomp GL_u;
	GL_u = greens(OMV1, params);
	dcomp fL_u;
	fL_u = Fgreens(OM, params, GL_u);

//initialise surface Green's function using mobius transformation for fixed layer 
	dcomp OMV0=params.E-v1_d-2.*T*F;
	dcomp GL_d;
	GL_d = greens(OMV0, params);
	dcomp fL_d;
	fL_d = Fgreens(OM, params, GL_d);

	dcomp Fret = (1./M_PI)*(-fR_u*fL_u - fR_d*fL_d + fR_u*fL_d + fR_d*fL_u);

	/* cout<<Fret<<endl; */
	return Fret;
	/* return -0.04068619021-0.08651454501*i;// Maple v = -1.6295 out */
	/* return 0.1025450815+0.1034204094*i;//Maple v = -2.8 out */
}

dcomp integration(a_struct &params)
{
	return Rspace(params);
}

double dsi_dE(double psi, a_struct &params, dcomp delta_E_p, dcomp delta_E_2p, dcomp delta_E_n, dcomp delta_E_2n, double delta)
{
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
	return Hes;
}
		
int main() 
{
	a_struct params;
	double N = 50;

	dcomp i;
	i = -1.;
	i = sqrt(i);
	params.E = 0.0 + i*1e-4;

	//dy/dE = -1/sqrt(-(E+v2-F)^2+1
	//or more obviously.. since E = u + 2t(cos(x)+cos(y)+cos(z), dE/dy = -2tsin(y)..
	double F = cos(params.x*params.a)+cos(params.z*params.a);
	double y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	double dy_dE = -1./(2.*params.t*sin(y));

	// see SPA check for J(N).mw... d derived from this worksheet by solving determinant 
	// of Hessian matrix, only valid for Ef and v2 as given in worksheet
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
	int fourier_start = 1;
	int fourier_end = 6;
	int size = (fourier_end-fourier_start);
	VectorXcd C_s;
	VectorXd dPsi_dE;
	C_s.resize(size);
	dPsi_dE.resize(size);
	
	dcomp c_s = integration(params);
	dcomp c_ss;

	const double delta = 1e-2;
	params.E -= delta;
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
	y = (1./params.a)*acos((real(params.E)-params.v2)/(2.*params.t)-F);
	params.P = M_PI/y;

	dcomp dp, d2p, dn, d2n;

	int index;
	for (params.s = fourier_start; params.s != fourier_end; params.s++){
		index = params.s - fourier_start;
		if (params.s == 0){
			C_s[index] = 0;
			dPsi_dE[index] = 0;
			continue;}
		c_ss = (1./params.s)*pow(c_s,params.s);
		dp = (1./params.s)*pow(delta_E_p,params.s);
		d2p = (1./params.s)*pow(delta_E_2p,params.s);
		dn = (1./params.s)*pow(delta_E_n,params.s);
		d2n = (1./params.s)*pow(delta_E_2n,params.s);

		double psi = atan2(imag(c_ss),real(c_ss));
		double dpsi_dE = dsi_dE(psi, params, dp, d2p, dn, d2n, delta);
		C_s[index] = c_ss;
		dPsi_dE[index] = dpsi_dE;
		cout<<c_ss<<endl;
	}

	double spincurrent, result;

	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );
	Myfile<<"N , J(N)"<<endl;

	for (double j = 1; j<=N; j = j + 0.1){
		spincurrent = 0.;
		// s here is the summation of Fourier terms
		
		for (params.s = fourier_start; params.s != fourier_end; params.s++){
			if (params.s == 0)
				continue;
			index = params.s - fourier_start;
			result = real(tau*C_s[index]*exp(i*static_cast<double>(j*params.s)*2.*y)//*exp(-i*2.*static_cast<double>(params.s)*y)
				/(2*det_hes*sinh(M_PI*params.kT*(2.*j*params.s*dy_dE//removed s in denominator 1/8/17 as duplicated in loop above
							/* )))); */
						       	+ dPsi_dE[index]))));
			spincurrent += result;
		}
		Myfile<< j <<" , "<< -(params.kT/(2.*j))*spincurrent << endl;

	}

	cout<<"finished!"<<endl;


return 0;
}
