#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <blaze/Math.h>

using namespace blaze;
using namespace std;
typedef complex<double> dcomp;

//this is finished and working. Utilises Cunningham points integration, KS method with array return.
DynamicVector<double> f(double x, double y, int N, double a) {
	double Pi(2.*asin(1));
	double F,T,V1,e,h,V2,V3;
	dcomp mu,GL,i, nu,GN,om,om3,mu3,nu3,GR;
	e = 1.6e-19;
	h = 6.626e-34;

	T=.5;  V1=-1.9; V2=-2.9; V3 = -1.1;
	double t = T;
	F = cos(x*a)+cos(y*a);
	i=-1;
	i=sqrt(i);
	dcomp E=(1e-8)*i;
	om=E-V1-2.*T*F;

	double test=4.*t*t-real(om)*real(om);
	if (test >= 0){
		nu=sqrt(om*om-4.*t*t);
	}
	if (test < 0){
		if (real(om) >= 0){

			nu=sqrt(om*om-4.*t*t);
		}
		if (real(om) < 0){

			nu=-sqrt(om*om-4.*t*t);
		}
	}

	GL=2./(om + nu);

	om3=E-V3-2.*T*F;

	double test3=4.*t*t-real(om3)*real(om3);
	if (test3 >= 0){
		nu3=sqrt(om3*om3-4.*t*t);
	}
	if (test3 < 0){
		if (real(om3) >= 0){

			nu3=sqrt(om3*om3-4.*t*t);
			}
		if (real(om3) < 0){

			nu3=-sqrt(om3*om3-4.*t*t);
		}
	}

	GR=2./(om3 + nu3);

	dcomp om2 = E - V2 - 2.*T*F;
	DynamicVector<dcomp> G_Vector(N);
	
	GN = 1./(om2-T*T*GL);
	G_Vector[0] = GN;
	for (int it=1; it < N; ++it){

		GN=1./(om2-T*T*GN);
		G_Vector[it] = GN;
	}

	DynamicVector<double> conductance(N);
	for (int J=0; J<N;++J){
	conductance[J] = e*e*(1./h)*real((1./(Pi*Pi))*T*T*imag(G_Vector[J])*imag(GR)/((1.-T*T*G_Vector[J]*GR)*(1.-T*T*conj(G_Vector[J])*conj(GR))));
	}

	return conductance;
	
}

int main() {

	double x,y,a;
	double Pi=2*asin(1);
	//number of atomic planes
	int N=50;
	a = 1.;
	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".csv";
	Myfile.open( Mydata.c_str(),ios::app );

	//begin integration
	double A = Pi/a;

	DynamicVector<double> integral(N);
	integral = 0.;
//integral width
	int n = 1000;

	for (int k=1; k<=n+1; k++) {
		x = A*(k-1)/n;
		for (int l=1; l<=k; l++) {
			y = A*(l-1)/n;
			if (k==l){
				integral = integral + 0.5*f(x,y,N,a);}
			if (k!=l){
				integral = integral + f(x,y,N,a);}
		}
	}
	integral = (8.*Pi*Pi/(n*n))*integral;


	Myfile<<"N , Gamma"<<endl;
	for (int J=1;J<abs(N);++J) {

		Myfile<<J+1<<" , "<<integral[J]<<endl;
	}
	cout<<"finished!"<<endl;

return 0;
}

