/* this program is a test which contains the best optimised integration
 * developed to date (5/10/16) by myself. It utilises double integration
 * over the two dimensional Brillouin zone by the Cunningham Points (CP)
 * method over the irreducible segment making it only useful for the
 * simple cubic case. The integral is adaptive with obvious lower bound, 
 * upper bound of sampling points, and absolute error observable from 
 * comments below */
//addendum 30-10-16 integration algorithm is now moved to a header file
//which optionally takes error as an input
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
/* #include "cunningham_points_integration.h" */
/* #include "cunningham_points_adaptive.h" */
/* #include "cunningham.h" */
#include "cunningham_spawn.h"

using namespace std;
typedef complex<double> dcomp;
int power;

/* struct a_struct{ */
/* 	const double a = 1.; */
/* 	double j; */
/* }; */

double f(double x, double y, const double a, const double j) {
	/* double E = params.j; */
	double Pi(2.*asin(1));
	double u,W,U0,T,om;
	dcomp mu,G,i, nu;
	u=0.; U0=0., T=.5;
	W=u+2.*T*(cos(x*a)+cos(y*a));
	i=-1;
	i=sqrt(i);
	om=j-W;
	mu=sqrt(4.*T*T-om*om);
	double test=4.*T*T-om*om;
	if (test >= 0){
		G=(i/mu)*(1.+pow((om+i*mu)/(2.*T),power*2)*((mu-i*(om-2.*U0))/(mu+i*(om-2.*U0))));
	}
	if (test < 0){
		if (om >= 0){
			nu=sqrt(om*om-4.*T*T);}
		if (om < 0){
			nu=-sqrt(om*om-4.*T*T);}
		G=(1./nu)*(1.+pow((om-nu)/(2.*T),power*2)*((nu-(om-2.*U0))/(nu+(om-2.*U0))));
	}
	return imag(G);
}

int main() {

	double E_start,E_end;
	/* a_struct params; */
	E_start = -2.99;
	E_end = 2.99;
	const double a = 1.;

	/* cout<<"Enter the atomic plane n, for the DOS you wish to calculate"<<endl; */
	/* cin>>power; */
	/* cout<<"\nName the data file\n"; */

	power = 8;

	string Mydata = "test";

	/* cin.ignore(); getline(cin, Mydata); */

	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );

	double result = 0.;
	for (double j=E_start;j<=E_end;j=j+0.01) {
		/* result = kspace(params,&f,0,0); */	
		/* result = kspace(&f, 3, 0, 48, a, j); */
		result = kspace(&f, 3, 0, a, j);
		result =a*a*result/(4.*M_PI*M_PI*M_PI);
		Myfile<<j<<" , "<<result<<endl;
	}

	cout<<"finished!"<<endl;

return 0;
}

