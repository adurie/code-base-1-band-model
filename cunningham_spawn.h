#ifndef CUNNINGHAM_H
#define CUNNINGHAM_H

#include <cmath>
#include <eigen3/Eigen/Dense>
//this is the first TRUE cunningham points algorithm I have made
//with automatic convergence checks. Parameters include the max
//number of points in one meridian (iterations), the relative 
//error value, and the number of points in one meridian to start
//with. 27-06-17

using namespace std;

double Condition(std::complex<double> result, std::complex<double> part_result){
	double condition;
	if ((std::abs(part_result) == 0) && (std::abs(result) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(part_result)/std::abs(result)-1.);
	return condition;
}
double Condition(double result, double part_result){
	double condition;
	if ((std::abs(part_result) == 0) && (std::abs(result) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(std::abs(part_result)/(result))-1.);
	return condition;
}
double Condition(const Eigen::VectorXd &result, const Eigen::VectorXd &part_result){
	double condition;
	if ((std::abs(part_result.sum()) == 0) && (std::abs(result.sum()) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(part_result.sum())/std::abs(result.sum())-1.);
	return condition;
}
double Condition(const Eigen::VectorXcd &result, const Eigen::VectorXcd &part_result){
	double condition;
	if ((std::abs(part_result.sum()) == 0) && (std::abs(result.sum()) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(part_result.sum())/std::abs(result.sum())-1.);
	return condition;
}

template <typename func, typename... Args, typename ret>// a square can only spawn a square
ret aux_square(func&& predicate, int depth, double error, int n, double v, double w, const ret& f, double A, const double a, Args&&... params)
{
	ret f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9;
	double x, z;
	x = v;		z = w;	
	f_1 = (1/9.)*f;

	x = v;	z = w + 2*A;	
	/* cout<<x<<", "<<z<<endl; */
	f_2 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v;	z = w - 2*A;	
	/* cout<<x<<", "<<z<<endl; */
	f_3 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v + 2*A;	z = w;	
	/* cout<<x<<", "<<z<<endl; */
	f_4 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v + 2*A;	z = w + 2*A;	
	/* cout<<x<<", "<<z<<endl; */
	f_5 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v + 2*A;	z = w - 2*A;	
	/* cout<<x<<", "<<z<<endl; */
	f_6 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v - 2*A;	z = w;	
	/* cout<<x<<", "<<z<<endl; */
	f_7 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v - 2*A;	z = w + 2*A;	
	/* cout<<x<<", "<<z<<endl; */
	f_8 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v - 2*A;	z = w - 2*A;	
	/* cout<<x<<", "<<z<<endl; */
	f_9 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	ret f_total;
	f_total = f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8 + f_9;
	/* cout<<f_total<<" "<<f<<" square"<<endl; */

	double rel_error;
	rel_error = Condition(f_total, f);//condition overloaded for all data types this header is designed for
		/* cout<<rel_error<<" tri"<<endl; */

	if (n>90){//avoids early convergence
		if ((rel_error <= error) || (depth <= 0))
			return f_total;
	}

	if ((rel_error > 0.7) && (f_total > 1e-4))//if the result is too erroneous, the function keeps iterating until convergence is achieved.
		depth++;// This is because if neighbouring points are zero, the result will appear to 
			//halve each time as each recursion is scaled, until convergence.
			
	//begin spawning points around points already integrated over
	ret g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9;

	g_1 = aux_square(predicate, depth-1, error, n*9, v, w, f_1, A/3., a, params...);
	g_2 = aux_square(predicate, depth-1, error, n*9, v, w + 2*A, f_2, A/3., a, params...);
	g_3 = aux_square(predicate, depth-1, error, n*9, v, w - 2*A, f_3, A/3., a, params...);
	g_4 = aux_square(predicate, depth-1, error, n*9, v + 2*A, w, f_4, A/3., a, params...);
	g_5 = aux_square(predicate, depth-1, error, n*9, v + 2*A, w + 2*A, f_5, A/3., a, params...);
	g_6 = aux_square(predicate, depth-1, error, n*9, v + 2*A, w - 2*A, f_6, A/3., a, params...);
	g_7 = aux_square(predicate, depth-1, error, n*9, v - 2*A, w, f_7, A/3., a, params...);
	g_8 = aux_square(predicate, depth-1, error, n*9, v - 2*A, w + 2*A, f_8, A/3., a, params...);
	g_9 = aux_square(predicate, depth-1, error, n*9, v - 2*A, w - 2*A, f_9, A/3., a, params...);

	ret g_total;
	g_total = g_1 + g_2 + g_3 + g_4 + g_5 + g_6 + g_7 + g_8 + g_9;
	return g_total;
}

template <typename func, typename... Args, typename ret>// a triangle can spawn square and triangles, but only triangles when x = z
ret aux_tri(func&& predicate, int depth, double error, int n, double v, double w, const ret& f, double A, const double a, Args&&... params)
{
	ret f_1, f_2, f_3, f_4, f_5, f_6;
	double x, z;
	x = v;		z = w;	
	f_1 = (1/9.)*f;

	x = v + 2*A;	z = w;	
	f_2 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v + 2*A;	z = w + 2*A;	
	f_3 = (0.5/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v;	z = w - 2*A;	
	f_4 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v + 2*A;	z = w - 2*A;	
	f_5 = (1/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = v - 2*A;	z = w - 2*A;	
	f_6 = (0.5/(n*9.))*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	ret f_total;
	f_total = f_1 + f_2 + f_3 + f_4 + f_5 + f_6;
	/* cout<<f_total<<" "<<f<<" tri"<<endl; */
	/* cout<<f_1<<" "<<f_2<<" "<<f_3<<" "<<f_4<<" "<<f_5<<" "<<f_6<<endl; */

	double rel_error;
	rel_error = Condition(f_total, f);//condition overloaded for all data types this header is designed for
		/* cout<<rel_error<<" tri"<<endl; */

	if (n>90){//avoids early convergence
		if ((rel_error <= error) || (depth <= 0))
			return f_total;
	}

	if ((rel_error > 0.7) && (f_total > 1e-4))//if the result is too erroneous, the function keeps iterating until convergence is achieved.
		depth++;// This is because if neighbouring points are zero, the result will appear to 
			//halve each time as each recursion is scaled, until convergence.

	//begin spawning points around points already integrated over
	ret g_1, g_2, g_3, g_4, g_5, g_6;

	g_1 = aux_tri(predicate, depth-1, error, n*9, v, w, f_1, A/3., a, params...);
	g_2 = aux_square(predicate, depth-1, error, n*9, v + 2*A, w, f_2, A/3., a, params...);
	g_3 = aux_tri(predicate, depth-1, error, n*9, v + 2*A, w + 2*A, f_3, A/3., a, params...);
	g_4 = aux_square(predicate, depth-1, error, n*9, v, w - 2*A, f_4, A/3., a, params...);
	g_5 = aux_square(predicate, depth-1, error, n*9, v + 2*A, w - 2*A, f_5, A/3., a, params...);
	g_6 = aux_tri(predicate, depth-1, error, n*9, v - 2*A, w - 2*A, f_6, A/3., a, params...);

	ret g_total;
	g_total = g_1 + g_2 + g_3 + g_4 + g_5 + g_6;
	return g_total;
}

template <typename func, typename... Args>
auto kspace(func&& predicate, int depth, double rel, const double a, Args&&... params)
    -> decltype(std::forward<func>(predicate)(std::declval<double>(), std::declval<double>(), std::declval<double>(),std::forward<Args>(params)...)) {
	    typedef decltype(std::forward<func>(predicate)(std::declval<double>(), std::declval<double>(), std::declval<double>(),std::forward<Args>(params)...)) ret;

	double error;
	if (rel == 0)
		error = 0.02;
	else
		error = rel;
	int max_width;

	/* string Mydata = "test.txt"; */
	/* ofstream Myfile; */	
	/* Myfile.open( Mydata.c_str(),ios::trunc ); */

	double A = M_PI/a;
	double x,z;
	ret f;
	double rel_error;
	
	x = A/2.;
	z = A/2.;
	//calculate the mean value point only	
	f = 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	//calculate the 5 remaining points in the 6 point grid
	ret f_1, f_2, f_3, f_4, f_5, f_6;
	A /= 6.;

	x = A;		z = A;	
	f_1 = (0.5/9.)*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = 3*A;	z = A;	
	f_2 = (1/9.)*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = 5*A;	z = A;	
	f_3 = (1/9.)*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = 5*A;	z = 3*A;	
	f_4 = (1/9.)*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	x = 3*A;	z = 3*A;	
	f_5 = (1/9.)*f;

	x = 5*A;	z = 5*A;	
	f_6 = (0.5/9.)*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);

	ret f_total;
	f_total = f_1 + f_2 + f_3 + f_4 + f_5 + f_6;

	/* cout<<f_total<<" "<<f<<endl; */
	/* cout<<f_1<<" "<<f_2<<" "<<f_3<<" "<<f_4<<" "<<f_5<<" "<<f_6<<endl; */
	/* rel_error = Condition(f_total, f);//condition overloaded for all data types this header is designed for */

	/* if (rel_error <= error){ */
	/* 	A = M_PI/a; */
	/* 	return f_total*8*A*A; */
	/* } */

	//begin spawning points around points already integrated over
	ret g_1, g_2, g_3, g_4, g_5, g_6;

	g_1 = aux_tri(predicate, depth, error, 9, A, A, f_1, A/3., a, params...);
	g_2 = aux_square(predicate, depth, error, 9, 3*A, A, f_2, A/3., a, params...);
	g_3 = aux_square(predicate, depth, error, 9, 5*A, A, f_3, A/3., a, params...);
	g_4 = aux_square(predicate, depth, error, 9, 5*A, 3*A, f_4, A/3., a, params...);
	g_5 = aux_tri(predicate, depth, error, 9, 3*A, 3*A, f_5, A/3., a, params...);
	g_6 = aux_tri(predicate, depth, error, 9, 5*A, 5*A, f_6, A/3., a, params...);

	ret g_total;
	g_total = g_1 + g_2 + g_3 + g_4 + g_5 + g_6;
	A = M_PI/a;
	return g_total*8*A*A;
}
#endif
