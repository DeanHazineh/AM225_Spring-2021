# include "discontinuousTwoComponent.hh"
#include "sol_rk4d.hh"
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <cmath>

const double EulerConstant = std::exp(1.0);
double tmax = 48. + 1./EulerConstant;

void logspace(double* array, int numberPoints, double minp, double maxp, double sign){
	double step = (abs(maxp)-abs(minp))/(double)numberPoints;
	for(int i=0; i<numberPoints; i++){
		array[i] = pow(10.,sign*(minp + i*step));
	}
}


// driver code for problem 7 
int main(){
	
	//Demo the classic rk4 code in solving the ODE
	discontinuous_twocomponent_rk4classic* test = new discontinuous_twocomponent_rk4classic();
	test->solve_fixed(tmax, 1000, true);
	delete test;
	puts("\n");
	
	// Demo the FSAL RK4(3) code in solving the ODE
	discontinuous_twocomponent_rk4d* rk4d = new discontinuous_twocomponent_rk4d();
	rk4d->solve_fixed(tmax, 0.01, true, false, 0., 1e-6, 0.);
	delete rk4d;
	puts("\n");	
	
	// Compute the work precision plot for classic rk4 method
	int numberPoints = 100;
	double* numsteps = new double[numberPoints];
	logspace(numsteps, numberPoints, 3, 7, 1.);
	for(int i=0; i<numberPoints; i++){
		//printf("%g \n", numsteps[i]);
		discontinuous_twocomponent_rk4classic* test = new discontinuous_twocomponent_rk4classic();
		test->solve_fixed(tmax, numsteps[i], false);
		int fevals = test->fcount;
		double errorx = test->q[0] - test->xanalytic(test->t);
		double errory = test->q[1] - test->yanalytic(test->t);
		printf("%d %g ", fevals, sqrt(errorx*errorx + errory*errory));
		puts("");
		delete test;	
	}
	puts("\n"); // split the line up
	
	// Use our FSAL Dense RK4(3) method to get work-precision curve
	int numberPoints2 = 100;
	double* Atol = new double[numberPoints2];
	logspace(Atol, numberPoints2, 2, 12, -1);
	for(int i=0; i<numberPoints2; i++){
		//printf("%g \n", Atol[i]);
		discontinuous_twocomponent_rk4d* rk4d = new discontinuous_twocomponent_rk4d();
		rk4d->solve_fixed(tmax, 0.01, false, false, 0., Atol[i], 0.);
 		double fevals = rk4d->fcount;
		double errorx = rk4d->q[0] - rk4d->xanalytic(rk4d->t);
		double errory = rk4d->q[1] - rk4d->yanalytic(rk4d->t);
		printf("%g %g %g ", Atol[i], fevals, sqrt(errorx*errorx + errory*errory));
		puts("");
		delete rk4d;	
	}
	
}
























