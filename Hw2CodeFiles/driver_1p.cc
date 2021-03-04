#include "sol_rk4d.hh"
#include "sol_rk4classic.hh"
#include "brusselator.hh"
#include "TwoComponentSystem.hh"
#include <stdio.h>
#include <cstdio>


void logspace(double* array, int numberPoints, double minp, double maxp){
	double step = (abs(maxp)-abs(minp))/(double)numberPoints;
	for(int i=0; i<numberPoints; i++){
		array[i] = pow(10.,-1.*(minp + i*step));
	}
}


int main(){

// Demo the FSAL RK4(3) code on brusselators problem
	brus_rk4* demoRK4D = new brus_rk4();
	// print step output only with no dense
	demoRK4D->solve_fixed(20.0, 0.01, true, false, 0.0, 1e-5, 1e-5);
	puts("\n"); 
	// print dense output only
	demoRK4D = new brus_rk4();
	demoRK4D->solve_fixed(20.0, 0.01, false, true, 1.0, 1e-5, 1e-5);
	delete demoRK4D;
	puts("\n");

// Generate Precision Work Plot
	double tmax = 20.0; 
	brus_rk4_reference br;
	br.solve_fixed(tmax, 200000);
 	double ref0=br.q[0];
	double ref1=br.q[1];
	
 	double dt0 = 0.01;
	double dt_dens = 0.0; // value of 0 turns off dense interpolation
	int numberPoints = 100;
	double* lambda = new double[numberPoints];
	logspace(lambda, numberPoints, 3, 13);
	for(int i=0; i<numberPoints; i++){
		double thisLambda = lambda[i];
		brus_rk4* brk4 = new brus_rk4();
		brk4->solve_fixed(tmax, dt0, false, false, dt_dens, thisLambda, thisLambda);
		double dy0 = ref0 - brk4->q[0];
		double dy1 = ref1 - brk4->q[1];
		double err = sqrt(dy0*dy0 + dy1*dy1);
		printf("%g %d %g", thisLambda, brk4->fcount, err);
		puts("");
		delete brk4;	
	}
	puts("\n");

	
// Run the rk4d dense solver on the two component system
	tmax = 8.0;
	double Atol_ = 3e-3; 
	double Rtol_ = 3e-3;
	dt_dens = 8./1200.;
	dt0 = 0.01;
	// first print just integration steps
	twoComponent_rk4* rk4sol_steps = new twoComponent_rk4();
	rk4sol_steps->solve_fixed(tmax, dt0, true, false, 0., Atol_, Rtol_);
	delete rk4sol_steps;
	puts("\n");
	// print dense output only
	twoComponent_rk4* rk4sol_dens = new twoComponent_rk4();
	rk4sol_dens->solve_fixed(tmax, dt0, false, true, dt_dens, Atol_, Rtol_);
	delete rk4sol_dens;			
}
















