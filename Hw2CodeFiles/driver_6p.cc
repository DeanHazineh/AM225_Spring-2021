
#include "TwoComponentSystem.hh"
#include <stdio.h>
#include <cmath>
#include "starOrbit.hh"

int main(){

// Run the Geng implicit RK 5 Method on the two-component problem with the usual test conditions
	twoComponent_irk5* solve_irk5 = new twoComponent_irk5();
	solve_irk5->solve_fixed(10.0, 10000, true);
	delete solve_irk5;
	puts("\n");
	
// Run Convergence Analysis
	for(int i=0; i<200; i++){
		int steps=int(100.*pow(100.,i*0.005)); // line mimiced via Rycroft extrap_conv.cc code
		twoComponent_irk5* wpSolve = new twoComponent_irk5();
		wpSolve->solve_fixed(10.0, steps, false);
		double t_ = wpSolve->t;
		double err = sqrt( pow(wpSolve->q[0]-cos(t_*t_/2.),2.) + pow(wpSolve->q[1]-sin(t_*t_/2.),2.) );
		printf("%d %g", wpSolve->fcount, err);
		puts("");
		delete wpSolve;
	}
	puts("\n");
	

// Simulate the galaxy ODE System
	double simtime = 2000.;
	double stepsize = 1./20.;
	orbit_irk5* orbitSolver_irk5 = new orbit_irk5();
	orbitSolver_irk5->solve_fixed(simtime, simtime/stepsize, true);
	puts("\n");
	
// Simulate for longer time and generate Poincare Map Data
	simtime = 1e5;
	stepsize = 1./20.;
	orbitSolver_irk5 = new orbit_irk5();
	orbitSolver_irk5->solve_Poincaire(simtime, simtime/stepsize);
	delete orbitSolver_irk5;
	//
}
