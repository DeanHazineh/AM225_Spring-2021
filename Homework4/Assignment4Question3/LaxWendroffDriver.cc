#include <cstdio>
#include <cstring>
#include <math.h>
#include "LaxWendroffSolver.hh"

void run2DLaxWendroff(const char* filename);
void run1DLinesLaxWendroff(const char* filename);
void runConvergence();


int main(){
   //run2DLaxWendroff("Out/2DSpaceTime");
   //run1DLinesLaxWendroff("Out/4LinesSolution");
   runConvergence();
}

void runConvergence(){
	for(int m=256; m<1248; m=m+50){
		const double T = 3.*M_PI/sqrt(5.);	
		LaxWendroffSol* test = new LaxWendroffSol(m);
		test->init_step_function();
		//test->init_ProblemStatement();
		
		double L2 = test->computeError(T,1);
		printf("%d %g \n", m, L2);
		delete test;
	}
}


void run2DLaxWendroff(const char* filename){
	const int snaps=50;    
	const int m=512;
	const double sf=1.0;
	const double T = 3.*M_PI/sqrt(5.);
	
	LaxWendroffSol* test = new LaxWendroffSol(m);
	test->init_step_function();
	//test->init_ProblemStatement();
	test->solve(filename, snaps, T, sf);
	delete test; 
}

void run1DLinesLaxWendroff(const char* filename){
	const int snaps=4;    
	const int m=512;
	const double sf=1.0;
	const double T = 3.*M_PI/sqrt(5.);
	
	LaxWendroffSol* test = new LaxWendroffSol(m);
	test->init_step_function();
	//test->init_ProblemStatement();
	test->solve(filename, snaps, T, sf);
	//test->AnalyticPeriodPrediction_init_ProblemStatement(true);
	test->AnalyticPeriodicPrediction_init_StepFunction(true);
	delete test;
}


