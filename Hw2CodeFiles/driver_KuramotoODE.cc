#include "KuramotoModel.hh"
#include <cstdio>
#include <random> 


int main(){	
	//testInitializeFunction();		
	
	int N = 1250;
	double t_max = 200.;
	double dt_dense = 1.;
	double timeShow[5] = {0.,10., 20., 50., 200.};
	int len_timeShow = 5;
	double J = 1.0;
	double K = 0.0;
	
	KuramotoModel_Solver* test = new KuramotoModel_Solver(N, timeShow, len_timeShow, true, J, K);
	test->solve_fixed(t_max, .01, false, true, dt_dense, 1e-6, 1e-6);
	printf("Total Steps: %d , Total Fevals: %d", test->stepCounter, test->fcount);
	delete test;

}
