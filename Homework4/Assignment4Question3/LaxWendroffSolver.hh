#ifndef LAXWENDROFF_HH
#define LAXWENDROFF_HH
#include <cstdio>
#include <cmath>

class LaxWendroffSol{
	public: 
		const int m; // Define the interval into m domains ci 
		const double dx; // the length of the domains 
		double dt; // timestep
      double *a; // primary solution grid
      double *b; // secondary solution grid
      
		LaxWendroffSol(int m_);
		~LaxWendroffSol();
		void init_step_function();
      void init_ProblemStatement();
		void solve(const char* filename, int snaps, double duration, double safe_fac);
		double* AnalyticPeriodPrediction_init_ProblemStatement(bool print);
		double* AnalyticPeriodicPrediction_init_StepFunction(bool print);
		double computeError(double duration, int initSwitch);
	private:
		double Afunc(int j); // defines the advection velocity spatial formula for problem
		void step(double dt);
		void print_line(FILE *fp,double x,double *zp,int snaps);

};


#endif
