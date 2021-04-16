#include <cstdlib>
#include <cstring>
#include <math.h> 
#include "LaxWendroffSolver.hh"


LaxWendroffSol::LaxWendroffSol(int m_) : m(m_), dx(2*M_PI/m), a(new double[m]), b(new double[m]){
}
 
LaxWendroffSol::~LaxWendroffSol(){
	delete[] a;
	delete[] b;
}

void LaxWendroffSol::init_step_function(){
	/** Initializes the solution to be a step function. */
	for(int i=0;i<m;i++) {
		double x=dx*(i+0.5);
		a[i]=fmax(M_PI/2.-abs(x-M_PI), 0.);
	}
}

void LaxWendroffSol::init_ProblemStatement(){
	/** Initializes solution to exp(sinx + 1/2 sin(4x)). */
	for(int i=0;i<m;i++) {
		//double x= i==m-1?dx*0.5:dx*(i+0.5);
		double x= dx*(i+0.5);
		a[i]=exp(sin(x) + 0.5*sin(4.*x));
	}
}	

double LaxWendroffSol::Afunc(int j){
	return 2. + 4./3.*sin((j+0.5)*dx);
}

void LaxWendroffSol::step(double dt){
	for(int j=0; j<m; j++){
		// define indices
		int jl= j==0?m-1+j:j-1;
      int jr= j==m-1?1-m+j:j+1;
		// Define Flux Expression explicitly to focus on clarity not speed
		double Fm = 0.5*(Afunc(jl)*a[jl] + Afunc(j)*a[j]) - Afunc(jl+0.5)*dt/2./dx*(Afunc(j)*a[j] - Afunc(jl)*a[jl] );
		double Fp = 0.5*(Afunc(j)*a[j] + Afunc(jr)*a[jr]) - Afunc(jr-0.5)*dt/2./dx*(Afunc(jr)*a[jr] - Afunc(j)*a[j] );
		
		// perform update 
		b[j] = a[j] - dt/dx*(Fp-Fm);
	}
	// After updating the state, swap the pointers
	double *c=a;a=b;b=c;
}


void LaxWendroffSol::print_line(FILE *fp,double x,double *zp,int snaps){
	/** Prints a line of stored snapshots to a file.
	 * \param[in] fp a pointer to the file to write to.
	 * \param[in] x the position in the domain corresponding to this line.
	 * \param[in] zp a pointer to the first snapshot data point to print.
	 * \param[in] snaps the number of snapshots (not including the starting
	 *                  snapshot). */
    fprintf(fp,"%g",x);
    for(int i=0;i<=snaps;i++) fprintf(fp," %g",zp[i*m]);
    fputc('\n',fp);
}


void LaxWendroffSol::solve(const char* filename, int snaps, double duration, double safe_fac){
	/** Solves the transport equation, storing snapshots of the solution to a file.
	 * \param[in] filename the name of the file to write to.
	 * \param[in] snaps the number of snapshots to save (not including the initial
	 *                  snapshot).
	 * \param[in] duration the number of iterations to step the solution forward by
	 *                     between snapshots.
	 * \param[in] safe_fac a safety factor to apply to the CFL timestep restriction.*/
	 
	double c = 2. + 4./3.;
	double interval=duration/snaps;
	double dt=dx/c/3.*safe_fac;
	int iters=static_cast<int>(interval/dt)+1;
	dt=interval/iters;
	double* trackTime = new double[snaps+1];
	trackTime[0] = 0.0;
	double timeElapsed = 0.0;	
	
	// Allocate memory to store solution snapshots. Integrate the system and // store snapshots after fixed time intervals
	double *z=new double[m*(snaps+1)];
	memcpy(z,a,m*sizeof(double));
	for(int i=1;i<=snaps;i++) {
		for(int k=0;k<iters;k++){
			step(dt);
			timeElapsed+=dt;
		}
		// Store the snapshot
		printf("elapsed Time: %g \n", timeElapsed);
		memcpy(z+i*m,a,m*sizeof(double));
		trackTime[i] = timeElapsed;
	}

	// Open the output file to store the snapshots
	FILE *fp=fopen(filename,"w");
	if(fp==NULL) {
	  fputs("Can't open output file\n",stderr);
	  exit(1);
	}
	
	// Save to file but make the first row the time vector
	fprintf(fp,"%g", 999.);
	for(int i=0;i<=snaps;i++) fprintf(fp," %g",trackTime[i]);
	fputc('\n',fp);
	// Now Save the data with the first column the grid x
	print_line(fp,-0.5*dx,z+(m-1),snaps);
	for(int j=0;j<m;j++) print_line(fp,(j+0.5)*dx,z+j,snaps);
	
	// Delete snapshots and close file
	fclose(fp);
	delete [] z;   
	delete [] trackTime;
}


double LaxWendroffSol::computeError(double duration, int initSwitch){
	double c = 2. + 4./3.;
	double interval=duration;
	double dt=dx/c/3.;
	int iters=static_cast<int>(interval/dt)+1;
	dt=interval/iters;
	
	double timeElapsed = 0.0;	
	for(int k=0;k<iters;k++){
		step(dt);
		timeElapsed+=dt;
	}
	
	// Compute the L2 Error and return it
	double* analytic;
	if (initSwitch==0) analytic = AnalyticPeriodPrediction_init_ProblemStatement(false);
	if (initSwitch==1) analytic = AnalyticPeriodicPrediction_init_StepFunction(false);

	double L2=0;
	for(int j=0; j<m; j++) L2 += pow(a[j]-analytic[j],2.);
	delete[] analytic; 
	return L2;
}


double* LaxWendroffSol::AnalyticPeriodPrediction_init_ProblemStatement(bool print){
	// analytic solution
	// Note my analytic derivation is incorrect and I couldn't figure out the way to fix it
	// as a test/to make the most of things, I will use the initial state as the analytic solution
	
	double* analyticSolT = new double[m];
	
	for(int i=0;i<m;i++) {
		double x=dx*(i+0.5);
		//analyticSolT[i] = exp(sin(x) + 0.5*sin(4.*x))*exp(-4./3.*cos(x)*3.*M_PI/sqrt(5.));
		analyticSolT[i] = exp(sin(x) + 0.5*sin(4.*x));
	}
	
	if(print){
		FILE *fp=fopen("Out/AnalyticSolution","w");
		if(fp==NULL) {
		  fputs("Can't open output file\n",stderr);
		  exit(1);
		}

		for(int j=0;j<m;j++){
			// print x in first collumn and value in the second
			fprintf(fp,"%g %g ",dx*(j+0.5), analyticSolT[j]); 
		   fputc('\n',fp);
		}
	}
	
	return analyticSolT;
}


double* LaxWendroffSol::AnalyticPeriodicPrediction_init_StepFunction(bool print){
	// analytic solution
	// Note my analytic derivation is incorrect and I couldn't figure out the way to fix it
	// as a test/to make the most of things, I will use the initial state as the analytic solution
	
	double* analyticSolT = new double[m];
	
	for(int i=0;i<m;i++) {
		double x=dx*(i+0.5);
		//analyticSolT[i] = fmax(M_PI/2.-abs(x-M_PI), 0.)*exp(-4./3.*cos(x)*3.*M_PI/sqrt(5.));
		analyticSolT[i] = fmax(M_PI/2.-abs(x-M_PI), 0.);
	}
	
	if(print){
		FILE *fp=fopen("Out/AnalyticSolution","w");
		if(fp==NULL) {
		  fputs("Can't open output file\n",stderr);
		  exit(1);
		}

		for(int j=0;j<m;j++){
			// print x in first collumn and value in the second
			fprintf(fp,"%g %g ",dx*(j+0.5), analyticSolT[j]); 
		   fputc('\n',fp);
		}
	}
	
	return analyticSolT;
}
















