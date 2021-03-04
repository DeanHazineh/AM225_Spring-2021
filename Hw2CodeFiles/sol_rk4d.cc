// Dean Hazineh
// Defines the RK4 Solver Class

#include "sol_rk4d.hh"
#include <cstdio>
#include <cstring>
#include <math.h>  

rk4d::rk4d(int dof_) : dof(dof_), fcount(0), t(0.),stepCounter(0),
	 q(new double[dof]),	dq(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k4(new double[dof]), k5(new double[dof]),
    facmax(3.0), facmin(1./3.), fac(0.9){
		// Initializes the fourth-order Runge-Kutta solver
		// with fixed parameterization for facx
    }

rk4d::~rk4d() {
    delete [] k5;
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
    delete [] q;
}

// This function solves the ODE
// dt_den = 0 turns off dense interpolation
// input is the maximum time to step through
// dt is the initial step size
// output boolean true prints integration steps and output boolean false prints dense if dt_den !=0 
// and the Atolerance and Rtolerance for errors
void rk4d::solve_fixed(double tmax, double dt, bool printSteps, bool printDense, double dt_den,
	 double Atol_, double Rtol_){
	 Atol = Atol_; Rtol = Rtol_;
	
    init();
    double t_den = 0.; // t_den tracks time for dense output state
    if(printSteps||printDense) print(t,q,error);
    
    ff(t,q,k1);
    // step through and solve the ODE
    while(t<tmax){ 
    	  if(t>0.){
    	  	  // update at loop start instead of end to interpolate tmax easier
		     // Move the required data for the next step into position
		     memcpy(q,dq,dof*sizeof(double));
		     memcpy(k1,k2,dof*sizeof(double));
		     if(printSteps) print(t,q, error);
	        dt = updateStep(dt); // defines the dt to use next
	        stepCounter+=1;
    	  }
    	
        step(dt);
        while(error>=1.0){
        		// If error is too high, then discard and redo the step with new stepsize
				//printf("Failed Step, Error: %g \n ", error);
				dt = updateStep(dt);
				step(dt);	
		  }
		  // final step is committed by updating time before dense output 
		  // see rycrofts sol_rk4d example code where time is updated in step originally
	  	  t += dt;
	  	  
        // Do any dense output interpolation
        // Taken from rycrofts example code sol_rk4d
        if(dt_den>0) {
            while(t_den+dt_den<t && t_den+dt_den<tmax){
                t_den+=dt_den;
                dense_output(1.+(t_den-t)/dt,dt);
                if(printDense) print(t_den,k3, 0.0);
            }
        }
          
    }
    // we have solved an accepted (low error) dt which jumps t past tmax
    // explicitly solve the tmax time at the end with interpolation
    // this is printed as the final point
    t_den = tmax;
    dense_output(1.+(t_den-t)/dt,dt);
    t = tmax;
    memcpy(q,k3,dof*sizeof(double));
    if(printSteps||printDense) print(t,q,0.0);
}


void rk4d::dense_output(double theta,double dt) {
	/** Pulled from Rycrofts example code 
	 *	Computes a Hermite interpolation of the solution, for dense output. The
	 * result is stored into the k3 array.
	 * \param[in] theta the fraction of the timestep at which to evaluate the
	 *                  interpolation.
	 * \param[in] dt the length of the current timestep. */
	 // In summary:
    // The function assumes that the current solution is in q, the new solution
    // is in dq, the current derivative is in k1, and the new derivative is in
    // k2
    double mth=1-theta;
    for(int i=0;i<dof;i++)
        k3[i]=mth*q[i]+theta*dq[i]
             -theta*mth*((1-2*theta)*(dq[i]-q[i])+dt*(theta*k2[i]-mth*k1[i]));
}

double rk4d::updateStep(double dt){
	// update the timestep based on current dt and error.
	// The formula is hnew = h min{facmax, max{facmin, fac(q/error)^1/()q+1))}} given on slide 25 of unit 1 lectures- rycroft
	double dtopt = dt*fmin(facmax, fmax(facmin, fac*pow(1./error,1./(3.+1.))) );
	return dtopt;
}

void rk4d::step(double dt) { 
    // Five step FSAL is implemented 
    // The butcher tableu is given on slide 22 of unit 1 lectures - Rycroft
    
    // Second RK Step
    for(int i=0;i<dof;i++) dq[i]=q[i] + (1./3.)*dt*k1[i];
    ff(t+(1./3.)*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i] + dt*(-(1./3.)*k1[i]+1.*k2[i]);
    ff(t+(2./3.)*dt,dq,k3);

    // Fourth RK step
    for(int i=0;i<dof;i++) dq[i]=q[i] + dt*(1.*k1[i]-1.*k2[i]+1.*k3[i]);
    ff(t+1.*dt,dq,k4);
	 
	 // Fifth Rk step
	 for(int i=0;i<dof;i++) dq[i]=q[i] + dt*((1./8.)*k1[i]+(3./8.)*k2[i]+(3./8.)*k3[i]+(1./8.)*k4[i]);
    ff(t+1.*dt,dq,k5); // new y (y1) is exactly the current dq state 
    fcount+=4;
    
    // compute error terms and update global error
    error = 0.0;
    for(int i=0; i<dof; i++){
    	double yhati = q[i] + dt*((1./12.)*k1[i]+(1./2.)*k2[i]+(1./4.)*k3[i]+(1./6.)*k5[i]);
    	double sci = Atol + Rtol*fmax(abs(q[i]),abs(dq[i]));
    	error += pow((dq[i]-yhati)/sci,2.)*(1./dof);
    }
    error = sqrt(error);
    
    // move k5 which is the next derivative into k2 for dense output calc
	 memcpy(k2,k5,dof*sizeof(double));
}













