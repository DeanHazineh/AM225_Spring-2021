#ifndef SOL_RK4D_HH
#define SOL_RK4D_HH

/** Class for solving an ODE IVP using the fourth-order Runge-Kutta method. */
class rk4d {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        int stepCounter; 
        
        rk4d(int dof_);
        virtual ~rk4d();
        
        void dense_output(double theta,double dt);
		  void solve_fixed(double tmax, double dt, bool printSteps, bool printDense, double dt_den,
				 double Atol_, double Rtol_);
        void step(double dt);
        
        virtual void print(double t_,double *in, double error)=0;     
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    
		  double error;
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
        double *k5;
        double Atol;
        double Rtol;
        double facmax;
        double facmin;
        double fac;
        
        double updateStep(double dt);

};

#endif
