#ifndef SOL_RK4Classic_HH
#define SOL_RK4Classic_HH

/** Class for solving an ODE IVP using the classic fourth-order Runge-Kutta method. */
class classic_rk4{
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        
        classic_rk4(int dof_);
        virtual ~classic_rk4();
       
        void solve_fixed(double t_end, int steps, bool output=false);
        void step(double dt);
        
        virtual void print(double t_,double *in)=0;     
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
		 
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
        
};

#endif
