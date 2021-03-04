#ifndef SOL_IRK5_HH
#define SOL_IRK5_HH


/** Class for solving an ODE IVP using the Geng  5th order IRK */
class irk5 {
	public:
		/** The total number of degrees of freedom in the ODE system. */
		int dof;
		/** A counter for the number of function evaluations. */
		int fcount;
		/** The current time. */
		double t;
		/** The solution vector. */
		double *q;
				 
		irk5(int dof_);
		virtual ~irk5();
		bool step(double dt);
		bool solve_fixed(double duration, int steps, bool output);
		bool solve_Poincaire(double duration, int steps);
		virtual void init() = 0;
      virtual void ff(double t_,double *in,double *out) = 0;
	   virtual void print(double t_,double *in)=0;     

	private:
		double *dq;
		double *k1;
		double *k2;
		double *k3;
		double *k1b;
		double *k2b;
		double *k3b;
};

#endif
