#include "rbfAnalysis.hh"

/* constructor function-- simply creates Rycroft rbf and conj grad routine 
so as to not alter the source code */
rbf_test::rbf_test(int n_,int type_) : rbf(n_,type_){
	set_length_scale(5./sqrt(n_));
}

/* initalize points without sorting so we can compare */
void rbf_test::init_unsortedrandom() {
    // Create random positions and resort y positions
    for(int i=0;i<n;i++) {
        px[i]=urand();
        py[i]=urand();
    }
    // Create function values
    for(int i=0;i<n;i++) {
        pf[i]=exp(-2*(px[i]*px[i]+py[i]*py[i]));
        rs[i]=pf[i];
    }
}

/* initialize points calling inherited rycrofts random_init or the local class functions 
including init_hilbert */
void rbf_test::initPoints(int mode, int hilbertM){
	if(mode==0) init_unsortedrandom();
	if(mode==1) init_random(); // chris sort by Y from rbf class
	if(mode==2) init_Hilbert(hilbertM);
}

/* check demo of RBF solver with lapack 
-- this has same functionality as inherited 
rbf::solve_weights_lapack() but with prints*/
void rbf_test::demo_LapackdenseSolver(){
	printf("\n Demo/Test of lapack RBF Solve\n ");
	assemble_matrix();
	printf("Assembled Matrix A\n");
	printMatrix(A, n);
	
	printf("RHS Vector\n");
	printVector(pf, n);
	
   solve_sym_matrix(n,A,rs);
	printf("RBF Solution \n");
	printVector(rs, n);	
}

void rbf_test::demo_JacobiPCG(int bls){
	printf("Demo/Test of JacobiPCG RBF Solve\n ");
   make_table();copy(pf,b);
   printf("Show Jacobi Blocks; bls: %d\n", bls);
	countJacobiP(true, bls);
	
   preconditioning_table(bls);
	printf("RHS Vector\n");
	printVector(pf, n);
	printf("Solve_pre Print Output:\n");
   solve_pre(true);
   copy(x,rs);
	printf("RBF Solution \n");
	printVector(rs, n);	
}
double rbf_test::timeSolve_JacobiPCG(int mode){
 	int l=0;
   double t0=omp_get_wtime();
   double t1;	
   do{
   	make_table();copy(pf,b);
   	preconditioning_table(sqrt(n)); //using sqrt(k) bls
  		solve_pre(false);
  		t1 = omp_get_wtime();
  		l++;
   }while(t1<t0+0.5);
   return (t1-t0)/l;
}

void rbf_test::printMatrix(double* mat, int n_){
	for(int j=0;j<n_;j++){
		for(int i=0;i<n_;i++){
			printf("%g ", mat[j+i*n_]);
		}printf("\n");
	}printf("\n");
}

void rbf_test::printVector(double* vec, int n_){
	for(int j=0;j<n_;j++){
		printf("%g ", vec[j]);
	}printf("\n");
}

/* following are codes to count elements */
int rbf_test::countNonzero(double* mat, int n_){
	int c=0;
	for(int ki=0; ki<n_*n_; ki++){
		if(mat[ki]!=0.0) c+=1;	
	}	
	return c;
}

int rbf_test::countDenseT(){
	assemble_matrix();
	return countNonzero(A,n);
}

int rbf_test::countJacobiP(bool print, int bls){
	//assemble_matrix();
   make_table();
	int lbls = n%bls;
	int fb=n/bls;
	
    // Compute each diagonal block and add up number of elements in blocks
    int k;
    int totalCounter=0;
    int blockcounter=0;
    for(k=0;k<=n-bls;k+=bls) {
        double* Ap = new double[bls*bls];
        fill_matrix_entries(Ap, k, bls);
        int c = countNonzero(Ap, bls);
        if(print){
		     printf("Jacobi Block #%d; # Non-zero: %d \n", blockcounter, c);
		     printMatrix(Ap, bls);
		     blockcounter+=1;}
        totalCounter+=c;
        delete[] Ap;
    }
    
    if(k<n) {
        double* Ap = new double[bls*bls];
        fill_matrix_entries(Ap,k,lbls);
        int c = countNonzero(Ap, bls);
        if(print){
		     printf("Jacobi Block #%d; # Non-zero: %d \n", blockcounter, c);
   	     printMatrix(Ap, bls);
   	     blockcounter+=1;}
        totalCounter+=c;
        delete[] Ap;
    }
	return totalCounter;
}


// vector approach motivated by the sorting solution 
// discussed on https://www.geeksforgeeks.org/sorting-vector-tuple-c-ascending-order/
// I could write it differently but I think their solution is very elegent and makes use
// of the same std::sort function in init_random giving nice continuity in implementation
void rbf_test::init_Hilbert(int orderM){
    // create a vector tuple to hold the coordinates and sort
    std::vector<std::tuple<int, double, double>> v;
    for(int i=0;i<n;i++) {
    	double thisx = urand();
    	double thisy = urand();
		double q = calculateq(orderM, thisx, thisy);
		v.push_back(std::make_tuple(q, thisx, thisy)); 
    }    
    std::sort(v.begin(), v.end()); // sort default to based on first component q
    
    // put sorted coordinates into rbf array
    // instead of the tuple and compute pf 
    for(int i=0;i<n;i++) {
	     px[i]= std::get<1>(v[i]);
        py[i]= std::get<2>(v[i]);
        pf[i]=exp(-2*(px[i]*px[i]+py[i]*py[i]));
        rs[i]=pf[i];
    }
}

// code motivated almost entirely via discussion from site 
// https://en.wikipedia.org/wiki/Hilbert_curve
// which I realized was a copy from the original source
// https://people.sc.fsu.edu/~jburkardt/cpp_src/hilbert_curve/hilbert_curve.cpp
// function returns the hilbert coordinate/cell q which is on N=2^M x N=2^M grid
// M is the order of the hilbert curve and really should be >0
int rbf_test::calculateq(int M, double xval_, double yval_){
	// xval and yval is passed in between -1 and 1 but must be 
	// adjusted to 0 to N
	 int N = pow(2,M);
	 int xval = (int) ((xval_+1.)/2.*N);
	 int yval = (int) ((yval_+1.)/2.*N);
	 int rx, ry, s, d=0;
    for (s=N/2; s>0; s/=2) {
        rx = (xval & s) > 0;
        ry = (yval & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(n, xval, yval, rx, ry);
    }
    return d;
}

void rbf_test::rot (int N, int &xval, int &yval, int rx, int ry){
	int t;
	if ( ry == 0 ){//reflect
		if ( rx == 1 ){
      xval = N - 1 - xval;
      yval = N - 1 - yval;
		}
//  Flip.
    t = xval;
    xval = yval;
    yval = t;
	}
  return;
}

void rbf_test::print_HilbertData(){
	for(int i=0;i<n;i++){
		printf("%d %g %g %g", i, px[i], py[i], pf[i]);
		puts(" ");
    }puts(" ");
}









