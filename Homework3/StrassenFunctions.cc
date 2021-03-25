#include "StrassenFunctions.hh"


double* genRandomMat(int N){
	double* mat_ = new double[N*N];	
	for(int j=0; j<N*N; j++){
		mat_[j] = doubleRandGen(10.);
		//mat_[j] = j;
	}
	return mat_;
}


void printMat(int N, double* mat){
	for(int j=0; j<N; j++){
		for(int i=0; i<N; i++){
			printf(" %g ", mat[i+j*N]);			
		};puts("");
	};puts("");
}


double* mat_mul(int N, double* A, double* B){
	double* C = new double[N*N]; // output matrix
	for(int j=0; j<N; j++){// row enumerate
		for(int i=0; i<N; i++){ // col enumerate
			double holdsum=0;
			for(int k=0; k<N; k++){
				holdsum += A[j*N + k]*B[i+k*N];
			}
			C[i+j*N] = holdsum;
		}	
	}
	return C;
}
double* mat_mul_menManage(int N, double* A, double* B){
	double* C = new double[N*N]; // output matrix
	
	// I found with a quick and non-thorough test that this 
	// is a low cost operation that can speed up potentially
	double* BTranspose = new double[N*N];
	for(int j=0; j<N; j++){
		for(int i=0; i<N; i++){
			BTranspose[i+j*N] = B[j+i*N];
		}
	}
	
	for(int j=0; j<N; j++){// row enumerate
		double Arow[N]; 
		std::copy(A+j*N, A+(j+1)*N, Arow);
		double Crow[N];
		for(int i=0; i<N; i++){// col enumerate
			double Bcol[N]; 
			std::copy(BTranspose+i*N, BTranspose+(i+1)*N, Bcol);
			Crow[i]=0.0;
			for(int k=0; k<N; k++) Crow[i] += Arow[k]*Bcol[k];
		}
		std::copy(Crow, Crow+N, C+j*N);	
	}
	delete BTranspose;
	return C;
}


void matAdd(int N, double* out, double* A, double* B, double alpha){
	for(int i=0; i<(N*N); i++){
		out[i] = A[i] + alpha*B[i];
	}
}
void matAdd_memManage(int N, double* out, double* A, double* B, double alpha){
	for(int j=0; j<N; j++){
		double Achunk[N], Bchunk[N], outChunk[N];
		std::copy(A+j*N, A+(j+1)*N, Achunk);
		std::copy(B+j*N, B+(j+1)*N, Bchunk);		
		for(int i=0;i<N;i++) outChunk[i] = Achunk[i] + alpha*Bchunk[i];
		std::copy(outChunk, outChunk+N, out+j*N);
	}
}


void splitStrassen(int N, double* in, double* zerozero, double* zeroone, double* onezero, double* oneone){
	for(int j=0; j<N/2; j++){//iterate over row chunks top half
		for(int i=0; i<N/2; i++){//iterate over columns first half
			zerozero[i+j*(N/2)] = in[i+j*N];
		}
		for(int i=N/2; i<N; i++){//iterate over columns second half
			zeroone[(i-N/2)+j*(N/2)] = in[i+j*N];
		}
	}
	
	for(int j=N/2; j<N; j++){//iterate over row chunks bottom half
		for(int i=0; i<N/2; i++){//iterate over columns first half
			onezero[i+(j-N/2)*(N/2)] = in[i+j*N];
		}
		for(int i=N/2; i<N; i++){//iterate over columns second half
			oneone[(i-N/2)+(j-N/2)*(N/2)] = in[i+j*N];
		}
	}
}

void recombineStrassen(int N, double* out, double* Czerozero, double* Czeroone, double* Conezero, double*Coneone){
	// This is nearly the same code as splitStrassen but with in and out swapped
	for(int j=0; j<N/2; j++){//iterate over row chunks top half
		for(int i=0; i<N/2; i++){//iterate over columns first half
			out[i+j*N] = Czerozero[i+j*(N/2)];
		}
		for(int i=N/2; i<N; i++){//iterate over columns second half
			out[i+j*N] = Czeroone[(i-N/2)+j*(N/2)];
		}
	}
	
	for(int j=N/2; j<N; j++){//iterate over row chunks bottom half
		for(int i=0; i<N/2; i++){//iterate over columns first half
			out[i+j*N] = Conezero[i+(j-N/2)*(N/2)];
		}
		for(int i=N/2; i<N; i++){//iterate over columns second half
			out[i+j*N] = Coneone[(i-N/2)+(j-N/2)*(N/2)];
		}
	}
}

/* Strassen Function with No Recursion*/
void writeStrassen(int N, double* out, double* A, double* B){
	// Get Strassen Sub matrices
	int Nhalf = N/2;
	double* Azerozero = new double[Nhalf*Nhalf];
	double* Azeroone = new double[Nhalf*Nhalf];
	double* Aonezero = new double[Nhalf*Nhalf];
	double* Aoneone = new double[Nhalf*Nhalf];
	splitStrassen(N, A, Azerozero, Azeroone, Aonezero, Aoneone);
	double* Bzerozero = new double[Nhalf*Nhalf];
	double* Bzeroone = new double[Nhalf*Nhalf];
	double* Bonezero = new double[Nhalf*Nhalf];
	double* Boneone = new double[Nhalf*Nhalf];
	splitStrassen(N, B, Bzerozero, Bzeroone, Bonezero, Boneone);
		
	// Do the Strassen Sub computation
	double* Q0 = new double[Nhalf*Nhalf];
	double* Q1 = new double[Nhalf*Nhalf];
	double* Q2 = new double[Nhalf*Nhalf];
	double* Q3 = new double[Nhalf*Nhalf];
	double* Q4 = new double[Nhalf*Nhalf];
	double* Q5 = new double[Nhalf*Nhalf];
	double* Q6 = new double[Nhalf*Nhalf];
	double* temp1 = new double[Nhalf*Nhalf];
	double* temp2 = new double[Nhalf*Nhalf];
	
	matAdd(Nhalf, temp1, Azerozero, Aoneone, 1.0);
	matAdd(Nhalf, temp2, Bzerozero, Boneone, 1.0);
	
	Q0 = mat_mul_menManage(Nhalf, temp1, temp2);
	
	matAdd(Nhalf, temp1, Aonezero, Aoneone, 1.0);
	Q1 = mat_mul_menManage(Nhalf,temp1, Bzerozero);
	
	matAdd(Nhalf, temp1, Bzeroone, Boneone, -1.0);
	Q2 = mat_mul_menManage(Nhalf, Azerozero, temp1);
	
	matAdd(Nhalf, temp1, Bonezero, Bzerozero, -1.0);
	Q3 = mat_mul_menManage(Nhalf, Aoneone, temp1);
	
	matAdd(Nhalf, temp1, Azerozero, Azeroone, 1.0);
	Q4 = mat_mul_menManage(Nhalf, temp1, Boneone);
	
	matAdd(Nhalf, temp1, Aonezero, Azerozero, -1.0);
	matAdd(Nhalf, temp2, Bzerozero, Bzeroone, 1.0);
	Q5 = mat_mul_menManage(Nhalf, temp1, temp2);

	matAdd(Nhalf, temp1, Azeroone, Aoneone, -1.0);
	matAdd(Nhalf, temp2, Bonezero, Boneone, 1.0);
	Q6 = mat_mul_menManage(Nhalf, temp1, temp2);	
	
	
	// output stored in A submatrix then stiched into out matrix
	matAdd(Nhalf, temp1, Q0, Q3, 1.0);
	matAdd(Nhalf, temp2, Q6, Q4, -1.0);
	matAdd(Nhalf, Azerozero, temp1, temp2, 1.0);
	matAdd(Nhalf, Aonezero, Q1, Q3, 1.0);
	matAdd(Nhalf, Azeroone, Q2, Q4, 1.0);
	matAdd(Nhalf, temp1, Q0, Q2, 1.0);
	matAdd(Nhalf, temp2, Q5, Q1, -1.0);
	matAdd(Nhalf, Aoneone, temp1, temp2, 1.0);
	
	recombineStrassen(N, out, Azerozero, Azeroone, Aonezero, Aoneone);
	
	delete Boneone;
	delete Bonezero;
	delete Bzeroone;
	delete Bzerozero;
	
	delete Aoneone;
	delete Aonezero;
	delete Azerozero;
	delete Azeroone;
	
	delete temp1;
	delete temp2;		
	
	delete Q0; 
	delete Q1; 
	delete Q2; 
	delete Q3;
	delete Q4;
	delete Q5;
	delete Q6;
}
/* Strassen With Recursion*/
double* writeStrassen2(int N, double* A, double* B){
	// Get Strassen Sub matrices
	int Nhalf = N/2;
	double* Azerozero = new double[Nhalf*Nhalf];
	double* Azeroone = new double[Nhalf*Nhalf];
	double* Aonezero = new double[Nhalf*Nhalf];
	double* Aoneone = new double[Nhalf*Nhalf];
	splitStrassen(N, A, Azerozero, Azeroone, Aonezero, Aoneone);
	double* Bzerozero = new double[Nhalf*Nhalf];
	double* Bzeroone = new double[Nhalf*Nhalf];
	double* Bonezero = new double[Nhalf*Nhalf];
	double* Boneone = new double[Nhalf*Nhalf];
	splitStrassen(N, B, Bzerozero, Bzeroone, Bonezero, Boneone);
		
	// Do the Strassen Sub computation
	double* Q0 = new double[Nhalf*Nhalf];
	double* Q1 = new double[Nhalf*Nhalf];
	double* Q2 = new double[Nhalf*Nhalf];
	double* Q3 = new double[Nhalf*Nhalf];
	double* Q4 = new double[Nhalf*Nhalf];
	double* Q5 = new double[Nhalf*Nhalf];
	double* Q6 = new double[Nhalf*Nhalf];
	double* temp1 = new double[Nhalf*Nhalf];
	double* temp2 = new double[Nhalf*Nhalf];
	
	if(Nhalf<=256){
		// Do the computation explicitly
		matAdd(Nhalf, temp1, Azerozero, Aoneone, 1.0);
		matAdd(Nhalf, temp2, Bzerozero, Boneone, 1.0);
		
		Q0 = mat_mul_menManage(Nhalf, temp1, temp2);
		
		matAdd(Nhalf, temp1, Aonezero, Aoneone, 1.0);
		Q1 = mat_mul_menManage(Nhalf,temp1, Bzerozero);
		
		matAdd(Nhalf, temp1, Bzeroone, Boneone, -1.0);
		Q2 = mat_mul_menManage(Nhalf, Azerozero, temp1);
		
		matAdd(Nhalf, temp1, Bonezero, Bzerozero, -1.0);
		Q3 = mat_mul_menManage(Nhalf, Aoneone, temp1);
		
		matAdd(Nhalf, temp1, Azerozero, Azeroone, 1.0);
		Q4 = mat_mul_menManage(Nhalf, temp1, Boneone);
		
		matAdd(Nhalf, temp1, Aonezero, Azerozero, -1.0);
		matAdd(Nhalf, temp2, Bzerozero, Bzeroone, 1.0);
		Q5 = mat_mul_menManage(Nhalf, temp1, temp2);

		matAdd(Nhalf, temp1, Azeroone, Aoneone, -1.0);
		matAdd(Nhalf, temp2, Bonezero, Boneone, 1.0);
		Q6 = mat_mul_menManage(Nhalf, temp1, temp2);	
		
		// output stored in A submatrix then stiched into out matrix
		matAdd(Nhalf, temp1, Q0, Q3, 1.0);
		matAdd(Nhalf, temp2, Q6, Q4, -1.0);
		matAdd(Nhalf, Azerozero, temp1, temp2, 1.0);
		matAdd(Nhalf, Aonezero, Q1, Q3, 1.0);
		matAdd(Nhalf, Azeroone, Q2, Q4, 1.0);
		matAdd(Nhalf, temp1, Q0, Q2, 1.0);
		matAdd(Nhalf, temp2, Q5, Q1, -1.0);
		matAdd(Nhalf, Aoneone, temp1, temp2, 1.0);
		
		double* out = new double[N*N];
		recombineStrassen(N, out, Azerozero, Azeroone, Aonezero, Aoneone);
		delete Boneone;
		delete Bonezero;
		delete Bzeroone;
		delete Bzerozero;
		
		delete Aoneone;
		delete Aonezero;
		delete Azerozero;
		delete Azeroone;		
		delete temp1;
		delete temp2;		
		delete Q0; 
		delete Q1; 
		delete Q2; 
		delete Q3;
		delete Q4;
		delete Q5;
		delete Q6;
		return out;
	}else{
		// Do the computation with recursive Strassen calls
		matAdd(Nhalf, temp1, Azerozero, Aoneone, 1.0);
		matAdd(Nhalf, temp2, Bzerozero, Boneone, 1.0);
		Q0 = writeStrassen2(Nhalf, temp1, temp2);
		
		matAdd(Nhalf, temp1, Aonezero, Aoneone, 1.0);
		Q1 = writeStrassen2(Nhalf,temp1, Bzerozero);
		
		matAdd(Nhalf, temp1, Bzeroone, Boneone, -1.0);
		Q2 = writeStrassen2(Nhalf, Azerozero, temp1);
		
		matAdd(Nhalf, temp1, Bonezero, Bzerozero, -1.0);
		Q3 = writeStrassen2(Nhalf, Aoneone, temp1);
		
		matAdd(Nhalf, temp1, Azerozero, Azeroone, 1.0);
		Q4 = writeStrassen2(Nhalf, temp1, Boneone);
		
		matAdd(Nhalf, temp1, Aonezero, Azerozero, -1.0);
		matAdd(Nhalf, temp2, Bzerozero, Bzeroone, 1.0);
		Q5 = writeStrassen2(Nhalf, temp1, temp2);

		matAdd(Nhalf, temp1, Azeroone, Aoneone, -1.0);
		matAdd(Nhalf, temp2, Bonezero, Boneone, 1.0);
		Q6 = writeStrassen2(Nhalf, temp1, temp2);	
		
		// output stored in A submatrix then stiched into out matrix
		matAdd(Nhalf, temp1, Q0, Q3, 1.0);
		matAdd(Nhalf, temp2, Q6, Q4, -1.0);
		matAdd(Nhalf, Azerozero, temp1, temp2, 1.0);
		matAdd(Nhalf, Aonezero, Q1, Q3, 1.0);
		matAdd(Nhalf, Azeroone, Q2, Q4, 1.0);
		matAdd(Nhalf, temp1, Q0, Q2, 1.0);
		matAdd(Nhalf, temp2, Q5, Q1, -1.0);
		matAdd(Nhalf, Aoneone, temp1, temp2, 1.0);
		
		double* out = new double[N*N];
		recombineStrassen(N, out, Azerozero, Azeroone, Aonezero, Aoneone);
		delete Boneone;
		delete Bonezero;
		delete Bzeroone;
		delete Bzerozero;
		
		delete Aoneone;
		delete Aonezero;
		delete Azerozero;
		delete Azeroone;		
		delete temp1;
		delete temp2;		
		delete Q0; 
		delete Q1; 
		delete Q2; 
		delete Q3;
		delete Q4;
		delete Q5;
		delete Q6;
		return out;
	}
}


double* BlasMultiply(int N, double* A, double* B){
  // Test corresponding BLAS routine. Due to C++/Fortran calling
  // conventions, all entries are passed as pointers, and there is an
  // underscore after the name.
	char trans='n';
	double alpha=1.0, beta=0.0;
	double* C = new double[N*N];
	
	dgemm_(&trans,&trans,&N,&N,&N,&alpha,A,&N,B,&N,&beta,C,&N);
	return C;
}





















