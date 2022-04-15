// A sequential version of matrix multiplacatino with matrix and matrix multiplaction with vector
// This solution adapts code from Dr. Sam Siewert's 'vmult.c' and 'gewpp.c' programs
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ICHAR 80  // Length of array holding description of the problem
#define MAX_DIM 1000 // Length of array holding description of the problem

// Function Prototypes
void  matrix_print (int nr, int nc, double **A);
void  vector_print (int nr, double *x);
//void seqvmult(int n, double **A, double *x);
void seqvmult(int n, double A[][MAX_DIM], double *x);
void seqmmult(int n, double **A, double **B);
void instantiate(int n, double **default_mat, double *defaul_vect, double **C);
void verify(double **a, double *b, double *x, int n);
//void instantiate(int n, double **default_mat, double **alt_mat, double *defaul_vect, double **C);
//void  matrix_print (int nr, int nc, double **a,  double A[][MAX_DIM]) {
//void  second_matrix_print (int nr, int nc, double A[][MAX_DIM]);
//void  second_matrix_print (int nr, int nc, int temp, double A[][MAX_DIM]);
//void  vector_print (int nr, double x[MAX_DIM]);
//void seqvmult(int n, double A[][MAX_DIM], double x[MAX_DIM]);
//void seqmmult(int n, double A[][MAX_DIM], double B[][MAX_DIM]);

double **default_mat;
//double **alt_mat;
double *default_vect;
double **C;
//double default_mat[MAX_DIM][MAX_DIM];
double alt_mat[MAX_DIM][MAX_DIM];
//double default_vect[MAX_DIM];
//double C[MAX_DIM][MAX_DIM];

int debug=0; // debug is normally off, but turned on for simple examples

int main (int argc, char *argv[]) {
  int count=0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double  **a, *b, *x;        // a=coefficients, b=RHS, x=solution vector
    double  **asave, *bsave;    // asave=coefficients original, bsave=RHS original
    char   desc[ICHAR];         // file header description 
    int    row_idx, col_jdx, n; // indexes into arrays and dimention of linear system
    FILE   *finput;             // file pointer to description, dimension, coefficients, and RHS

  // Timing declarations
  struct timespec start, end;
  double fstart=0.0, fend=0.0;

    if(argc > 1)
    {
       printf("Using custom input file %s: argc=%d, argv[0]=%s, argv[1]=%s\n", argv[0], argc, argv[0], argv[1]); 
       finput = fopen(argv[1],"r");
    }
    else
    {
        printf("Using gauss.dat DEFAULT example input file\n");
        finput = fopen("gauss.dat","r");
    }

    if (finput == NULL) 
    { 
    	printf("Data file gauss.dat not found\n");
    	return(-1);
    }

    // Get a one line description of the matrix problem and
    // then the dimension, n, of the system A[n][n] and b[n]  
    fgets(desc, ICHAR , finput);      
    fscanf(finput, "%d", &n);
    //instantiate(n, default_mat, alt_mat, default_vect, C);
    instantiate(n, default_mat, default_vect, C);
    printf("%s", desc);  
    printf("\nDimension of matrix = %d\n\n", n);


    // Dynamic allocation of the arrays and vectors
    //
    // Allocate n pointers that point to each row in a 2D array
    //
    a = calloc(n, sizeof(double *));
    asave = calloc(n, sizeof(double *));

    // Now allocate n rows with n columns each, one row at a time
    //
    // A 2D array could alternatively be used, but this approach allows us
    // to use C pointers, for better or worse, which are efficient.
    //
    // possible loop pragma here
    //
    for(row_idx = 0; row_idx < n; ++row_idx) 
    { 
        // Allocate n column values for this row
    	a[row_idx] =  calloc(n, sizeof(double));
    	asave[row_idx] =  calloc(n, sizeof(double));
    }

    // Allocate space for the RHS vector of dimension n
    b = calloc(n, sizeof(double));
    bsave = calloc(n, sizeof(double));

    // Allocate space for the solution vector (unknowns) of dimension n
    x = calloc(n, sizeof(double));

    printf("Memory allocation done\n");

    // Read the elements of coefficient array A
    for (row_idx=0; row_idx < n; row_idx++)
    {
    	for (col_jdx=0; col_jdx < n; col_jdx++) 
        {
    		fscanf(finput,"%lf ",&a[row_idx][col_jdx]);
            asave[row_idx][col_jdx]=a[row_idx][col_jdx];
    	}
    }

    printf("Coefficient array read done\n");

    // Read the elements of RHS vector "b"
    for (row_idx=0; row_idx < n; row_idx++)
    {
       fscanf(finput,"%lf ",&b[row_idx]);
       // RHS b gets rearranged by pivoting, so save original
       bsave[row_idx]=b[row_idx];
    }
       
    printf("RHS vector read done\n");

    fclose(finput); // Close the input file

    // Now print out the problem to be solved in vector matrix form for n equations, n unknowns
    printf("\nMatrices read from input file\n");

    printf("\nCoefficient Matrix A\n\n");
    matrix_print(n, n, a);

    printf("\nRHS Vector b\n\n");
    vector_print(n, bsave);


    // Sequential tests
    clock_gettime(CLOCK_MONOTONIC, &start);
    seqvmult(n, alt_mat, bsave);
    clock_gettime(CLOCK_MONOTONIC, &end);
    fstart=start.tv_sec + (start.tv_nsec / 1000000000.0);
    fend=end.tv_sec + (end.tv_nsec / 1000000000.0);
    printf("SEQUENTIAL vmult %d took %lf seconds\n", n, (fend-fstart));

    // Sequential tests
    clock_gettime(CLOCK_MONOTONIC, &start);
//    seqmmult(n, default_mat, default_mat);
    clock_gettime(CLOCK_MONOTONIC, &end);
    fstart=start.tv_sec + (start.tv_nsec / 1000000000.0);
    fend=end.tv_sec + (end.tv_nsec / 1000000000.0);
    printf("SEQUENTIAL mmult %d took %lf seconds\n", n, (fend-fstart));

  printf("\nSolution x\n\n");
  vector_print(n, bsave);

  // Multiply solution "x" by matrix "a" to verify we get RHS "bsave"
	verify(asave, bsave, b, n); 

  return(0);
}//main 

//void seqvmult(int n, double A[][MAX_DIM], double x[MAX_DIM]) {
//void seqvmult(int n, double **A, double *x){
void seqvmult(int n, double A[][MAX_DIM], double *x){
  int row_idx, col_jdx;
  //double rhs[n];
  double *rhs;

  rhs = calloc(n, sizeof(double));
  for (row_idx=0; row_idx < n; ++row_idx) {
    rhs[row_idx] = 0.0;
    for (col_jdx=0; col_jdx < n; ++col_jdx) {
      ///////////////////////////////////////DOES NOT LIKE INDEXING IN A array
      //rhs[row_idx] += A[row_idx][col_jdx] * x[col_jdx];
      rhs[row_idx] += A[row_idx][col_jdx] * x[col_jdx];
      rhs[row_idx] += A[row_idx][col_jdx] * x[col_jdx];
      printf("X at %d = %lf\n",col_jdx, x[col_jdx]);
    }
  }

  if(debug) {
    printf("Computed RHS is:\n");
   // vector_print(n, rhs);
  }
}//seqvmult

//void seqmmult(int n, double A[][MAX_DIM], double B[][MAX_DIM]) {
void seqmmult(int n, double **A, double **B){
  int row_idx, col_jdx, coeff_idx;

  for (row_idx=0; row_idx < n; ++row_idx) {
    for (col_jdx=0; col_jdx < n; ++col_jdx) {
      if(debug) printf("C[%d][%d]:\n", row_idx, col_jdx);

      for(coeff_idx=0; coeff_idx < n; ++coeff_idx) {
        C[row_idx][col_jdx] += A[row_idx][coeff_idx] * B[coeff_idx][col_jdx];
        if(debug) printf("A[%d][%d]=%lf, b[%d][%d]=%lf\n", row_idx, coeff_idx, A[row_idx][coeff_idx], coeff_idx, col_jdx, B[coeff_idx][col_jdx]);
      }
    }
  }

//  if(debug) {
    printf("Computed C is:\n");
    printf("A=\n"); matrix_print(n, n, A); printf("\n");
    printf("B=\n"); matrix_print(n, n, B); printf("\n");
    printf("C=\n"); matrix_print(n, n, C); printf("\n");
 // }
}//seqmmult

//void  matrix_print (int nr, int nc, double **a,  double A[][MAX_DIM]) {
void  matrix_print (int nr, int nc, double **A) {
  int row_idx=0, col_jdx=0;

  for (row_idx = 0; row_idx < nr; row_idx++) {
    for (col_jdx = 0; col_jdx < nc; col_jdx++) {
      printf ("%9.4f  ", A[row_idx][col_jdx]);
    }

    printf("\n"); // Insert a new line at end of each row
  }
}//matrix_print

//void  vector_print (int nr, double x[MAX_DIM]) {
void  vector_print (int nr, double *x) {
  int row_idx;

  for (row_idx = 0; row_idx < nr; row_idx++) {
    printf ("%9.4f  \n", x[row_idx]);
  }

  printf("\n");  // Insert a new line at the end
}//vector_print

//void instantiate(int n, double **default_mat, double **alt_mat, double *defaul_vect, double **C){
void instantiate(int n, double **default_mat, double *defaul_vect, double **C){
    default_mat = calloc(n, sizeof(double *));
   // alt_mat = calloc(n, sizeof(double *));
    C = calloc(n, sizeof(double *));
    default_vect = calloc(n, sizeof(double));

    for(int row_idx = 0; row_idx < n; ++row_idx) 
    { 
        // Allocate n column values for this row
    	default_mat[row_idx] =  calloc(n, sizeof(double));
   // 	alt_mat[row_idx] =  calloc(n, sizeof(double));
    	C[row_idx] =  calloc(n, sizeof(double));
    }
}



void verify(double **a, double *b, double *x, int n){
  int row_idx, col_jdx;
  double rhs[n];
  
  for(row_idx=0; row_idx < n; ++row_idx){
    rhs[row_idx]=0.0;

    for(col_jdx=0; col_jdx < n; ++col_jdx){
      rhs[row_idx] += a[row_idx][col_jdx] * x[col_jdx];
      printf("\nRHS value at %d is equal to %lf\n", row_idx, a[row_idx][col_jdx]);
    }
  }

  printf("Computed RHS is: \n");
  vector_print(n, rhs);

  printf("Original RHS is: \n");
  vector_print(n, b);
}
