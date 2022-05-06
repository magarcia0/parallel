// An omp version of Matrix x Matrix and Matrix x Vector multiplication
// This solution adapts code from Dr. Sam Siewert's 'vmult.c' 
// and 'gewpp.c' programs
//
// Date: 
// 4/15/2022
//
// Author:
// Marco A. Garcia
// 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define ICHAR 80  // Length of array holding description of the problem
#define MAX_DIM 1000

// Function Prototypes
void matrix_print (int nr, int nc, double **A);
void vector_print (int nr, double *x);
void ompvmult(int n, double **A, double *x);
void ompmmult(int n, double **A, double **B);

double **default_mat;
double *default_vect;
double **C;
double alt_mat[MAX_DIM][MAX_DIM];

int debug=0; // debug is normally off, but turned on for simple examples
int thread_count=8; // default thread count of 4. Can be changed if desired


int main (int argc, char *argv[]) {

  if(argc < 3){
    printf("Usage: ./serial_matrixmult <file.dat> <'matrix' or 'vector'>\n");
    exit(-1);
  }

  int count=0;
  double  **a, **mat_two, *b, *x;        // a=coefficients, b=RHS, x=solution vector, mat_two=matrix two
  double  **asave, *bsave;    // asave=coefficients original, bsave=RHS original
  char   desc[ICHAR];         // file header description 
  int    row_idx, col_jdx, n; // indexes into arrays and dimention of linear system
  FILE   *finput;             // file pointer to description, dimension, coefficients, and RHS
  int ret;

  // Timing declarations
  struct timespec start, end;
  double fstart=0.0, fend=0.0;
  char *mult_type;


  if(argc > 1)
  {
    finput = fopen(argv[1],"r");
    mult_type=argv[2];

    char* vect = "vector\0";
    char* mat = "matrix\0";
    ret = strncmp(mult_type, vect, 6);
    if(ret!=0){
      ret = strncmp(mult_type, mat, 6);
      if(ret!=0){
        printf("ERROR: Must choose either 'vector' or 'matrix' for third arg\n");
        exit(-1);
      }
    }

    printf("\nUser chose matrix mult with --->%s<---\n", mult_type);
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

  printf("%s", desc);  
  printf("\nDimension of matrix = %d\n\n", n);


  // Dynamic allocation of the arrays and vectors
  a = calloc(n, sizeof(double *));
  asave = calloc(n, sizeof(double *));
  mat_two = calloc(n, sizeof(double *));
  C = calloc(n, sizeof(double *));

  for(row_idx = 0; row_idx < n; ++row_idx) 
  { 
    // Allocate n column values for this row
    a[row_idx] =  calloc(n, sizeof(double));
    asave[row_idx] =  calloc(n, sizeof(double));
    mat_two[row_idx] =  calloc(n, sizeof(double));
    C[row_idx] =  calloc(n, sizeof(double));
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

  if( (ret=strncmp(argv[2], "vector", 6) == 0) ){
    // Read the elements of RHS vector "b"
    for (row_idx=0; row_idx < n; row_idx++){
      fscanf(finput,"%lf ",&b[row_idx]);
      // RHS b gets rearranged by pivoting, so save original
      bsave[row_idx]=b[row_idx];
    }

    printf("RHS vector read done\n");

  }
  else {
    // Read the elements of second matrix "b"
    for (row_idx=0; row_idx < n; row_idx++)
    {
      for (col_jdx=0; col_jdx < n; col_jdx++) 
      {
        fscanf(finput,"%lf ",&mat_two[row_idx][col_jdx]);
      }
    }

    printf("Second matrix read done\n");
  }

  fclose(finput);

  // Now print out the problem to be solved in vector matrix form for n equations, n unknowns
  printf("\nMatrices read from input file\n");

  if( (ret=strncmp(argv[2], "vector", 6) == 0) ){
    printf("Coefficient Matrix A\n\n");
    matrix_print(n, n, a);
    printf("\nRHS Vector b\n\n");
    vector_print(n, b);
  }
  else{
    printf("First Matrix A\n\n");
    matrix_print(n, n, a);
    printf("Second Matrix B\n\n");
    matrix_print(n, n, mat_two);
  }


  if( (ret=strncmp(argv[2], "vector", 6) == 0) ){
    // OMP tests
    clock_gettime(CLOCK_MONOTONIC, &start);
    ompvmult(n, a, bsave);
    clock_gettime(CLOCK_MONOTONIC, &end);
    fstart=start.tv_sec + (start.tv_nsec / 1000000000.0);
    fend=end.tv_sec + (end.tv_nsec / 1000000000.0);
    printf("PARALLEL OpenMP mmult %d took %lf seconds, with %d workers\n", n, (fend-fstart), thread_count);
  }
  else if( (ret=strncmp(argv[2],"matrix", 6) == 0) ){
    // OMP tests
    clock_gettime(CLOCK_MONOTONIC, &start);
    ompmmult(n, a, mat_two);
    clock_gettime(CLOCK_MONOTONIC, &end);
    fstart=start.tv_sec + (start.tv_nsec / 1000000000.0);
    fend=end.tv_sec + (end.tv_nsec / 1000000000.0);
    printf("PARALLEL OpenMP mmult %d took %lf seconds, with %d workers\n", n, (fend-fstart), thread_count);
  }

  return(0);
}//main 


void ompvmult(int n, double **A, double *x){
  int row_idx, col_jdx;
  double *rhs;
  rhs = calloc(n, sizeof(double));

#pragma omp parallel for num_threads(thread_count) private(row_idx, col_jdx) shared(n)
  for (row_idx=0; row_idx < n; ++row_idx) {
    rhs[row_idx] = 0.0;
    for (col_jdx=0; col_jdx < n; ++col_jdx) {
      rhs[row_idx] += A[row_idx][col_jdx] * x[col_jdx];
    }
  }

  printf("Computed Answer is:\n");
  vector_print(n, rhs);
}//ompvmult


void ompmmult(int n, double **A, double **B){
  int row_idx, col_jdx, coeff_idx;

#pragma omp parallel for num_threads(thread_count) private(row_idx, col_jdx, coeff_idx) shared(n)
  for (row_idx=0; row_idx < n; ++row_idx) {
    for (col_jdx=0; col_jdx < n; ++col_jdx) {
      if(debug) printf("C[%d][%d]:\n", row_idx, col_jdx);

      for(coeff_idx=0; coeff_idx < n; ++coeff_idx) {
        C[row_idx][col_jdx] += A[row_idx][coeff_idx] * B[coeff_idx][col_jdx];
        if(debug) printf("A[%d][%d]=%lf, b[%d][%d]=%lf\n", row_idx, coeff_idx, A[row_idx][coeff_idx], coeff_idx, col_jdx, B[coeff_idx][col_jdx]);
      }
    }
  }

  printf("C=\n"); matrix_print(n, n, C); printf("\n");
}//ompmmult


void  matrix_print (int nr, int nc, double **A) {
  int row_idx=0, col_jdx=0;

  for (row_idx = 0; row_idx < nr; row_idx++) {
    for (col_jdx = 0; col_jdx < nc; col_jdx++) {
      printf ("%9.4f  ", A[row_idx][col_jdx]);
    }

    printf("\n");
  }
}//matrix_print


void  vector_print (int nr, double *x) {
  int row_idx;

  for (row_idx = 0; row_idx < nr; row_idx++) {
    printf("%9.4f  \n", x[row_idx]);
  }

  printf("\n");
}//vector_print
