//LAST I DID: Move main do-while loop before MPI_INIT
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DIM (5)

#define f1(c1,c2,c3,c4,c5)  (50.0 - 0.0*c2 + 1.0*c3 - 0.0*c4 - 0.0*c5)/6.0
#define f2(c1,c2,c3,c4,c5)  (0.0 + 1.0*c1 - 0.0*c3 - 0.0*c4 + 0.0*c5)
#define f3(c1,c2,c3,c4,c5)  (160 - 0.0*c1 + 1.0*c2 - 0.0*c4 - 0.0*c5)/9.0
#define f4(c1,c2,c3,c4,c5)  (  0 - 0.0*c1 + 1.0*c2 + 8.0*c3 + 2.0*c5)/11.0
#define f5(c1,c2,c3,c4,c5)  (  0 + 3.0*c1 + 1.0*c2 - 0.0*c3 - 0.0*c4)/4.0

double a[DIM][DIM] = {{6.0, 0.0, -1.0, 0.0, 0.0}, {-3.0, 3.0, 0.0, 0.0, 0.0}, {0.0, -1.0, 9.0, 0.0, 0.0}, {0.0, -1.0, -8.0, 11.0, -2.0}, {-3.0, -1.0, 0.0, 0.0, 4.0} };
double b[DIM] = {50.0, 0.0, 160.0, 0.0, 0.0};
double sol[DIM]={11.509, 11.509, 19.057, 16.998, 11.509};
//double local_x[DIM]={0.0, 0.0, 0.0, 0.0, 0.0};
double local_x[1];
int n = DIM;
void vector_print(int nr, double *x);

int main(void) {
  double x0=0.0, y0=0.0, z0=0.0, a0=0.0, b0=0.0;
  double x1=0.0, y1=0.0, z1=0.0, a1=0.0, b1=0.0;
  double e1=0.0, e2=0.0, e3=0.0, e4=0.0, e5=0.0, e=0.0;
  double local_m=0.0;
  int m=0, rc=0, count=1;
  int my_rank=0, comm_sz=0, size=0;
  int local_start=0, local_end=0;
  const int root = 0;

  MPI_Comm comm;
  MPI_Init(NULL, NULL);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &my_rank);
  m=5;
  local_m = (double)m/(double)comm_sz;
  size = local_m*(double)n;
  double local_a[DIM];
  double local_rhs[n];
  double x[DIM];
  double RHS[DIM];


  if(my_rank == 0){
    printf("Enter tolerable error:\n");
    rc = scanf("%lf", &e);
    //for(int i =0; i < size; i++)
    // local_rhs[i]=0.0;
    if(rc < 0) 
      perror("input");

    do
    {

      /* Calculation */
      x1 = f1(x0,y0,z0,a0,b0);
      y1 = f2(x1,y0,z0,a0,b0);
      z1 = f3(x1,y1,z0,a0,b0);
      a1 = f4(x1,y1,z1,a0,b0);
      b1 = f5(x1,y1,z1,a1,b0);

      /* Error */
      e1 = fabs(x0-x1);
      e2 = fabs(y0-y1);
      e3 = fabs(z0-z1);
      e4 = fabs(a0-a1);
      e5 = fabs(b0-b1);

      count++;

      /* Set value for next iteration */
      x0 = x1;
      y0 = y1;
      z0 = z1;
      a0 = a1;
      b0 = b1;

    } while( (e1>e) || (e2>e) || (e3>e) || (e4>e) || (e5>e) );

    printf("\nGSIT Solution: x=%0.4f, y=%0.4f, z = %0.4f  a = %0.4f  b = %0.4f\n\n",x1,y1,z1,a1,b1);
    printf("Math Tool Solution:\n");
    vector_print(n, sol);

    x[0] = x1;
    x[1] = y1;
    x[2] = z1;
    x[3] = a1;
    x[4] = b1;
    // Additional verification steps to assess error in answer.
  }

  MPI_Barrier(comm);
  //Scatter matrix A by rowns
  MPI_Scatter(a, size, MPI_DOUBLE, local_a, size, MPI_DOUBLE, root, comm);
  //Scatter matrix C by columns
  MPI_Scatter(x, 1, MPI_DOUBLE, local_x, 1, MPI_DOUBLE, root, comm);

  for(int idx=0; idx < DIM; idx++)
  {
    ////        a = B             local_a = Bpart
    //local_x == Cpart C = x
    local_rhs[idx] = local_x[0] * local_a[idx];
  }

  MPI_Reduce(local_rhs, RHS, DIM ,MPI_DOUBLE, MPI_SUM, 0, comm);

  if (my_rank==0){
    //Compare original RHS "b" to computed RHS
    printf("Computed RHS is:\n");
    vector_print(n, RHS);

    printf("Original is:\n");
    vector_print(n, b);
  }




  //  verify(a, local_a, local_rhs, b, x, local_m, my_rank, comm, local_start, local_end);
  //MPI_Gather(&local_rhs, size, MPI_DOUBLE, &soln, size, MPI_DOUBLE, 0, comm);

  MPI_Barrier(comm);
  MPI_Finalize();

  return 0;
}

void vector_print(int nr, double *x) {
  int row_idx;

  for (row_idx = 0; row_idx < nr; row_idx++)
  {
    printf ("index = %d, %9.4f  \n",row_idx, x[row_idx]);
  }

  printf("\n");  // Insert a new line at the end
}
