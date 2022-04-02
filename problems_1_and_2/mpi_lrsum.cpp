//This program adapts from L Muller's program
// 'Reimann sum in C++' from the video linked below
// https://www.youtube.com/watch?v=NWCnABJHqIQ&ab_channel=LMuller
#include <iostream>
#include <cmath>
#include <mpi.h>

using std::cout;
using std::endl;

//functon prototypes
double area(double, double, double);
void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p,double* n_p);

int main()
{
  int my_rank=0;
  int comm_sz=0;
  double a=0.0;
  double b=0.0; 
  double n=0.0;
  double integrated_area=0.0;
  double total_int_area=0.0;
  double local_a=0.0;
  double local_b=0.0;
  double local_n=0.0;

  /* Let the system do what it needs to start up MPI */
  MPI_Init(NULL, NULL);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  Get_input(my_rank, comm_sz, &a, &b, &n);

 // if(my_rank == 0) 
//    printf("my_rank=%d, a=%15.14lf, b=%15.14lf, number of total steps=%lf\n", my_rank, a, b, n);

  local_n = comm_sz;
  local_a = my_rank*(b/comm_sz);//used formula form sumdigits()
  local_b = my_rank*(b/comm_sz) + (b/comm_sz);

  integrated_area = area(local_a,local_b,local_n);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&integrated_area, &total_int_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(my_rank==0)
    cout << "The area is " << total_int_area <<endl;
  MPI_Finalize(); 
  return 0;

}
void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, double* n_p) {
  int rc=0;

  if (my_rank == 0) {
    printf("Enter a, b, and n\n");
    rc=scanf("%lf %lf %lf", a_p, b_p, n_p); 
    if(rc < 0) 
      perror("Get_input");
  }//if

  MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}  /* Get_input */

double area(double a, double b, double n){
  double dX = 0.01/ n;
  double area = 0.0;
  double x = a;
  while(x <= b){
    area += sin(x)*dX;
    x += dX;
  }

  return area;
}
