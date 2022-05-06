#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <math.h>

using std::cin;
using std::cout;
using std::endl;
using std::vector;


//function prototypes
double ex4_accel(double time);
double ex4_vel(double time);
double ex4_pos(double time, double prevTime, double velocity, double prevPos);
double Trap(double left_endpt, double right_endpt, double trap_count, double base_len);    
void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, double* n_p);
double funct_to_integrate(double x);
//limited precision of real values

int main() {
  int my_rank=0, comm_sz=0;
  double local_n=0;   
  double local_a=0.0, local_b=0.0;
  double local_int_area=0.0, total_int_area=0.0;
  double a=0.0, b=0, n=0.0;
  double step_size=0.0;
  vector<double> accelValues;
  vector<double> velValues;
  vector<double> posValues;

  /* Let the system do what it needs to start up MPI */
  MPI_Init(NULL, NULL);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  Get_input(my_rank, comm_sz, &a, &b, &n);

  if(my_rank == 0) 
    printf("-----my_rank=%d, a=%15.14lf, b=%15.14lf, number of total steps=%lf\n", my_rank, a, b, n);

  step_size = (b-a)/n;  /* step is the same for all processes */
  local_n = n/comm_sz;  /* So is the number of trapezoids  */

  local_a = a + my_rank*local_n*step_size;
  local_b = local_a + local_n*step_size;

  printf("*****my_rank=%d, start a=%lf, end b=%lf, number of trapezoids=%lf, step_size=%lf\n",
      my_rank, local_a, local_b, step_size, local_n);
  local_int_area = Trap(local_a, local_b, step_size, local_n);

  printf("my_rank=%d, integrated area = %lf, step_size * number trapezoids=%lf\n", 
      my_rank, local_int_area, (step_size*local_n));
  MPI_Reduce(&local_int_area, &total_int_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    printf("/////With n=%lf trapezoids, our estimate\n", n);
    printf("of the integral from %f to %f= %15.14lf\n", a, b, total_int_area);
  }

  /* Shut down MPI */
  MPI_Finalize();


  return 0;
}

double Trap(double left_endpt, double right_endpt, double trap_count, double base_len)
{
   double estimate, x; 
   int i;

   estimate = (funct_to_integrate(left_endpt) + funct_to_integrate(right_endpt))/2.0;

   for (i = 1; i <= trap_count-1; i++) 
   {
      x = left_endpt + i*base_len;
      estimate += funct_to_integrate(x);
   }
   estimate = estimate*base_len;

   return estimate;
}

void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, double* n_p) {
  int rc=0;

  if (my_rank == 0) {
    printf("Enter a, b, and n\n");
    rc=scanf("%lf %lf %lf", a_p, b_p, n_p); if(rc < 0) perror("Get_input");
  } 
  MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}  /* Get_input */

double funct_to_integrate(double x) 
{
  //return sin(x);
  return(ex4_accel(x));
  //return(ex4_vel(x));
}

double ex4_accel(double time)
{
  // computation of time scale for 1800 seconds
  static double tscale=1800.0/(2.0*M_PI);
  //static double tscale=1800.0/(2.0*3.141592654);
  // determined such that acceleration will peak to result in translation of 122,000.0 meters
  static double ascale=0.2365893166123;

  return (sin(time/tscale)*ascale);
}

// determined based on known anti-derivative of ex4_accel function
double ex4_vel(double time)
{
  // computation of time scale for 1800 seconds
  static double tscale=1800.0/(2.0*M_PI);
  //static double tscale=1800.0/(2.0*3.141592654);
  // determined such that velocity will peak to result in translation of 122,000.0 meters
  static double vscale=0.2365893166123*1800.0/(2.0*M_PI);
  //static double vscale=0.2365893166123*1800.0/(2.0*3.141592654);

  return ((-cos(time/tscale)+1)*vscale);
}

double ex4_pos(double time, double prevTime, double velocity, double prevPos)
{
  double result=0.0;
  result=((velocity*(time-prevTime))+prevPos);

  return result;
}
