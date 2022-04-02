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
float ex4_accel(float time);
float ex4_vel(float time);
float ex4_pos(float time, float prevTime, float velocity, float prevPos);
float LeftRiemann(float left_endpt, float right_endpt, int rect_count, float base_len);
void Get_input(int my_rank, int comm_sz, float* a_p, float* b_p, float* n_p);
float funct_to_integrate(float x);
//limited precision of real values

int main() {
  int my_rank=0, comm_sz=0;
  float local_n=0;   
  float local_a=0.0, local_b=0.0;
  float local_int_area=0.0, total_int_area=0.0;
  float a=0.0, b=0, n=0.0;
  float step_size=0.0;
  vector<float> accelValues;
  vector<float> velValues;
  vector<float> posValues;

  /* Let the system do what it needs to start up MPI */
  MPI_Init(NULL, NULL);

  /* Get my process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  Get_input(my_rank, comm_sz, &a, &b, &n);

  if(my_rank == 0) 
    printf("-----my_rank=%d, a=%7.6f, b=%7.6f, number of total steps=%f\n", my_rank, a, b, n);

  step_size = (b-a)/n;  /* step is the same for all processes */
  local_n = n/comm_sz;  /* So is the number of rectangles  */

  local_a = a + my_rank*local_n*step_size;
  local_b = local_a + local_n*step_size;

  printf("*****my_rank=%d, start a=%f, end b=%f, number of rectangles=%f, step_size=%f\n",
      my_rank, local_a, local_b, step_size, local_n);
  local_int_area = LeftRiemann(local_a, local_b, step_size, local_n);

  printf("my_rank=%d, integrated area = %f, step_size * number rectangles=%f\n", 
      my_rank, local_int_area, (step_size*local_n));
  MPI_Reduce(&local_int_area, &total_int_area, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    printf("/////With n=%f rectangles, our estimate\n", n);
    printf("of the integral from %f to %f= %7.6f\n", a, b, total_int_area);
  }

  /* Shut down MPI */
  MPI_Finalize();


  return 0;
}

float LeftRiemann(float left_endpt, float right_endpt, int rect_count, float base_len) 
{
  float left_value, x, area=0.0; 
  int i;

  // estimate of function on left side to forward integrate
  left_value = funct_to_integrate(left_endpt);
  x = left_endpt;

  for (i = 1; i <= rect_count; i++) 
  {
    area += left_value * base_len;

    // new values to add to area
    x += base_len;
    left_value = funct_to_integrate(x);
  }

  return area;

} /*  LeftRiemann  */


void Get_input(int my_rank, int comm_sz, float* a_p, float* b_p, float* n_p) {
  int rc=0;

  if (my_rank == 0) {
    printf("Enter a, b, and n\n");
    rc=scanf("%f %f %f", a_p, b_p, n_p); if(rc < 0) perror("Get_input");
  } 
  MPI_Bcast(a_p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
}  /* Get_input */



float funct_to_integrate(float x) 
{
  //return sin(x);
  return(ex4_accel(x));
  //return(ex4_vel(x));
}

float ex4_accel(float time)
{
  // computation of time scale for 1800 seconds
  static float tscale=1800.0/(2.0*M_PI);
  //static float tscale=1800.0/(2.0*3.141592654);
  // determined such that acceleration will peak to result in translation of 122,000.0 meters
  static float ascale=0.2365893166123;

  return (sin(time/tscale)*ascale);
}

// determined based on known anti-derivative of ex4_accel function
float ex4_vel(float time)
{
  // computation of time scale for 1800 seconds
  static float tscale=1800.0/(2.0*M_PI);
  //static float tscale=1800.0/(2.0*3.141592654);
  // determined such that velocity will peak to result in translation of 122,000.0 meters
  static float vscale=0.2365893166123*1800.0/(2.0*M_PI);
  //static float vscale=0.2365893166123*1800.0/(2.0*3.141592654);

  return ((-cos(time/tscale)+1)*vscale);
}

float ex4_pos(float time, float prevTime, float velocity, float prevPos)
{
  float result=0.0;

  result=((velocity*(time-prevTime))+prevPos);

  return result;
}
