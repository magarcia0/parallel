//This program adapts from Dr. Sam Siewert's regula.c program
//
//Author: Marco Garcia, 4/29/22
//
//https://www.codewithc.com/c-program-for-regula-falsi-method/
//
// Refactored to improve variable names and comments
// 
// Sam Siewert, 4/26/22
//
// Regula Falsi does not require the derivative of the function to estimate a new interval
// that contains the root like Newton-Raphson, but rather uses the SECANT, which can simply be
// determined by evaluating the function at each end of the search interval.
//
// The intersection of the SECANT based upon slope, will ALWAYS contain the ZERO either on the 
// LEFT side of the search interval or the RIGHT side.
//
// We simply pick the side that contains the ZERO and iterate.
//
// Theory can be understood by starting with Wikipedia - https://en.wikipedia.org/wiki/Regula_falsi
//
// The double false position is the fastest converging and most reliable.
//
#include<stdio.h>
#include<math.h>
#include<time.h>

double f(double x)
{
  // As can be seen with Desmos, this example equation has roots near: -2.851, 4.758, and 7.792
  return (-(x*x*x) + 9.7*(x*x) -1.3*x -105.7);
}

void regula (double *c, double a, double b, double fb, double fa, int *itr)
{
  *c = b - ( fb*((b - a) / (fb - fa)) );

  ++(*itr);

  printf("Iteration no. %3d X = %7.5f \n", *itr, *c);
}


int main(void)
{
  int itr = 0, maxmitr;
  double x3, allerr;
  double a, b, c;
  struct timespec start, end;
  double fstart=0.0, fend=0.0;


  printf("\nEnter the values of x0, x1, allowed error and maximum iterations:\n");
  scanf("%lf %lf %lf %d", &a, &b, &allerr, &maxmitr);

  clock_gettime(CLOCK_MONOTONIC, &start);
  // Get the first value for the intersection of the SECANT on the interval
  regula (&c, a, b, f(b), f(a), &itr);

  for(int i=0; i < 3; i++) {
    do
    {
      // Recall that a ZERO crossing is anywhere where the value of the function at 2
      // different x values is negative.
      //
      // If we have a ZERO crossing between new x2 and x0, ZERO crossing between x1 and x2
      if (f(a)*f(c) < 0)
        b=c;

      // ELSE we have a ZERO crosssing between x0 and x2
      else
        a=c;

      if(i==2){
        //a=7;
        a=8;
        //b=8;
        b=7;
      }

      // Get the new value for the intersection of the SECANT on the new interval
      regula (&x3, a, b, f(b), f(a), &itr);


      if (fabs(x3-c) < allerr)
      {
        printf("After %d iterations, root = %20.15lf\n", itr, x3);
        a=itr;
        b=0;
        itr=0;
        if(i==2){
          clock_gettime(CLOCK_MONOTONIC, &end);
          fstart=start.tv_sec + (start.tv_nsec / 1000000000.0);
          fend=end.tv_sec + (end.tv_nsec / 1000000000.0);
          printf("\nSerial Regula Falsi took %lf seconds\n", (fend-fstart));

          return 0;
        }
        break;
      }
      c=x3;
    }
    while (itr<maxmitr);
  }

  printf("Solution does not converge or iterations not sufficient:\n");
  return 1;
}
