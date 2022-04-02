//Adapted from Dr. Sam Sewiert 'timeprofiles.c' file
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using std::cout;
using std::endl;

// direct generation of acceleration at any time with math library and arithmetic
double ex4_accel(double time);
// direct generation of acceleration at any time with right riemann sum
double Right_rsum(double a, double b, double n);

int main(int argc, char *argv[])
{
  double res = 0.0;

  for(int i=0; i < 1801; i++){
    res=ex4_accel((double)i);
    printf("%lf", (double)i);
    printf(",%08.7lf\n", res);
  }

  return 0;
}

double Right_rsum(double a, double b, double n){
  double dX = abs(b-a) / n;
  double area = 0.0;
  double x = a+dX;
  while(x <= b){
        area += sin(x)*dX;
    //    function to integrate below
    //area += pow(x, 2)*dX;
    x += dX;
  }

  return area;
}  /* integrate */

double ex4_accel(double time)
{
  // computation of time scale for 1800 seconds
  static double tscale=1800.0/(2.0*M_PI);
  // determined such that acceleration will peak to result in translation of 122,000.0 meters
  static double ascale=0.2365893166123;

  return (sin(time/tscale)*ascale);
}

