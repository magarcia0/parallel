// This Program adapts from Dr. Sam Siewert's porgram 'newton.c'
//
// https://www.codewithc.com/c-program-for-newton-raphson-method/
//
// Refactored to improve variable names and commenting
//
// Sam Siewert, 4/26/2022
//
// For theory overview, start with Wikipedia - https://en.wikipedia.org/wiki/Newton's_method
//
// Draw some examples and convince yourself that the slope of the tangent (derviative of the function) at x, 
// helps with the search for the ZERO crossing, which is a root of the function.
//
#include<stdio.h>
#include<math.h>

double f(double x) {
  return ( (-1*pow(x, 3)) + (9.7*pow(x, 2)) - (1.3*x) - 105.7 );
}

double df (double x) {
  // As can be seen with Desmos, this example equation has roots near: -2.851, 4.758, and 7.792
  return ( (-3.0*pow(x, 2)) + (19.4*x) - 1.3 );
}

int main(void) {
  int itr=0, maxmitr=0;
  double h=0.0, x0=0.0, x1=0.0, allerr=0.0;

  printf("\nEnter x0, allowed error and maximum iterations\n");
  scanf("%lf %lf %d", &x0, &allerr, &maxmitr);

  for (itr=1; itr<=maxmitr; itr++)
  {
    h=f(x0)/df(x0);

    x1=x0-h;

    printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);

    if (fabs(h) < allerr)
    {
      printf("After %3d iterations, root = %20.15f\n", itr, x1);
      return 0;
    }
    x0=x1;
  }

  printf(" The required solution does not converge or iterations are insufficient\n");

  return 1;
}
