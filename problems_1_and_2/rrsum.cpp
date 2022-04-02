//This program adapts from L Muller's program
// 'Reimann sum in C++' from the video linked below
// https://www.youtube.com/watch?v=NWCnABJHqIQ&ab_channel=LMuller
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

double Right_rsum(double a, double b, double n);

void Get_input(double* a_p, double* b_p,double* n_p);

int main()
{
  double a=0.0;
  double b=0.0; 
  double n=0.0;
  double integrated_area=0.0;

  Get_input(&a, &b, &n);

  integrated_area = Right_rsum(a, b, n);

  cout << "The area is " << integrated_area <<endl;
  return 0;
}

void Get_input(double* a_p, double* b_p, double* n_p) {
  int rc=0;

  printf("Enter a, b, and n\n");
  rc=scanf("%lf %lf %lf", a_p, b_p, n_p); 
  if(rc < 0) 
    perror("Get_input");

}  /* Get_input */

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

