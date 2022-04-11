//Last I did and figured out was how to use the right riemann sum and ex4_vel.h to get velocity and the other method to get acceleration
//
#include <iostream>
#include <cmath>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

#define TENTH 0.1
#define HUNDREDTH 0.01

//function prototypes
double ex4_accel(double time);
double ex4_vel(double time);
double ex4_pos(double time, double prevTime, double velocity, double prevPos);
double Right_rsum(double a, double b, double step_size);
double LeftRiemann(double left_endpt, double right_endpt, int rect_count, double base_len);
void Get_input(double* a_p, double* b_p,double* n_p);
double funct_to_integrate(double x);
//limited precision of real values

int main() {
  double a=0.0;
  double b=1800; 
  double n=0.1;
  double step_size=0.0;
  double integrated_area=0.0;
  vector<double> accelValues;
  vector<double> velValues;
  vector<double> posValues;

  //calculate acceleration values
 /* for(int i = 0; i < 1801; i++){
    double res=0.0;
    res=ex4_accel((double)i);
    accelValues.push_back(res);
    //printf(",%08.7lf\n", res);
    //cout<<accelValues[i]<<endl;
  }*/
   double res=0.0;
   step_size = (b-a)/n;  /* step is the same for all processes */
   //res=LeftRiemann(a, b, step_size, n);
   res=LeftRiemann(a, b, step_size, n);
   printf("%15.14lf\n", res);
/*
 *
  for(int i = 0; i < 1801; i++){
    double res=0.0;
    res=ex4_vel((double)i);
    velValues.push_back(res);
    cout<<velValues[i]<<endl;
  }

  for(int i=0; i < 1801; i++){
    if(i == 0){
      double res=0.0;
      res=ex4_pos((double)i, 0, velValues[i], 0);
      posValues.push_back(res);
    }else if(i != 0){
      double res=0.0;
      res=ex4_pos((double)i, (double)(i-1), velValues[i], posValues[i-1]);
      posValues.push_back(res);
    }
    cout<<posValues[i]<<endl;
  }
  */

  return 0;
}

void Get_input(double* a_p, double* b_p, double* n_p)
{
  int rc=0;

  printf("Enter a, b, and n\n");
  rc=scanf("%lf %lf %lf", a_p, b_p, n_p); 
  if(rc < 0) 
    perror("Get_input");

}  /* Get_input */

double Right_rsum(double a, double b, double step_size)
{
  
  double dX = b-a / step_size;
  //double dX=n;
  double area = 0.0;
  double x = a + dX;

  while(x <= b){
    //function to integrate below
    area += funct_to_integrate(x)*dX;
    x += dX;
    printf("%lf\n", area);
  }

  return area;
}  /* integrate */

/*
double LeftRiemann(double left_endpt, double right_endpt, int rect_count, double base_len) 
{
   double left_value, x, area=0.0; 
   int i;

   // estimate of function on left side to forward integrate
   left_value = funct_to_integrate( (left_endpt+base_len) );
   x = left_endpt+base_len;

   cout<< "rect count: "<< rect_count <<endl;
   for (i = 1; i < rect_count+1; i++) 
   {
      area += left_value * base_len;

      // new values to add to area
      x += (left_endpt+base_len) + base_len;
      left_value = funct_to_integrate(x);
   }

   return area;

}
*/
double LeftRiemann(double left_endpt, double right_endpt, int rect_count, double base_len) 
{
   double left_value, x, area=0.0; 
   int i;
   cout<< "rect count: "<< rect_count <<endl;

   // estimate of function on left side to forward integrate
   left_value = funct_to_integrate(left_endpt);
   x = left_endpt;

   for (i = 1; i <= rect_count-1; i++) 
   {
      area += left_value * base_len;

      // new values to add to area
      x += left_endpt + base_len;
      left_value = funct_to_integrate(x);
   }

   return area;

} /*  LeftRiemann  */




double funct_to_integrate(double x) 
{
    //return sin(x);
 //   return(ex4_accel(x));
    return(ex4_vel(x));
}

double ex4_accel(double time)
{
  // computation of time scale for 1800 seconds
  static double tscale=1800.0/(2.0*M_PI);
  // determined such that acceleration will peak to result in translation of 122,000.0 meters
  static double ascale=0.2365893166123;

  return (sin(time/tscale)*ascale);
}

// determined based on known anti-derivative of ex4_accel function
double ex4_vel(double time)
{
  // computation of time scale for 1800 seconds
  static double tscale=1800.0/(2.0*M_PI);
  // determined such that velocity will peak to result in translation of 122,000.0 meters
  static double vscale=0.2365893166123*1800.0/(2.0*M_PI);

  return ((-cos(time/tscale)+1)*vscale);
}

double ex4_pos(double time, double prevTime, double velocity, double prevPos)
{
  double result=0.0;

  result=((velocity*(time-prevTime))+prevPos);

  return result;
}
