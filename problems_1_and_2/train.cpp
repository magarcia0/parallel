#include <iostream>
#include <cmath>
#include <vector>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

//function prototypes
double ex4_accel(double time);
double ex4_vel(double time);
double ex4_pos(double time, double prevTime, double velocity, double prevPos);
double Right_rsum(double a, double b, double n);
void Get_input(double* a_p, double* b_p,double* n_p);

int main() {
  double a=0.0;
  double b=0.0; 
  double n=0.0;
  double integrated_area=0.0;
  vector<double> accelValues;
  vector<double> velValues;
  vector<double> posValues;

  //calculate acceleration values
  for(int i = 0; i < 1801; i++){
    double res=0.0;
    res=ex4_accel((double)i);
    accelValues.push_back(res);
    //printf(",%08.7lf\n", res);
    //cout<<accelValues[i]<<endl;
  }

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
  //take derivative of those values by calling
  //  ex4_vel() function to get velocities
  //Get_input(&a, &b, &n);

  //integrated_area = Right_rsum(a, b, n);

  // cout << "The area is " << integrated_area <<endl;

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

double Right_rsum(double a, double b, double n)
{
  double dX = abs(b-a) / n;
  double area = 0.0;
  double x = a+dX;
  while(x <= b){
    //function to integrate below
    area += sin(x)*dX;
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
