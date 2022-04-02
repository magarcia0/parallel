//Adapted from Dr. Sam Sewiert 'timeprofiles.c' file
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// direct generation of acceleration at any time with math library and arithmetic
double ex4_accel(double time);
double ex4_vel(double time);

// Implement methods of integration
double Local_Riemann(double a, double b, int n, double func(double));
//double Local_Trap(double a, double b, int n, double func(double));


int main(int argc, char *argv[])
{


   return 0;
}


double Local_Riemann(double a, double b, int n, double funct(double))
{
    double dt, interval_sum=0.0, local_a, local_b, time;
    int my_rank;
    int rank_count;

    dt = (b-a)/((double)n);

    int idx, local_n;

    local_n = n / rank_count;

    local_a = a + my_rank*local_n*dt;
    local_b = local_a + local_n*dt;


    for(idx=1; idx <= local_n; idx++)
    {
        time = local_a + idx*dt;
        interval_sum += (funct(time) * dt);
        //printf("Step for my_rank=%d at time=%lf, f(t)=%lf, sum=%lf\n", my_rank, time, funct(time), interval_sum);
    }

    //printf("Local Riemann = %lf for my_rank=%d of rank_count %d with dt=%lf, on a=%lf to b=%lf for %d steps\n",
    //        interval_sum, my_rank, rank_count, dt, local_a, local_b, local_n);

    return interval_sum;
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

