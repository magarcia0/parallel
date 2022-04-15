#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double getcoeff(void);

void main(int argc, char *argv[])
{
    time_t elapsedsec;
    double coeff;
    int dim=10, idx, jdx;

    if(argc > 1)
    {
        sscanf(argv[1], "%d", &dim);
    }

    // initialize random coeff generator
    elapsedsec=time((time_t *)0);
    srand((unsigned int)elapsedsec);

    printf("#Auto generated linear test\n");
    printf("%d\n", dim); 
     
    for(idx=0; idx < dim; idx++)
    {
        for(jdx=0; jdx < dim; jdx++)
        {
           printf("%lf ", getcoeff()); 
        }
        printf("\n");
    }

    for(idx=0; idx < dim; idx++)
    {
       printf("%lf\n", getcoeff()); 
    }

}

double getcoeff(void)
{
    int randnum1, randnum2;
    double coeff;

    randnum1=rand();
    randnum2=rand();
    coeff = ((double)randnum1 / (double)RAND_MAX) + (100.0 * ((double)randnum2 / (double)RAND_MAX));

    return(coeff);
}

