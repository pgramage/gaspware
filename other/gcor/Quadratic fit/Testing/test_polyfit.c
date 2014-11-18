#include <stdio.h>
#include "polynomial_fit.h"

//compile with $gcc -o testpoly.o test_polyfit.c -lm -lgsl -lgslcblas                                                                                                                          
                                                                                    

#define NP 11
double x[] = {0,  1,  2,  3,  4,  5,  6,   7,   8,   9,   10};
double y[] = {1,  6,  17, 34, 57, 86, 121, 162, 209, 262, 321};
double e[] = {0.1,  0.6,  0.17, 0.34, 0.57, 0.86, 0.121, 0.162, 0.209, 0.262, 0.321};
#define DEGREE 2
double coeff[DEGREE];
 
int main()
{
  int i;
  double chisq=0;
 
  polynomialfit(NP, DEGREE, x, y, e, chisq, coeff);
  for(i=0; i<DEGREE; i++)   printf("%g\t", coeff[i]);
  printf("%g\n", chisq);
  return 0;
}