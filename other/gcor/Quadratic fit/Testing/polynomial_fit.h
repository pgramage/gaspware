#include <math.h>
#include <gsl/gsl_multifit.h>

//Modified by R. Lica. Taken from
//        http://rosettacode.org/wiki/Polynomial_regression
//        http://www.gnu.org/software/gsl/manual/html_node/Fitting-Examples.html

void polynomialfit(int n, int degree, 
		   double *xi, double *yi, double *ei, double *chisq, double *coeff)
{
  int i, j;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  X = gsl_matrix_alloc (n, degree);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);

  c = gsl_vector_alloc (degree);
  cov = gsl_matrix_alloc (degree, degree);

  for (i = 0; i < n; i++)
    {
      
      
      gsl_matrix_set (X, i, 0, 1.0);
      for(j=0; j < degree; j++)
      {
	gsl_matrix_set(X, i, j, pow(xi[i], j));
	
      }
      
      
      gsl_vector_set (y, i, yi[i]);
      gsl_vector_set (w, i, 1.0/(ei[i]*ei[i]));
    }

  {
    gsl_multifit_linear_workspace * work 
      = gsl_multifit_linear_alloc (n, degree);
    gsl_multifit_wlinear (X, w, y, c, cov,
                          chisq, work);
    gsl_multifit_linear_free (work);
  }
/*
#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {
    printf ("# best fit: Y = %g + %g X + %g X^2\n", 
            C(0), C(1), C(2));

    printf ("# covariance matrix:\n");
    printf ("[ %+.5e, %+.5e, %+.5e  \n",
               COV(0,0), COV(0,1), COV(0,2));
    printf ("  %+.5e, %+.5e, %+.5e  \n", 
               COV(1,0), COV(1,1), COV(1,2));
    printf ("  %+.5e, %+.5e, %+.5e ]\n", 
               COV(2,0), COV(2,1), COV(2,2));
    printf ("# chisq = %g\n", chisq);
  }
*/

  for (i=0; i<degree; i++) coeff[i]=gsl_vector_get(c,i);
  
  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (w);
  gsl_vector_free (c);
  gsl_matrix_free (cov);

}