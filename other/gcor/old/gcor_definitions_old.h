#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include <gsl/gsl_multifit.h>

typedef struct Data2Fit { 
  int ch;
  int chShift;
  double err;
} Data2Fit;


void initialize(FILE *settings, int *chNum,   int *detNum, char *name, \
                  int *runstart,  int *runstop, int *minWIDTH, int *maxWIDTH,  int *SHIFT, \
                  int *SWEEP,     int *degree,  int *low,    int *high) {
  
  if (fopen("gcor_settings.txt", "rt")) 
  {
    settings = fopen("gcor_settings.txt", "rt");
    printf("Settings taken from 'gcor_settings.txt'.\n");
    
  }
  else
  {
    printf("Printing default settings in 'gcor_settings.txt'.\n");
    if(fopen("gcor_settings.txt", "wt"))
    {
      settings = fopen("gcor_settings.txt", "wt");
      
    }
    else { printf("ERROR: Cannot create 'gcor_settings.txt'\n"); exit(0); }
    fprintf(settings, \
"Channels=8192\n\
DetNum=1\n\
RefFile=L0.0001\n\
EndFile=L0.0100\n\
minWIDTH=140\n\
maxWIDTH=240\n\
SHIFT=150\n\
SWEEP=120\n\
Degree=1\n\
Low-Limit=50\n\
High-Limit=3000\n\
");
    
    fclose(settings);
  }
  
  
  settings = fopen("gcor_settings.txt", "rt");
  fscanf(settings, \
"Channels=%d\n\
DetNum=%d\n\
RefFile=%2s.%04d\n\
EndFile=%2s.%04d\n\
minWIDTH=%d\n\
maxWIDTH=%d\n\
SHIFT=%d\n\
SWEEP=%d\n\
Degree=%d\n\
Low-Limit=%d\n\
High-Limit=%d\n", chNum, detNum, name, runstart, name, runstop, minWIDTH, maxWIDTH, SHIFT, SWEEP, degree, low, high);
  
  
  
  printf("------------\nChan\t%d\n\
DetNum\t%d\n\
RefFile\t%2s.%04d\n\
EndFile\t%2s.%04d\n\
minWIDTH\t%d\n\
maxWIDTH\t%d\n\
SHIFT\t%d\n\
SWEEP\t%d\n\
Degree\t%d\n\
Low\t%d\n\
High\t%d\n-------------\n", *chNum, *detNum, name, *runstart, name, *runstop, *minWIDTH, *maxWIDTH, *SHIFT, *SWEEP, *degree, *low, *high);

    
}

void polynomialfit(int n, int degree, double *xi, double *yi, double *ei, double *chisq, double *coeff) {

  //Modified by R. Lica. Taken from
  //        http://rosettacode.org/wiki/Polynomial_regression
  //        http://www.gnu.org/software/gsl/manual/html_node/Fitting-Examples.html

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

void readAscii(FILE *fi, int *spec, int chNum)  {    
    int ii;
    for (ii=0; ii<chNum; ii++) if( fscanf (fi, "%*d %d", &spec[ii]) == 0) break;
    fclose(fi);
  }
  
void readBin(FILE *fi, double *data, int idet, int chNum, int lowlimit, int highlimit) {
  
  int temp[chNum];
  int i=0;
  for(i=0; i<idet+1; i++)  fread(temp, sizeof(int),chNum,fi);
  
  //Applying low and high limits;
  
  for (i=0; i<chNum; i++) {
    data[i]=temp[i];
    if(i<lowlimit || i>highlimit) data[i]=0;
  }
  
  fclose(fi);
}

void deriv(double *spec1, double *spec2, int low, int high, int width) {
  
  //FILE *out;
  //out=fopen("Spec-der.dat", "wt");
  int i;
  double der1[high-width], der2[high-width];
  for (i=low; i<high-width; i++) {
    der1[i]=(spec1[i+width]-spec1[i])/width;
    der2[i]=(spec2[i+width]-spec2[i])/width;
  }
  
  for (i=low; i<high-width; i++) {
    //fprintf(out, "%d\t%lf\t%lf\t%lf\t%lf\n", i, spec1[i], der1[i], spec2[i], der2[i]);
    spec1[i]=der1[i];
    spec2[i]=der2[i];
    
    
  }
  
  
}

void smooth(double *spec1, double *spec2, int minW, int maxW, int low, int high) {
  
  int i, j;
  int ilow=low+minW/2, ihigh=high-maxW/2;
  double sum1=0, sum2=0, width=0;
  
  for (i=ilow; i<ihigh; i++)
  { 
    sum1=0, sum2=0, width=0;
    for (j=-(minW+(maxW-minW)*i/(ihigh-ilow))/2; j<(minW+(maxW-minW)*i/(ihigh-ilow))/2; j++)
    {
      sum1+=spec1[i+j];
      sum2+=spec2[i+j];
      width++;
    }
    spec1[i]=sum1/width;
    spec2[i]=sum2/width;
  
    //fprintf(fo, "%d\t%d\n", i, spec1[i]); 
  }
  
  
  
}

void normalize(double *spec1, double *spec2, int n, double *norm) {
  
  
  int i;
  double sum1=0, sum2=0;
  *norm=0;
  for (i=0; i<n; i++) { sum1+=spec1[i]; sum2+=spec2[i]; }
  if (sum1>sum2)
    for (i=0; i<n; i++)
    {
      *norm = sum1/(sum2+1);
      spec2[i]=spec2[i]**norm;
    }
  else
    for (i=0; i<n; i++)
    {
      *norm = sum2/(sum1+1);
      spec1[i]=spec1[i]**norm;
    }
  //printf("Norm. = %.3lf\t", norm);
  
}

double chisq(double *spec1, double *spec2, int l1, int l2, int n) { //chisq of two regions
  
  int i;        
  double sum=0;
  double c2=0;
  
  for (i=0; i<n; i++)
  {
    c2 += pow(( spec1[l1+i]- spec2[l2+i]), 2);
    sum+= abs(spec1[l1+i])+abs(spec2[l2+i]);
  }
  
  return c2/sum;  
}

int area(double *spec, int l, int n) {
  int i;
  double a=0;
  for (i=l; i<l+n; i++) a+=spec[i];
  return a;
}
    
double ratio(double *spec, int l, int n) { //checks if there is a significant change in the region; ration=max/min
  int ii;
  double min=spec[l], max=spec[l];
  double ratio=0;
  for (ii=l; ii<l+n; ii++)
  {
    if (spec[ii] < min ) min = spec[ii];
    if (spec[ii] > max ) max = spec[ii];
  }
  
  
  ratio = max - min; //not to divide by zero
  if (max - min > 2*sqrt(max) ) return ratio;
  else return 0;
    
}

void autoshift(double *spec1, double *spec2, struct Data2Fit *shData, \
           int low, int high, int minWIDTH, int maxWIDTH, int SHIFT, int sweep, int regions) { // finds shift and error for each region by sweeping and minimizing chiSq
 
  int i, j;
  double minchi;
  int WIDTH[regions];    //for implementing variable width
  for (i=0; i<regions; i++)
    WIDTH[i]=minWIDTH+(maxWIDTH-minWIDTH)*i/regions;
  
  for(i=0; i<regions; i++) 
  {
    shData[i].ch=0;
    shData[i].chShift=0;
    shData[i].err=0;
  }
  
  for (i=(low + sweep)/SHIFT; i<(high-sweep)/SHIFT; i++) 
  {
    minchi=chisq(spec1, spec2, i*SHIFT, i*SHIFT, WIDTH[i]);
    if (ratio(spec2, i*SHIFT, WIDTH[i]))
    {
      for (j=-sweep; j<sweep; j++)
      {
	if( minchi >= chisq(spec1, spec2, i*SHIFT+j, i*SHIFT, WIDTH[i]))
	{
	  
	  minchi=chisq(spec1, spec2, i*SHIFT+j, i*SHIFT, WIDTH[i]);
	  
	  shData[i].ch = i*SHIFT + WIDTH[i]/2; 
	  shData[i].chShift = j;
	  //shData[i].err = 2;
	  shData[i].err = sqrt(minchi)*5;
	  //shData[i].err = sqrt(minchi)*area(spec2, low, high)/(area(spec2, i*SHIFT, WIDTH[i]));
	  //printf("%d\t%d\t%f\t%d\n", shData[i].ch, shData[i].chShift, shData[i].err, WIDTH[i]); //testing
	  if (shData[i].err > 20) shData[i].ch = -1;
	}
	
      }
      
	
    }
    else shData[i].ch = -1;
    
  }
  
}

int performFit(struct Data2Fit *shData, int degree, double *chisq, double *coeff, int regions) {
  
  int i, j=0, n=0;
  for (i=0; i<regions; i++)
    if(shData[i].ch > 0) n++;
  double xi[n], yi[n], ei[n];
  
  for (i=0; i<regions; i++) 
    if(shData[i].ch > 0)
    { 
      xi[j]=shData[i].ch;
      yi[j]=shData[i].chShift; 
      ei[j]=shData[i].err;      
      j++;
      
      //printf("%d\t%d\t%f\n", shData[i].ch, shData[i].chShift, shData[i].err); //testing
    
      
    }
  if (j<2) n=0; //not enough points to perform fit
  else polynomialfit(n, degree, xi, yi, ei, chisq, coeff);
  
  /*
  for(i=0; i < degree; i++) {
    printf("%lf\n", coeff[i]);
  }
  printf("%g\n", *chisq/n); // normalised chi sqared = chisq/degrees of freedom
  */
  return n;
    
}

void gnuplot(FILE * gnuplotPipe, int irun, int idet, struct Data2Fit *shData, int n, int degree, double chisq, double *coeff, int regions, double norm) {
  
  
  
  char title[200], plot1[200]; 
  sprintf(title, "set title \"DET#%02d RUN#%04d  Chisq=%.2lf  Norm=%.2lf\"", idet, irun, chisq/n, norm);
  if(degree == 2) sprintf(plot1, "plot 'data.temp' using 1:2:3 with yerrorbars, %lf + %lf*x", coeff[0], coeff[1]);
  else if(degree == 3) sprintf(plot1, "plot 'data.temp' using 1:2:3 with yerrorbars, %lf + %lf*x + %lf*x*x", coeff[0], coeff[1], coeff[2]);
  else if(degree == 4) sprintf(plot1, "plot 'data.temp' using 1:2:3 with yerrorbars, %lf + %lf*x + %lf*x*x %lf*x*x*x", coeff[0], coeff[1], coeff[2], coeff[3]);
  else { printf("ERROR: Degree must be max 3\n"); exit(0);}
  
  
  FILE * temp = fopen("data.temp", "wt");
  //FILE * gnuplotPipe = popen ("gnuplot", "w");
  int i;
  for (i=0; i < regions; i++)
    {
      if(shData[i].ch > 0)
	fprintf(temp, "%d %d %lf \n",  shData[i].ch, shData[i].chShift, shData[i].err); //Write the data to a temporary file
    }
   

    fprintf(gnuplotPipe, "%s \n %s \n", title, plot1); //Send commands to gnuplot one by one.
    fflush(gnuplotPipe);
  fclose(temp);
    
    
}

void writeCal(FILE *fo, int Degree, double *coeff, int irun, int idet) {
 
  int deg;
  fprintf(fo, "%5d%5d%5d", irun, idet, Degree);
  for (deg=0; deg<Degree; deg++)
      {
	if (deg==0) fprintf(fo, "%9.3f", coeff[deg]);
	else if (deg==1) fprintf(fo, "%10.6f", 1+coeff[deg]);
	else fprintf(fo, "%15.6E", coeff[deg]);
      }
  fprintf(fo, "\n");
    
	
}
      