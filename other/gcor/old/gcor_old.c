// GainCor for LaBr3:Ce spectra - R. Lica - Dec. 2013
// Automatic Gain Correction relative to a reference spectra based on minimised chi sqared (not a peak finder)


/*
Compile with     $gcc -o gcor gcor.c -lm -lgsl -lgslcblas
You must install gnuplot-x11


///////////////////////////////////////////////////////////////////

Changelog:

Version 1.1 (5 Dec 2013): variable width, smoothing of spectra.
Version 1.0 (29 Nov 2013): can read one detector from several runs, normalization of spectra.

///////////////////////////////////////////////////////////////////
*/
  


#include "gcor_definitions.h"


int main(int argc, char **argv) {

  
  
  printf("\n======= GCOR v1.1 - Automatic Gain Correction =======\n\tR. Lica, IFIN-HH, Dec2013 \n\n");
  
  
  
  int chNum, detNum, runstart, runstop, irun, idet;
  int minWIDTH, maxWIDTH, SHIFT, SWEEP, degree;  
  int low, high;
  char name[20];
  
  FILE *settings, *fi, *fo;
  initialize(settings, &chNum, &detNum, name, &runstart, &runstop, &minWIDTH, &maxWIDTH, &SHIFT, &SWEEP, &degree, &low, &high);
  
  degree++; // in the program degree represents the number of coefficients of the polynomial and NOT the degree of the polynomial. i know...
  int regions=(high-low)/SHIFT;
  struct Data2Fit shData[regions];
  int i, j, k, nData=0;
  double refSpec[chNum], Spec[chNum];
  double coeff[degree], chisq, norm;
  FILE * gnuplotPipe = popen ("gnuplot", "w");
  char outfile[20], infile[20], answer;
  sprintf(outfile, "gcor.cal");
  fo=fopen(outfile, "wt");
  
  

       
  for(idet=0; idet<detNum; idet++)
  {
    
  fprintf(fo, "%5d%5d%5d%9.3f%10.6f\n", runstart, idet, 2, 0.0, 1.0); //each detector from the first run is set as reference  
  
  for(irun=runstart+1; irun<=runstop; irun++)
  {
    
    for (i=0; i<regions; i++) 
    {
      shData[i].ch=0;
      shData[i].chShift=0;
      shData[i].err=0;
    }
  
  sprintf(infile, "%2s.%04d", name, runstart);  
  if(fopen(infile, "rb")) fi=fopen(infile, "rb");
  else {printf("Cannot open %2s.%04d\n", name, runstart); exit(0); }
  readBin(fi, refSpec, idet, chNum, low, high);
  
  sprintf(infile, "%2s.%04d", name, irun);
  if(fopen(infile, "rb")) fi=fopen(infile, "rb");
  else continue;
  readBin(fi, Spec, idet, chNum, low, high);
  
    
  normalize(refSpec, Spec, chNum, &norm);
  smooth(refSpec, Spec, minWIDTH/5, maxWIDTH/2, low, high);
  deriv(refSpec, Spec, low, high, 5);
  
  
  autoshift(refSpec, Spec, shData, low, high, minWIDTH, maxWIDTH, SHIFT, SWEEP, regions); 
  
  
  nData = performFit(shData, degree, &chisq, coeff, regions);
  if (nData == 0) {printf("Warning! %2s#%02d.%04d: Could not extract data suitable for fit. Change settings!\n", name, idet, irun); goto skip;}
  gnuplot(gnuplotPipe, irun, idet, shData, nData, degree, chisq, coeff, regions, norm);
  writeCal(fo, degree, coeff, irun, idet);
  
  if (answer=='a') goto skip;
  printf("Going to %2s#%02d.%04d ([y]/n)? (Type 'a' for automatic fit)", name, idet, irun+1);
  answer = getchar();
  skip:
  if (answer=='n') exit(0);
  
    
  }
  
  printf("---------------------\nGoing to Detector #%02d\n---------------------\n", idet+1);
  
    
  
  
  }
   
  exit(0); 
  
  
}