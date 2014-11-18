/*// GainCor for LaBr3:Ce spectra - R. Lica - Dec. 2013

 Automatic Gain Correction relative to a reference spectra based on minimised chi sqared (not a peak finder)

 Compile with: $gcc -o gcor2 gcor2.c -lm -lgsl -lgslcblas
 You must install the graphical interface for gnuplot 
   (Ubuntu) $sudo apt-get install gnuplot-x11


------------------------------
==Version 1.3 (20 Dec 2013)
-a region will be auto-centered on the highest channel inside it
-if the fit is not good, the program will skip the output

==Version 1.2 (14 Dec 2013): 
-the user will set which regions to consider -> autosearch is not reliable yet
   Disadvantage: the same regions will be considered for all the detectors 
                 -> it is recommended that the reference run to have all 
                    the detectors aligned
-added plotting of raw and derivative spectra
-any run can be set as reference
 
==Version 1.1 ( 5 Dec 2013): 
-variable width, smoothing of spectra, derivative.

==Version 1.0 (29 Nov 2013): 
-can read one detector from several runs, normalization of spectra.
------------------------------
*/
  


#include "gcor2_definitions.h"

/* Global variables
   
int chNum, detNum, runref, runstop;
int sweep, degree;  
int chan[100], width[100];
char name[20];
int regions;
double sens;
int low, high;
    
*/


int main(int argc, char **argv) {

  
  int i, j, k, onetime=0;
  printf("\n======= GCOR v1.2 - Automatic Gain Correction =======\n\tR. Lica, IFIN-HH, Dec2013 \n\n");
  
  
  //Initializing 
    
  regions=initialize();
  degree++; // in the program, 'degree' represents the number of coefficients of the polynomial, NOT the degree of the polynomial
  if (degree<2 || degree>5) {printf("ERROR: Degree must be between 1 and 4"); exit(0);}
  low=chan[0]-width[0]/2;
  high=chan[regions-1]+width[regions-1]/2+sweep;
  
  
  struct Data2Fit shData[regions];
  double refSpec[chNum], Spec[chNum], smoothSpec[high-low];
  int nData=0;
  double coeff[degree], chisq, norm;
  char answer;
  char title[200];
  FILE * gnuplotPipe = popen ("gnuplot", "w");
  FILE * gnuplotPipe2 = popen ("gnuplot", "w");
  FILE * gnuplotPipe3 = popen ("gnuplot", "w");
  FILE *fo;
  char outfile[20];
  sprintf(outfile, "gcor.cal");
  fo=fopen(outfile, "wt");
  
  //Reading the data
  int irun[runstop], irunref, idet, runCount;
  int **rawSpec; // same as  'int rawSpec[runstop][detNum*chNum];' but can have any size
  rawSpec = calloc(runstop, sizeof(int *));
  for(i = 0; i < runstop; i++) rawSpec[i] = calloc(detNum*chNum, sizeof(int));
  runCount = readData(rawSpec, irun, &irunref); 
   
  
  
  //GCOR main code
  
  for(idet=0; idet<detNum; idet++)
  {
    onetime=0;
      
  
  for(i=0; i<runCount; i++) {
    if(irun[i]!=runref)
  {
        
    for (j=0; j<regions; j++) 
    {
      shData[j].ch=0;
      shData[j].chShift=0;
      shData[j].err=0;
    }
  
  
    for (j=0; j<chNum; j++)
    {
      refSpec[j]=rawSpec[irunref][idet*chNum+j];
      Spec[j]=rawSpec[i][idet*chNum+j];
    }
    
  if (answer!='a') {
  sprintf(title, "Raw Spectra");
  gnuplot_spec(gnuplotPipe3, refSpec, Spec, title); 
  }
  
  //if (onetime==0) { centerMax(refSpec); onetime++;}
  smooth(refSpec, Spec, 10, 30);
  for (j=low; j<high; j++) smoothSpec[j]=Spec[j];
   
  norm = normalize(refSpec, Spec);
  deriv(refSpec, Spec, 3);
  shift(smoothSpec, refSpec, Spec, shData); 
  nData = performFit(shData, &chisq, coeff);
  if (nData == 0) 
  {
    printf("Warning! %2s#%02d.%04d: Could not extract data suitable for fit\n", name, idet, irun[i]);
    fprintf(fo, "%5d%5d%5d%9.3f%10.6f\n", irun[i], idet, 2, 0.0, 1.0);
    goto skip;
    
  }
   
  if (answer!='a') //disable the display in automatic mode
  {
    gnuplot(gnuplotPipe, irun[i], idet, shData, nData, chisq, coeff, norm);
    sprintf(title, "Derivative");
    gnuplot_spec(gnuplotPipe2, refSpec, Spec, title);
  }
  
  writeCal(fo, coeff, irun[i], idet);
  
  if (answer=='a') goto skip;
  printf("%2s#%02d.%04d\t Going to next? [y]/n/a\t", name, idet, irun[i]);
  answer = getchar();
  skip:
  if (answer=='n') exit(0);
  
    
  }
  
  else fprintf(fo, "%5d%5d%5d%9.3f%10.6f\n", runref, idet, 2, 0.0, 1.0); //each detector from the reference run is set as reference
  }
  
  printf("\n---------------------\nGoing to Detector #%02d\n---------------------\n", idet+1);
  
  }
   
  exit(0); 
  
  
}