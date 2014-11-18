//  Convlution Algorithm v2.0 (best fit) for decay spectra - R. Lica, April 2013

// R^2 = 1 - SSerr/SStot; SSerr = sumsq(fy-oy); SStot = sumsq(fy-fymean)

#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#define E 2.7182818284




int main(int argc, char **argv)
{

  int i, j, k;
  double step, hllim1, hllim2;
  double hl, hlfit, r2=0, maxr2=0.5, SSerr=0, SStot=0, fymean=0;
  int chan=0;
  double lambda, x[16384],fx[16384], iy[16384], fy[16384], oy[16384], boy[16384];
  char outfile[64];
  FILE *fi;
  FILE *ff;
  FILE *fo;
  
  
  if (argc < 3) {
    printf("\nConvolution v2.0 - ERROR:  [input file] [fit file] required as arguments:\n\n");
    exit(0);
  }
  fi=fopen(argv[1], "r");
  ff=fopen(argv[2], "r");
   
  
  printf("\n\n-- Convolution v2.0 - will calculate R^2 best fit of convoluted input data-\n");
  
  printf("Step:\t");
  scanf("%lf", &step);
  printf("Half-life limits (up, down):\t");
  scanf("%lf, %lf", &hllim1, &hllim2);
  
  hl=hllim1;
  
  
    
  while (fscanf (fi, "%lf", &x[chan])==1 && fscanf (fi, "%lf", &iy[chan])==1 && 
    fscanf (ff, "%lf", &fx[chan])==1 && fscanf (ff, "%lf", &fy[chan])==1) chan++;

      
  
  // R^2 = 1 - SSerr/SStot; SSerr = sumsq(fy-oy); SStot = sumsq(fy-fymean)
    
  
  for (i=0; i<chan; i++) fymean+=fy[i];
  fymean=fymean/chan;  
  for (i=0; i<chan; i++) SStot+=pow((fy[i]-fymean), 2);
    
    
  while(  hl<hllim2 )
  {
    
  hl+=step;  
  for (i=0; i<16384; i++) oy[i]=0;
  SSerr=0, fymean=0;
  
  lambda = log(2)/hl;
  for (i=1; i<chan; i++)
    for (j=0; j<=i; j++)
      oy[i]+=iy[j]*(lambda*pow(E,lambda*(x[j]-x[i])));    ////////  THE CONVOLUTION with exponential decay  
    
    
  
  for (i=0; i<chan; i++) SSerr+=pow((fy[i]-oy[i]), 2);
  r2=1-SSerr/SStot;
  if (r2>maxr2) 
  { 
    maxr2=r2; 
    hlfit=hl;
    for (i=0; i<chan; i++) boy[i]=oy[i];
  }
  
    
  }
  
  
  
  
        
  
      
  sprintf(outfile, "conv_%.2lf_%s", hlfit, argv[1]);
  fo=fopen(outfile, "wt");    
  for (i=0; i<chan; i++) fprintf(fo, "%.3lf\t%.3lf\n", x[i], boy[i]);
  printf("\nR^2:\t%.6lf\nHL:\t%.2lf\nConvoluted spectra written in ====== [%s] ======\n",maxr2, hlfit, outfile);
  
  exit(0);
}