//  Convlution Algorithm for decay spectra - R. Lica, April 2013



#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#define E 2.7182818284




int main(int argc, char **argv)
{

  int i, j, k;
  char outfile[64];
  FILE *fi;
  FILE *fo;
  
  
  if (argc < 2) {
    printf("\nConvolution v1.0 - ERROR:  input file required as argument:\n\n");
    exit(0);
  }
  fi=fopen(argv[1], "r");
  
   
  
  printf("\n\n-- Convolution v1.0 --\n");
  
  
  
  double hl;
  int chan=0;
  double lambda, x[16384], iy[16384], oy[16384];
  
        
  while (fscanf (fi, "%lf", &x[chan])==1 && fscanf (fi, "%lf", &iy[chan])==1) chan++;

      
  start:
  printf("\nDecay half-life (same units as input):");
  scanf("%lf", &hl);
  
  for (i=0; i<16384; i++) oy[i]=0;
  sprintf(outfile, "conv_%.2lf_%s", hl, argv[1]);
  fo=fopen(outfile, "wt");
  lambda = log(2)/hl;
  for (i=1; i<chan; i++)
    for (j=0; j<=i; j++)
      oy[i]+=iy[j]*(lambda*pow(E,lambda*(x[j]-x[i])));    ////////  THE CONVOLUTION with exponential decay
  
      
  for (i=0; i<chan; i++) fprintf(fo, "%.3lf\t%.3lf\n", x[i], oy[i]);
    
    
    
  printf("\nConvolution v1.0 - convoluted spectra written in ====== ['%s'] ======\n", outfile);
  printf("\nTry again? (y/n):");
  char ask[10];
  scanf("%s", ask);
  if (ask[0]=='n') exit(0); 
  else goto start;
  
}