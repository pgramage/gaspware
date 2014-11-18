//Cristina Nita

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#define DATATYPE float

double sqr(double x);

int main( int argc, char *argv[] )
{

struct stat fileStat;
DATATYPE *data1, *data2, *data;
DATATYPE *err1, *err2, *err_tot;
FILE *inFile1;
FILE *inFile2;
FILE *inFile3;
FILE *file;

double K_S; 
int ii = 0, i = 0, j = 0, m = 0;

if( argc != 4 )
   {
   fprintf(stderr,"\n\nMissing Arguments! [data.r2] [background.r2] [final.r2] \n\n");
       exit(0);
   }
   
   printf("Input normalization constant: ");
   scanf("%lf", &K_S);

if( stat( argv[1], &fileStat ) )
   {
   fprintf(stderr,"Cannot open %s !\n", argv[1]);
       exit(0);
   }

   printf("%s: File size %d bytes ==> %d channels\n", argv[1], fileStat.st_size, fileStat.st_size/sizeof( DATATYPE ));
   
   data1 = calloc(fileStat.st_size/sizeof( DATATYPE ), sizeof( DATATYPE ));
   err1 = calloc(fileStat.st_size/sizeof( DATATYPE ), sizeof( DATATYPE ));
  
   inFile1 = fopen( argv[1], "rb" );
   
   fread( data1, sizeof( DATATYPE ), fileStat.st_size/sizeof( DATATYPE ), inFile1 );

   fclose( inFile1 );

for( ii = 0; ii < fileStat.st_size/sizeof( DATATYPE ); ii++ )
   {

        if(data1[ii] == 0){ err1[ii] = 1.0; }
           else { err1[ii] = sqrt(data1[ii]); }
   }


if( stat( argv[2], &fileStat ) )
   {
   fprintf(stderr,"Cannot open %s !\n",argv[2]);
       exit(0);
   }

   printf("%s: File size %d bytes ==> %d channels - Norm = %lf \n", argv[2], fileStat.st_size, fileStat.st_size/sizeof( DATATYPE ), K_S);
   
   data2 = calloc(fileStat.st_size/sizeof( DATATYPE ), sizeof( DATATYPE ));
   err2 = calloc(fileStat.st_size/sizeof( DATATYPE ), sizeof( DATATYPE ));

   inFile2 = fopen( argv[2], "rb" );
   
   fread( data2, sizeof( DATATYPE ), fileStat.st_size/sizeof( DATATYPE ), inFile2 );

   fclose( inFile2 );

ii=0;
for( ii = 0; ii < fileStat.st_size/sizeof( DATATYPE ); ii++ )
   {

        if( data2[ii] == 0 ){ err2[ii] = 1.0; }
           else { err2[ii] = sqrt( data2[ii] ); }
   }

if( stat( argv[3], &fileStat ) )
   {
   fprintf(stderr,"Cannot open %s !\n",argv[3]);
       exit(0);
   }

   printf("%s: File size %d bytes ==> %d channels\n", argv[3], fileStat.st_size, fileStat.st_size/sizeof( DATATYPE ));

   data = calloc(fileStat.st_size/sizeof( DATATYPE ), sizeof( DATATYPE ));
   err_tot = calloc(fileStat.st_size/sizeof( DATATYPE ), sizeof( DATATYPE ));

   inFile3 = fopen( argv[3], "rb" );
     
   fread( data, sizeof( DATATYPE ), fileStat.st_size/sizeof( DATATYPE ), inFile3 );

   fclose( inFile3 );

   
   file = fopen("output.dat","wt");

if(file != NULL){
ii=0;
for( ii = 0; ii < fileStat.st_size/sizeof( DATATYPE ); ii++ )
   {
     err_tot[ii] = sqrt( sqr( err1[ii] ) + ( sqr( K_S ) * sqr( err2[ii] ) ) );
     fprintf(file,"%d\t%f\t%f\n",ii, data[ii], err_tot[ii]);
     }
   fclose(file);
                  }
                  
   printf("----------> Printing in 'output.dat'\n\n");               
   return 0;
}

double sqr(double x)
{
  x = x * x;
  return( x );
} 