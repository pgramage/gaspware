#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include <termios.h>

#define  MAXCOLS 3
#define  MAXPTS  16384
#define  NUMPARS 4

/* Carl Wheldon July 2009 */

/* To compile:
    gcc halflife.c -Wall -pedantic -o halflife -lm -O2
*/

/*%%%% A program to fit a Gaussian convoluted with an exponential decay %%%%*/

struct ft {
  double pars[NUMPARS], errs[NUMPARS], y[16384], dy[16384], x[16384];
  int    freepars[NUMPARS], nfp, ndp;
} ft;


int 	ask_yn(char *questn);
int 	ascii_read(char fname[], int *col);
void 	ascii_write(char name[], int numch, int col);
void 	chi2(double *chisq);
void 	col_determ(FILE *file, int *col);
void 	compress2(char inname[], char outname[]);
double 	erfpol(double x);
double 	erfpol_new(double x);
int 	eval(double *pars, double x, double *fit, double *derivs, int mode);
void 	file_status(char *filename);
int 	fitter(int, double *chisqr, int vb);
void 	get_ans(char ans[], int num);
void 	get_line(char ans[]);
int 	get_mode(int md);
void 	get_pars(char ans0[], double pars[], int num);
void 	itoa(int n, char s[]);
int 	matinv(double *array, int nip, int npars);
void 	read_par(double pars[], int offset, int num, char *s);
void 	reverse(char s[]);
void 	scaler(double *fit, double *derivs, int mode);
void 	set_ext(char fname[], char ext[]);
void 	skip_hash(FILE *file);
int 	supp_zeroes(int numch);
void 	write_output(char inname[], char outname[], double init[],
	double x[], double *fit, double chisqr, int mode);

double fit, derivs[NUMPARS], bkgnd;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int main(int argc, char *argv[])
{
    extern double fit, derivs[NUMPARS], bkgnd;
    double  chisq = 0.0, init[NUMPARS+1], bklo = 0.0, chisqo = 10000.0;
    double  x12[2], tmpd, tmpx[16384];
    int     f = 0, i, md = 0, col = 0, flg = 0, r = -1, strt = 1, flg2 = 0;
    char    ans[80] = "", inname[80] = "", nm[80] = "", outname[80] = "";
    char    nofit[80] = "";
    struct  stat statbuf;
    
    /*initialise background parameter (constant offset) to zero*/
    bkgnd = 0.0;
    /*no. of data points*/
    ft.ndp = 0;
    /*no. of fixed parameters*/
    ft.nfp = 0;
    /*set all parameters to be free*/
    for (i = 0; i < NUMPARS; i++) ft.freepars[i] = 1;

    x12[0] = 0.0; x12[1] = 0.0;
    
    printf("\n\t**** Welcome to program halflife ****\nThis fits"
    	   " a Gaussian convoluted with an exponential decay.\n"
	    " The data can be in 1,2 or 3 column ascii format (x y dy):\n"
	    "  1 col: y; 2 col: x y; 3 col: x, y, dy\n\n");

    /*argv[i] is the ith argument, i.e. first is the program name*/
    if ( argc > 2 || (argc == 2 && (stat(argv[1], &statbuf))) )
    {
	printf("\nUnrecognised arguments...usage: halflife\n"
		" or: halflife InputFileName\n");
	if (argc == 2) printf(" ***File %s does not exist\n",argv[1]);
	return -1;
    }
    else if (argc == 2)
    {
	strcpy(inname,argv[1]);
	printf("Input filename: %s\n",inname);
	flg = 1;
    }
        
    /*zero arrays*/
    for (i = 0; i < NUMPARS; i++)
    {
	ft.pars[i] = 0.0;
	init[i] = 0.0;
    }
    
    while (1)
    {
	if ( strt == 1 || (md = get_mode(md)) == 1 )
	{
    	    if (!flg)
	    {
		printf("Type ascii data filename\n");
		if (strlen(inname) > 0 ) printf(" [<Enter> for %s]\n",inname);
    	    	get_line(nm);
/*    	    	if (((nm[0] == 'D') || (nm[0] == 'd')) && (nm[1] == ((char)0)));*/
    	    	if (strlen(nm) == (int)0 && (nm[0] == ((char)0))) ;
    	    	else strcpy(inname, nm);
	    }
	    flg = 0;
	    strt = 0;
    	    /*get output filename and set extension*/
    	    strncpy(outname, inname, 80);
    	    set_ext(outname, ".fit");
    	    printf("Output filename for fit: %s\n", outname);
    	    file_status(outname);
   	    strncpy(nofit, outname, 80);
    
    	    /*read in ascii data*/
    	    ft.ndp = ascii_read(inname, &col);
	}
	else if (md == 2)
	{
	    /*read in initial parameter values*/
	    read_par(init, 0, 1, "Enter initial guess for halflife (-ve if to the left of prompt peak)\n");	    
	    read_par(init, 1, 1, "Enter initial guess for FWHM of prompt peak");
    	    read_par(init, 2, 1, "Enter initial guess for prompt peak centroid");
    	    read_par(init, 4, 1, "Enter value for background constant background level (not fitted)\n");
	    	    
	    /*convert t_1/2 to -1*tau for use in eval()*/
	    ft.pars[0] = -1.0*init[0]/log(2);
	    /*convert FWHM to sigma for use in eval()*/
      	    ft.pars[1] = init[1]*0.6005612;
	    ft.pars[2] = init[2];
	    ft.pars[3] = 1.0;
	    bkgnd = 0.0;
	    /*now call scaler function to get pars[3]*/
    	    scaler(&fit, derivs, 0);
	    init[3] = ft.pars[3];
	    if (flg2)
	    {
    	    	read_par(init, 3, 1, "Initial guess for scaling factor");
	    	ft.pars[3] = init[3];
	    }
	    init[3] = ft.pars[3];
	    bkgnd = init[4];
	    printf("Initial parameter values:\n"
		   "t_1/2 = %.3f\n"
		   "FWHM = %.3f\n"
		   "Centroid = %.2f\n"
		   "Scaling factor = %.2f\n"
		   "Background level (not fitted) = %.2f\n",
		    init[0],init[1],init[2],init[3],init[4]);
	    flg2 = 1;	    
	}
	else if (md == 3)
	{
	    /*fix and free parameters*/
	    while (1)
	    {
	    	printf("Type the number of the parameter to be fixed/freed: * = fixed\n"
		    " %s   1) t_1/2\n"
		    " %s   2) FWHM\n"
		    " %s   3) Centroid\n"
		    " %s   4) Scaling Factor\n"
		    "        0) Finish\n",
    	    	    (ft.freepars[0] == 0) ? "  * " : "free",
	    	    (ft.freepars[1] == 0) ? "  * " : "free",
	    	    (ft.freepars[2] == 0) ? "  * " : "free",
    	    	    (ft.freepars[3] == 0) ? "  * " : "free");
		
    	    	get_ans(ans,1);
 	    	if (ans[0] >= '0' && ans[0] <= '4') f = ans[0] - '0';
		
		if (f == 0) break;
		else if (f > 0 && f < 5 && ft.freepars[f-1] == 0)
		{
		    ft.freepars[f-1] = 1;
		    ft.nfp--;
		}
		else if (f > 0 && f < 5 && ft.freepars[f-1] == 1)
		{
		    
    	    	    read_par(init, f-1, 1, "Enter value for fixed parameter");
		    ft.freepars[f-1] = 0;
		    printf("Parameter %d fixed as = %f\n",f,init[f-1]);
		    ft.pars[f-1] = init[f-1];
		    ft.nfp++;
		}
	    }	    
	}
	else if (md == 4 && flg2)
	{
	    /*multiply tau by -1.0 as a +ve tau spectrum goes right to left*/
	    ft.pars[0] = -1.0*init[0]/log(2);
	    
	    /*convert FWHM to sigma as defined in function eval()*/
      	    ft.pars[1] = init[1]*0.6005612;
	    
      	    ft.pars[2] = init[2];
      	    ft.pars[3] = init[3];
	    
	    /*perform fit and write output file*/
	    if ( (r = fitter(20, &chisq, 1)) == 0 )
	    {
		/*write output file*/
	    	write_output(inname, outname, init, ft.x, &fit, chisq, 0);
	    }
	    else if (r > 0 && r < 3)
		printf("Bad fit, try changing initial parameters\n");
	    else printf("Try freeing some parameters\n");
	}
	else if (md == 5)
	{
	    printf("Data will now be compressed by factor of 2\n"
	    	" Initial guesses for coefficient will also be compressed\n");
	    compress2(inname, outname);
	    /*change initial guesses for parameters*/
	    init[0] /= 2;
	    init[1] /= 2;
	    init[2] /= 2;
	    init[3] *= 2;
	    init[4] *= 2;
	    bkgnd *= 2;
	    ft.pars[0] = -1.0*init[0]/log(2);
      	    ft.pars[1] = init[1]*0.6005612;
	}
	else if (md == 6 && r == 0)
	{
	    /*check background parameter (chisq) for bkgnd +/- 1*/
	    bklo = bkgnd;
	    chisqo = 10000.0;
	    if (bkgnd - 1.0 > 0) bkgnd -= (double)1.0;
	    printf("  Bkgnd   t1/2  Chisq/D.O.F \n");
	    for (i = 0; i < 21; i++)
	    {
	    	if ( (r = fitter(20, &chisq, 0)) == 0 )
	    	{
		    printf("  %.2f  %.2f  %.3f\n",
			    bkgnd,(-1.0*log(2)*ft.pars[0]),chisq);
	    	    if (chisq < chisqo) bklo = bkgnd;
		    chisqo = chisq;
		}
		else printf(" Bad fit!\n");
/*		chi2(&chisq);
		printf("For bkgnd = %f --> chisq/D.O.F = %.3f\n",bkgnd,chisq);*/
		bkgnd += (double)0.1;
	    }
	    bkgnd = bklo;
/*    	    read_par(init, 4, 1, "Enter value for background constant background level (not fitted)\n");*/
	    init[4] = bkgnd;
	    /*write output file*/    
	    /*perform fit and write output file*/
	    if ( (r = fitter(20, &chisq, 1)) == 0 )
	    {
		/*write output file*/
	    	printf(" ==> Writing output for bkgnd = %f\n",bkgnd);
	    	write_output(inname, outname, init, ft.x, &fit, chisq, 0);
	    }
	    else printf(" Bad fit! Could not write output file.\n");
    	}
	else if ( md == 6 && r != 00)
	    printf(" ***Must get a good fit before exploring background level!***\n");	    
	else if ( md == 4 && ! flg2)
	    printf(" ***Enter initial parameter values before fitting!***\n");
	else if (md == 7) write_output(inname, outname, init, ft.x, &fit, chisq, 1);
	else if (md == 8)
	{
	    read_par(x12, 0, 2, "Enter range for x (x1 x2)");	    
	    if ((int)(x12[0] + 0.5) > (int)(x12[1]+0.5))
	    {
		tmpd = x12[1];
		x12[1] = x12[0];
		x12[0] = tmpd;
	    }
	    for (i = 0; i <= (int)(x12[1]+0.5-x12[0]);i++)
		tmpx[i] = x12[0] + i;
	    
	    r = ft.ndp;
	    ft.ndp = i;
	    /*write output file*/
    	    file_status(nofit);
	    printf(" ==> Writing output for bkgnd = %f\n",bkgnd);
	    write_output(inname, nofit, init, tmpx, &fit, chisq, 2);
	    ft.ndp = r;
	}
	else if (md == 9)
	{
	    while (1)
	    {
	    	printf("Enter number of columns required 2 or 3\n");
    	    	get_ans(ans,1);
 	    	if (ans[0] >= '2' && ans[0] <= '3')
		{
		    f = ans[0] - '0';
		    break;
		}
	    }
	    strcpy(nm,inname);
	    ans[1] = ans[0];
	    ans[0] = '_';
	    strncpy(ans+2, "col", 3);
    	    /*if not equal to NULL, find '.' then copy "_?col"*/
    	    if ( strrchr(nm,'.') )
	    {
		/*check if filename contains _?col already*/
		if ( strncmp( (strrchr(nm,'.')-5), ans ,5) )
		    strncpy( ans+5, (strrchr(nm,'.')), 5);
	    }
    	    /*if equal to NULL just add _3col to the end*/
    	    else strcat( ans, ".txt" );
	    
	    strcpy( (strrchr(nm,'.')), ans );
    	    printf("Output filename for %d column data: %s\n",f,nm);
    	    file_status(nm);
            ascii_write(nm, ft.ndp, f);
	}
	else if (md == 0) break;
    }
    return 0;
} /*END main()*/

/*==========================================================================*/
/*ascii_read: read ASCII format data	    	    	    	    	    */
/****************************************************************************/
int ascii_read(char * fname, int *col)
{
    float   rd = 0;
    int     chan = 0, i, res = 0, lchan = MAXPTS;
    FILE    *file1;
    	    
    /*open ascii file*/
    if ((file1= fopen(fname, "r" )) == NULL)
    {
    	printf("Cannot open file: %s \n", fname);
    	return -1;	    	    	    
    }
    
    /*zero spectrum array*/
    for (i = 0; i < MAXCOLS; i++)
    {        
    	ft.x[i] = 0.0;
    	ft.y[i] = 0.0;
    	ft.dy[i] = 0.0;
    }
   
    /*get number of columns in ascii file*/
    col_determ(file1, col);
    printf("Ascii %d Column format....\n", *col);    
    if (*col == 0)
    {
	printf("No numbers found in file %s; Exiting....", fname);
	return -1;
    }
	    
    for (chan = 0; chan < MAXPTS; chan++)
    {
	/*only for 1 column data set the x value equal to chan*/
	if (*col == 1)
	{
	    ft.x[chan] = (double)(chan);
	    if (chan == 0)
	    	printf("\t**** Spectrum will start from x=0 ****\n");
	    res = fscanf(file1, "%f", &rd);	
	    ft.y[chan] = (double)rd;
	    ft.dy[chan] = sqrt(ft.y[chan]);
	}
	else if (*col == 2)
	{
	    res = fscanf(file1, "%f", &rd);
	    ft.x[chan] = (double)rd;
	    res = fscanf(file1, "%f", &rd);
	    ft.y[chan] = (double)rd;
	    ft.dy[chan] = sqrt(ft.y[chan]);
	}
	else if (*col == 3)
	{
	    res = fscanf(file1, "%f", &rd);
	    ft.x[chan] = (double)rd;
	    res = fscanf(file1, "%f", &rd);
	    ft.y[chan] = (double)rd;
	    res = fscanf(file1, "%f", &rd);
	    ft.dy[chan] = (double)rd;
	}
	switch (res)
	{
	    case 0:
	    {
	    	printf("\nchannel= %d, ft.y[channel]= %f\n", chan,
		    	ft.y[chan]);
    	    	printf("Read error occurred for file: %s\n", fname);
		/*close ascii file*/
	    	fclose(file1);
	    	return -1;
	    }
	    case EOF:
	    {
	    	if (chan < MAXPTS)
		{
/*            	    printf("\nMax. no. of channels expected = %d\n", MAXPTS);*/
	    	    printf("    reached EOF after reading " 
		    	    "channel %d \n", chan );
		    if (chan == 0)
		    {
			printf("\n*******Incorrect file format*******\n");
			printf("....Exiting....\n\n");
			return -1;
		    }
/*            	    printf("-->.spe file will still be written......\n"); */
		}
		lchan = chan;
	    	chan = MAXPTS;
    	    	break;
	    }
	    default:
	    {
/*	    	if (chan == 0) printf("Reading ascii spectrum.......\n"); */
	    	break;
    	    }
    	}
    }
    chan = lchan;
    /*check for zeroes in data and ask about suppression*/
    chan = supp_zeroes(chan);
    /*close ascii file*/
    fclose(file1); 
    return (chan);
} /*END ascii_read()*/

/*==========================================================================*/
/* ascii_write: write an ASCII format spectrum	    	    	    	    */
/****************************************************************************/
void ascii_write(char name[], int numch, int col)
{
    int i;
    FILE *fasc;
    /*open .txt file*/
    if ( (fasc = fopen(name, "w" )) == NULL)
    {
        printf("Cannot open file: %s \n", name);
        return ; 
    }
   
    /*write .txt file*/
    for (i = 0; i < numch; i++)
    {
    	fprintf(fasc, "%8.1f\t%8.3f", ft.x[i], ft.y[i]);
	if (col == 3) fprintf(fasc, "\t%8.3f\n",ft.dy[i]);
	else fprintf(fasc, "\n");
    }
   
    printf(" written ==> %s, %d chs %d columns.\n", name, numch, col);
    fclose(fasc);
} /*END ascii_write()*/

/*==========================================================================*/
/* ask yes/no. 0 means no, 1 means yes	    	    	    	    	    */
/****************************************************************************/
int ask_yn(char *questn)
{
    char ans[3] = "";
    
    while (1)
    {
    	printf("%s\n", questn);
	get_ans(ans,1);
	if (ans[0] == 'y' || ans[0] == 'Y') return 1;
	else if (ans[0] == 'n' || ans[0] == 'N') return 0;
    }
} /*END ask_yn()*/

/*==========================================================================*/
/* chisq: calculate chisq   	    	    	    	    	    	    */
/****************************************************************************/
void chi2(double *chisq)
{
    static int	npars = NUMPARS;
    extern double fit, derivs[NUMPARS];
    double  diff = 0.0, dat;
    int     i, nip, ndf;
    
    *chisq = 0.0;
   /*no. independent parameters*/
    nip = npars - ft.nfp;
    /*no. degrees of freedom*/
    ndf = ft.ndp - nip;
    
    for (i = 0; i < ft.ndp; i++)
    {
    	eval(ft.pars, ft.x[i], &fit, derivs, 0);
    	diff = ft.y[i] - fit;
    	dat = ft.dy[i] * ft.dy[i];
    	if (dat == 0.f) dat = 1.f;
    	
    	*chisq += diff * diff / dat;
    }
    *chisq /= (double)ndf;
} /*END chisq()*/

/*==========================================================================*/
/* col_determ: find if a file has 1, 2, 3 or 4 column format      	    	    */
/****************************************************************************/
void col_determ(FILE *file1, int *col)
{
    int     i = 0, hash = 0;
    fpos_t  pos;
    
    /*skip comments lines starting with #*/
    skip_hash(file1);

    /*Try to decide if spectrum is 1, 2 or 3 column ascii format*/
    *col = 0;
    i = 0;
    
    /*store current file position in pos*/
    fgetpos(file1, &pos);
    
    while (1)
    {
    	/*read until first digit of first number*/
	while ( isdigit(hash = fgetc(file1)) == 0 )
    	    ;
	
	*col = 1;
    	/*found first digit of first number. read rest of digits*/
    	/*ignore decimal points*/
	while ( (isdigit(hash = fgetc(file1)) != 0) || (hash == '.') ) ;

    	/*push the last character back on the stream*/
    	ungetc(hash, file1);
    	/*now next test characters until end of line*/
    	while ( (hash = fgetc(file1)) != '\n' )
    	{
    	    /*if blank space, increment blank counter and continue*/
    	    if ( (hash == ' ') || (hash == '\t') ) i++;
	    
    	    /*if character is a digit*/
    	    /*check how many columns of numbers are present*/
    	    else if ( (isdigit(hash)) != 0 )
    	    {
    		/*if number is second column set 2 col. flag*/
    		if ( i > 0)
		{
		    /*brackets necessary round *col to increment what it
		    	points to not the address!*/
		    (*col)++;
/*		    printf("col. = %d\n", *col); */
		}
		
    		/*otherwise if > MAXCOLS columns set 1 col. flag and hope for best*/
    		if ( (i != 0) && (*col > MAXCOLS) )
    		{
    		    *col = 1;
    		    break;
    		}
    		/*reset blank counter to zero*/
    		i = 0;
    		/*read rest of digits*/
		while ( (isdigit(hash = fgetc(file1)) != 0) || (hash == '.') ) ;
		
    		/*push the last character back on the stream*/
    		ungetc(hash, file1);
    	    }
    	}
    	break;
    }	    
    /*reset file position to start*/
    fsetpos(file1, &pos);   
} /*END col_determ()*/

/*==========================================================================*/
/*compress2: compress data by factor of 2	     	    	    	    */
/****************************************************************************/
void compress2(char inname[], char outname[])
{
    static int	fact = 2;
    int    i;
    char    ext[10] = "", next[10] = "-", nname[80] = "";
     
    /*do compression*/
    for (i = 0; i < (int)(ft.ndp/2); i++)
    {
	ft.x[i] = (int)(ft.x[(2*i)]/2);
	ft.y[i] = ft.y[(2*i)] + ft.y[(2*i)+1];
	ft.dy[i] = sqrt( (ft.dy[(2*i)]*ft.dy[(2*i)]) +
		(ft.dy[(2*i)+1]*ft.dy[(2*i)+1]) );
    }
    printf("\t ....done. %d channels\n",i);
    ft.ndp = i;
    /*get output filename extension*/
    strcpy(nname, inname);
    if ( (strrchr(nname,'.')) ) strcpy( ext, (strrchr(nname,'.')) );
    /*now add the compression factor to the input and output filenames*/
    itoa(fact, next+1);
    strcat(next, ext);
    set_ext(nname, next);
    strcpy(outname, nname);
    set_ext(outname, ".fit");
    if (ask_yn(" Write compressed data to file? (y/n)"))
    {
   	printf("Output filename for compressed data: %s\n", nname);
    	file_status(nname);
        ascii_write(nname, i, 2);
    	strcpy(inname, nname);	
    }
    fact *= 2;
}/*END compress2()*/

/*==========================================================================*/
/*erfpol: returns good approximation to erfc(x)ex*x                         */
/****************************************************************************/
double erfpol(double x)
{
    double erfpol, t;
    
    t = (double)1.0/((double)1.0 + (double)0.47047*fabs(x));
    erfpol = t*((double)0.3480242 + t*((double)-0.0958798 + (double)0.7478556*t));
    return erfpol;
}/*END erfpol()*/

/*==========================================================================*/
/*erfpol_new: returns good approximation to erfc(x)ex*x                     */
/****************************************************************************/
double erfpol_new(double x)
{
    double erfpol, t;
    
    t = (double)1.0/((double)1.0 + (double)0.3275911*fabs(x));
    erfpol = t*((double)0.254829592 + t*((double)-0.284496736 +
            t*((double)1.421413741 + t*((double)-1.453152027 +
            (double)1.061405429*t))));
    return erfpol;
}/*END erfpol_new()*/

/*==========================================================================*/
/*eval: calculate fit using present values of the pars 	    	    	    */
/****************************************************************************/
int eval(double *pars, double x, double *fit, double *derivs, int mode)
{
/*  a function corresponding to a Gaussian convoluted with an
    exponential decay

     computes an exponential with decay time tau convoluted with
     a Gaussian with width sig i.e.
     convo = integral (-oo to 0 ) {exp(t/tau)*exp(-(x-t)^2/sig^2)}dt
     if tau is positiv, a left tail convolusion is calculated, if
     tau is negative, a right tail convolusion
     input parameter:
     x: independent variable, any value
     sig: width of the gaussian defined by exp (-x^2/sig^2), sig >= 0.
     tau: decay constant of exponential, exp(x/tau), any value
     return value
     one of the following functions
     sig = 0 and tau=0 convo = 0.
     sig = 0           convo = exp (x/tau)/tau
     tau = 0           convo = exp( -(-x/sig)^2) /(sqrt(pi)*sig)
     else              convo = above integral normalized to area=1
     F. Riess,  March 1987

     Calculate the fit using present values of the pars
       x is the channel number variable
       x1 is the channel number with respect to the prompt peak centroid
       pars[0] = the mean life*(-1)
       pars[1] = the sigma*sqrt(2) of the prompt peak
       pars[2] = the prompt peak centroid position
       pars[3] = the scaling factor for the calculated curve*/

    double  x1, xsq, xt, e, ep;
    double  d = 0.001, ofit = 0.0, npars[NUMPARS];
    int     i, lp = 0;
    
    for (i = 0; i < NUMPARS; i++)
    {
	npars[i] = pars[i];
	derivs[i] = 0.0;
    }
    while (lp < NUMPARS)
    {
    	*fit = 0.0;
    	x1 = x - pars[2];
    	if (pars[0] == 0.0)
    	{
	    if (pars[1] > 0.0)
	    {
	    	xsq = (x1/pars[1])*(x1/pars[1]);
	    	if (xsq < 15.0) *fit = exp(-xsq)/((double)1.7724539*pars[1]);
	    }
    	}
    	else if (pars[1] == 0.0)
    	{
	    xt = x1/pars[0];
	    if (xt > -20.0 && xt < 70.0) *fit = exp(xt)/fabs(pars[0]);
	    
    	}
    	else if (pars[1] > 0.0)
    	{
	    ep = ((double)0.5*pars[1])/pars[0];
	    xsq = (x1/pars[1]) + ep;
	    if (pars[0] < 0.0) xsq *= (double)-1.0;
	
	    xt = x1/pars[0];
	    ep *= ep;
	    e = xt + ep - xsq*xsq;
	    if (fabs(e) < 20.0) *fit = (double)0.5*erfpol_new(xsq)*exp(e);
	    
	    if (xsq < 0.0)
	    {
	    	e = xt + ep;
	    	if (e > -20.0 && e < 70.0) *fit = exp(e) - *fit;
		
	    }
	    *fit /= fabs(pars[0]);
    	}
    	/*derivs[3] can be calculated analytically*/
    	if (lp == 0) derivs[3] = *fit;
    	*fit *= pars[3];
	/*add offset incase of constant background; bkgnd is not fitted*/
	*fit += bkgnd;
	
    	/* calculate derivs only for mode.ge.1 */
    	if (mode >= 1)
	{
	    if (lp == 0)
	    {
		pars[0] += d;
		ofit = *fit;
	    }
	    else if (lp == 1)
	    {
		derivs[0] = (*fit-ofit)/d;
		pars[0] = npars[0];
		pars[1] += d;
	    }		
	    else if (lp == 2)
	    {
		derivs[1] = (*fit-ofit)/d;
		pars[1] = npars[1];
		pars[2] += d;
	    }		
	    else if (lp == 3)
	    {
		derivs[2] = (*fit-ofit)/d;
		pars[2] = npars[2];
	    }
	    lp++;
	}
	else return 0;
    }
    *fit = ofit;
    return 0;
} /*END eval()*/

/*==========================================================================*/
/* file_status: check file status   	    	    	    	    	    */
/****************************************************************************/
void file_status(char *fname)
{
    char ans[3] = "";
    struct stat statbuf;      /*need include files sys/stat.h and sys/types.h*/
    
    while ( (stat(fname, &statbuf) == 0) )
    {
    	printf("\n*****Output file %s exists. Overwrite (y/n)?\n", fname);
	get_ans(ans,1);
	if (ans[0] == 'y' || ans[0] == 'Y') break;
	else if (ans[0] == 'n' || ans[0] == 'N')
	{
	    printf("Enter new file name inc. extension (eg %s):\n",
		    strrchr(fname,'.'));
	    get_line(fname);
    	}
    }
} /*END file_status()*/

/*==========================================================================*/
/*fitter: fit parameters. vb == 0 suppresses output info except warnings    */
/****************************************************************************/
int fitter(int maxits,  double *chisq, int vb)
{
    static int npars = NUMPARS;

    extern double fit, derivs[NUMPARS];
    double ddat, alpha[NUMPARS][NUMPARS], array[NUMPARS][NUMPARS];
    double  r1, diff, beta[NUMPARS], b[NUMPARS], delta[NUMPARS];
    double  chisq1, flamda, dat, ers[NUMPARS];
    int    conv, nits, test, i, j, k, l, m, nextp[NUMPARS], ndf, nip;

    /* this subroutine is a modified version of 'CURFIT', in Bevington */
    /* see page 237 */
  
    *chisq = 0.0;
    
    for (i = 0; i < npars; i++)
    {
    	derivs[i] = 0.0;
    	for (j = 0; j < npars; j++) array[i][j] = (double)0.0;
    }
    fit = 0.0;
  
    /*no. independent parameters*/
    nip = npars - ft.nfp;
    /*no. degrees of freedom*/
    ndf = ft.ndp - nip;
    if (ndf < 1)
    {
    	printf("No Degrees-of-Freedom.\n");
    	return 3;
    }
    if (nip < 2)
    {
    	printf("Too many fixed parameters.\n");
    	return 4;
    }
    k = 0;
    for (j = 0; j < npars; ++j)
    {
    	if (ft.freepars[j]) nextp[k++] = j;
    }
    if (k != nip)
    {
    	printf("No. of independent pars != sum(freepars)!!!\n");
    	return 5;
    }
    flamda = 0.001f;
    nits = 0;
    test = 0;
    derivs[0] = 1.0f;
    for (i = 0; i < npars; ++i)
    {
    	ft.errs[i] = 0.0f;
    	b[i] = ft.pars[i];
    }
    /* evaluate fit, alpha & beta matrices, & chisq */
    NEXT_ITERATION:
    for (j = 0; j < nip; ++j)
    {
    	beta[j] = (double)0.0;
    	for (k = 0; k <= j; ++k) alpha[k][j] = 0.0;
    }
    chisq1 = 0.0f;
    for (i = 0; i < ft.ndp; ++i)
    {
    	eval(ft.pars, ft.x[i], &fit, derivs, 1);
    	diff = ft.y[i] - fit;
    	dat = ft.dy[i] * ft.dy[i];
    	if (dat < 1.0f) dat = 1.0f;
	
    	ddat = (double)dat;
    	chisq1 += diff * diff / dat;
    	for (l = 0; l < nip; ++l)
	{
      	    j = nextp[l];
      	    beta[l] += diff * derivs[j] / dat;
      	    for (m = 0; m <= l; ++m)
	    {
	    	alpha[m][l] += (double) derivs[j] * (double) derivs[nextp[m]] / ddat;
      	    }
    	}
    }
    chisq1 /= (double)ndf;
    
    /* invert modified curvature matrix to find new parameters */
    INVERT_MATRIX:
    array[0][0] = flamda + 1.0f;
    for (j = 1; j < nip; ++j)
    {
    	for (k = 0; k < j; ++k)
	{
      	    if (alpha[j][j] * alpha[k][k] == 0.0)
	    {
	    	printf("Cannot - diag. element %i or %i eq. to 0.0\n", j, k);
	    	return 1;
      	    }
      	    array[k][j] = alpha[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
      	    array[j][k] = array[k][j];
    	}
    	array[j][j] = flamda + 1.0f;
    }
    matinv(array[0], nip, npars);
    if (!test)
    {
    	for (j = 0; j < nip; ++j)
	{
      	    if (alpha[j][j] * alpha[j][j] == 0.0)
	    {
	    	printf("Cannot - diag. element %i eq. to 0.0\n", j);
	    	return 1;
      	    }
      	    delta[j] = 0.f;
      	    for (k = 0; k < nip; ++k)
	    {
	    	delta[j] += beta[k] * array[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
      	    }
    	}
    	/* if chisq increased, increase flamda & try again */
    	*chisq = 0.0f;
    	for (l = 0; l < nip; ++l)
	{
      	    j = nextp[l];
      	    b[j] = ft.pars[j] + delta[l];
    	}
    	for (i = 0; i < ft.ndp; ++i)
	{
      	    eval(b, ft.x[i], &fit, derivs, 0);
      	    diff = ft.y[i] - fit;
      	    dat = ft.dy[i] * ft.dy[i];
      	    if (dat == 0.f) dat = 1.f;
      	    
	    *chisq += diff * diff / dat;
    	}
    	*chisq /= (double)ndf;
    	if (*chisq > chisq1 && flamda < 2.f)
	{
      	    flamda *= 10.f;
      	    goto INVERT_MATRIX;
    	}
    }
    
    /* evaluate parameters and errors */
    /* test for convergence */
    conv = 1;
    for (j = 0; j < nip; ++j)
    {
    	if (array[j][j] < 0.0) array[j][j] = 0.0;
	
    	ers[j] = sqrt(array[j][j] /  alpha[j][j]) * sqrt(flamda + 1.0f);
    	if ((r1 = delta[j], fabs(r1)) >= ers[j] / 1e3f) conv = 0;	
    }
    if (!test)
    {
    	for (j = 0; j < npars; ++j)
	{
      	    ft.pars[j] = b[j];
    	}
    	flamda /= 10.f;
    	++nits;
    	if (! conv && nits < maxits)
	{
	    goto NEXT_ITERATION;
	}

    	/* re-do matrix inversion with flamda=0 to calculate errors */
    	flamda = 0.f;
    	test = 1;
    	goto INVERT_MATRIX;
    }

    /* list data and exit */
    for (l = 0; l < nip; ++l) ft.errs[nextp[l]] = ers[l];
    
    if (vb) printf(" %i indept. pars    %i degrees of freedom\n", nip, ndf);
    if (conv)
    {
    	if (vb) printf(" %i iterations,  Chisq/D.O.F. = %.3f\n", nits, *chisq);
    	return 0;
    }
    printf(" Failed to converge after %i iterations,  Chisq/D.O.F. = %.3f\n"
    	"  WARNING - do not believe quoted errors.\n", nits, *chisq);
    return 2;
} /*END fitter()*/	

/*==========================================================================*/
/* get_ans: get answer without waiting for carriage return   	    	    */
/****************************************************************************/
void get_ans(char ans[], int num)
{
    int     i;
    struct  termios newt, oldt;
    
    while (1)
    {
    	tcgetattr(0, &oldt);
	newt = oldt;
    	newt.c_lflag &= ~ICANON;
    	newt.c_cc[VMIN] = 1;
    	newt.c_cc[VTIME] = 0;
	/*handle sigs*/
/*    	newt.c_lflag |= ISIG;*/
    	tcsetattr(0, TCSANOW, &newt);
    	i = 0;
    	while( (ans[i++] = (char)getchar()) != '\n' && i < num) ;
    	
	tcsetattr(0, TCSANOW, &oldt);
    	if (ans[i-1] != '\n') printf("\n");
	else if (ans[0] == '\n') continue;
	
    	ans[i] = '\0';
    	return ;
    }
} /*END get_ans()*/

/*==========================================================================*/
/* get_line: read a line into s, return a length    	    	    	    */
/****************************************************************************/
void get_line(char ans[])
{
    int     c, i = 0;

    memset(ans,'\0',sizeof(ans));
    /*note that scanf() leaves a carriage return in the keyboard buffer*/
    while ( (c = getchar()) != '\n' && c != EOF && i < 80) ans[i++] = c;
    
    ans[i] = '\0';
} /*END get_line()*/

/*==========================================================================*/
/* get_mode: get mode from user input	    	    	    	    	    */
/****************************************************************************/
int get_mode(int md)
{
    char    ans[10] = "";
    
    while(1)
    {
    	printf("\n 1) Enter filename and read data\n");
    	printf(" 2) Enter initial parameters\n");
    	printf(" 3) Fix or free parameters\n");
    	printf(" 4) Perform fit and write output to file\n");
    	printf(" 5) Compress data by factor of 2\n");
    	printf(" 6) Explore background value as function of chi^2\n");
    	printf(" 7) Print current values of parameters to screen\n");
    	printf(" 8) Write output with current values of parameters\n");
    	printf(" 9) Write input to file as 3 column data (x y dy)\n");
    	printf(" 0) Quit\n");
	
	get_ans(ans,1);
	
 	if (ans[0] >= '0' && ans[0] <= '9')
	{
	    md = ans[0] - '0';
	    break;
	}
    }
    return md;
} /*END get_mode()*/

/*==========================================================================*/
/* get_pars: extract comma or space separated numbers from string ans0	    */
/****************************************************************************/
void get_pars(char ans0[], double pars[], int num)
{
    int     i, j = 0, k, minus, p = 0;
    char    ans1[80] = "";

/*    get_line(ans0);*/
    
    for (i = 0; i < num; i++)
    {
	k = 0;
	minus = 0;
	memset(&ans1,'\0',sizeof(ans1));
	while ( ! isdigit(ans0[j]) )
	{
	    if (ans0[j] == '-') minus = 1;
	    j++;
	}
	if (minus == 1)
	{
	    j -= 1;
	    ans0[j] = '-';
	}	
	while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
		&& j < strlen(ans0) ) ans1[k++] = ans0[j++];
	
	if (k > 0) p++;
	j++;
	pars[i] = (double)atof(ans1);
    }
    printf("\n");
} /*END get_pars()*/

/*===========================================================================*/
/* convert integer n to string */
/*****************************************************************************/
void itoa(int n, char s[])
{
    int i = 0, sign;
    
    if ((sign = n) < 0)     	    	    	    	    	/*record sign*/
    n = -n;		    	    	    	    	    /*make n positive*/
	
    do {                    	    	/* generates digits in reverse order */
    	s[i++] = n % 10 + '0';      	    	    	     /*get next digit*/
    } while ((n /= 10) > 0);    	    	    	    	 /* delete it*/
    if (sign < 0)
    {
    	s[i++] = '-';
    }
    s[i] = '\0';
    reverse(s);
} /*END itoa*/

/*==========================================================================*/
/*matinv: invert matrix     	    	    	    	    	    	    */
/****************************************************************************/
int matinv(double *array, int norder, int dim)
{
    double amax, save, d1;
    int i, j, k, ik[1000], jk[1000];

    for (k = 0; k < norder; ++k)
    {
    	/* find largest element array(i,j) in rest of matrix */
    	amax = 0.f;
    	while (1)
	{
      	    for (i = k; i < norder; ++i)
	    {
	    	for (j = k; j < norder; ++j)
		{
	    	    if (fabs(amax) - (d1 = array[i + j*dim], fabs(d1)) <= 0.f)
		    {
	    	    	amax = array[i + j*dim];
	    	    	ik[k] = i;
	    	    	jk[k] = j;
	    	    }
	    	}
      	    }
      	    if (amax == 0.f) return 0;
      	    /* interchange rows and columns to put amax in array(k,k) */
      	    i = ik[k];
      	    if (i < k) continue;
      	    if (i > k)
	    {
	    	for (j = 0; j < norder; ++j)
		{
	    	    save = array[k + j*dim];
	    	    array[k + j*dim] = array[i + j*dim];
	    	    array[i + j*dim] = -save;
	    	}
      	    }
      	    j = jk[k];
      	    if (j >= k) break;
    	}
    	if (j > k)
	{
      	    for (i = 0; i < norder; ++i)
	    {
	    	save = array[i + k*dim];
	    	array[i + k*dim] = array[i + j*dim];
	    	array[i + j*dim] = -save;
      	    }
    	}
    	/* accumulate elements of inverse matrix */
    	for (i = 0; i < norder; ++i)
	{
      	    if (i != k) array[i + k*dim] = -array[i + k*dim] / amax;
    	}
    	for (i = 0; i < norder; ++i)
	{
      	    for (j = 0; j < norder; ++j)
	    {
	    	if (i != k && j != k)
	    	    array[i + j*dim] += array[i + k*dim] * array[k + j*dim];
      	    }
    	}
    	for (j = 0; j < norder; ++j)
	{
      	    if (j != k) array[k + j*dim] /= amax;
    	}
    	array[k + k*dim] = 1.f / amax;
    }
    /* restore ordering of matrix */
    for (k = norder-1; k >= 0; --k)
    {
    	j = ik[k];
    	if (j > k)
	{
      	    for (i = 0; i < norder; ++i)
	    {
	    	save = array[i + k*dim];
	    	array[i + k*dim] = -array[i + j*dim];
	    	array[i + j*dim] = save;
      	    }
    	}
    	i = jk[k];
    	if (i > k)
	{
      	    for (j = 0; j < norder; ++j)
	    {
	    	save = array[k + j*dim];
	    	array[k + j*dim] = -array[i + j*dim];
	    	array[i + j*dim] = save;
      	    }
    	}
    }
    return 0;
} /*END matinv()*/

/*==========================================================================*/
/*read_par: read user input	    	    	    	    	    	    */
/****************************************************************************/
void read_par(double pars[], int offset, int num, char s[])
{
    int     i = 0;
    char    ans[80] = "";
    
    printf("%s",s);
    printf(" [<Enter> for ");
    for (i = 0; i < num; i++) printf("%f ",pars[offset+i]);
    printf("]\n");
    get_line(ans);
    if (strlen(ans) == (int)0 && (ans[0] == ((char)0))) ;
    else get_pars(ans, pars+offset, num);
}/*END read_par()*/

/*==========================================================================*/
/* reverse: reverse string s in place	    	    	    	    	    */
/****************************************************************************/
void reverse(char s[])
{
	int c, i, j;
	
	for (i = 0, j = strlen(s)-1; i < j; i++, j--)
	{
	    c = s[i];
	    s[i] = s[j];
	    s[j] = c;

	}
} /*END reverse()*/

/*==========================================================================*/
/*scaler: workout approximate scaling factor	    	    	    	    */
/****************************************************************************/
void scaler(double *fit, double *derivs, int mode)
{
    double  step, dsqr, dnsqr, sign, osign;
    double  exp = 0.0, ex, sum = 0.0, thr = 0.0, st[16384][2];
    int     i, isl, isu, lp = 1, cnt = 0;
    
    for (i = 0; i < 16384; i++)
    {
	st[i][0] = 0.0;
	st[i][1] = 0.0;
    }
    
    isl = ft.x[0];
    isu = ft.x[ft.ndp - 1];
    for (i = 0; i < ft.ndp; i++)
    {
	ex = ft.x[i];
	exp = ft.y[i];
	sum += exp;
	st[(int)ex][0] = exp;
    	eval(ft.pars, ft.x[i], fit, derivs, 0);
	thr = *fit;
	st[(int)ex][1] = thr;
    }
    /*for mode=1 just calculate curve, don't attempt a fit*/
    if (mode == 1) return ;
    /*sum is total no. of experimental counts*/
    ft.pars[3] = sum;
    step = sum/10.0;
    sign = 0.0;
    while ( (step > 1.0 || lp == 1) && lp < 50)
    {
	cnt++;
	lp++; 
	osign = -1.0*sign;
	dsqr = 0.0;
	dnsqr = 0.0;
	for (i = isl; i <= isu; i++)
	{
	    exp = st[i][0];
	    thr = st[i][1];
	    dsqr += (exp - (thr*ft.pars[3]))*(exp - (thr*ft.pars[3]));
	    dnsqr += (exp-(thr*(ft.pars[3]+1)))*(exp-(thr*(ft.pars[3]+1)));
	}
	if (dnsqr > dsqr) sign = 1.0;
	else sign = -1.0;
	
	if (sign == osign) step *= 0.4;
	
	ft.pars[3] -= sign*step;
    }
    ft.pars[3] += sign*step; 
    printf("  %d iterations for scaling value\n",lp); 
} /*END scaler()*/

/*==========================================================================*/
/* set_ext: set file extension of string fname[] to ext[]   	    	    */
/****************************************************************************/
void set_ext(char fname[], char ext[])
{    
    /*if not equal to NULL, find '.' then copy ext*/
    if ( (strrchr(fname,'.')) ) strcpy( (strrchr(fname,'.')), ext );
    /*if equal to NULL just add ext to the end*/
    else strcat( fname, ext );
} /*END set_ext()*/

/*==========================================================================*/
/* skip_hash: skip comment lines starting with hash (#) at start of file    */
/****************************************************************************/
void skip_hash(FILE *file)
{
    int     i = 0, hash = 0;
    fpos_t  pos;
     	   
    for (i = 0; i < 1000; i++)
    {
    	/*store position corresponding to start of line*/
    	fgetpos(file, &pos);
    	if ( (hash = fgetc(file)) != '#' )
    	{
    	    /*reset position to start of line*/
    	    fsetpos(file, &pos);
    	    break;
    	}
    	/*read rest of line*/
    	while ( ( (hash = fgetc(file)) != '\n' ) && (hash != EOF) )
    	    ;
    }
} /*END skip_hash()*/

/*==========================================================================*/
/*supp_zeroes: check for zeroes y value for zeroes and optionally suppress  */
/****************************************************************************/
int supp_zeroes(int numch)
{
    int     i = 0, j = 0;
    char    ans[4] = "";
    
    /*check y value of data for zeros*/
    while (i < numch && ft.y[i++] != 0.0f) ;

    if (i == numch) return numch;
    
    while (1)
    {
    	printf("\n**** Zeroes found in y data: Suppress (y/n)?\n");
	get_ans(ans,1);
	if (ans[0] == 'y' || ans[0] == 'Y')
	{
    	    /*now remove zeroes*/
    	    for (i = 0; (i < numch) && ((i + j) < numch); i++)
    	    {
	    	while ( (ft.y[i+j] == 0.0f) && ((i + j) < numch) ) j++;
		
		    ft.x[i] = ft.x[i+j];
		    ft.y[i] = ft.y[i+j];
		    ft.dy[i] = ft.dy[i+j];
    	    }   
    	    if ( ((i + j) >= numch) && (ft.y[i+j] == 0.0f) ) i--;
	
	    return i;
	}
	else if (ans[0] == 'n' || ans[0] == 'N') return numch;
    }
} /*END supp_zeroes()*/

/*==========================================================================*/
/* write_output: write output to file and coefficients to screen    	    */
/****************************************************************************/
void write_output(char inname[], char outname[], double init[], double x[],
	double *fit, double chisqr, int mode)
{
    int     i, f;
    FILE    *fout = 0;
    
    /*open output file*/
    if (mode != 1)
    {
    	if ((fout = fopen(outname, "w" )) == NULL)
    	{
    	    printf("Cannot open file: %s \n", outname);
    	    return ;
	}			
    }
    /*write each initial guesses to output file*/    
    if (ft.pars[0] < 0.0) f = -1;
    else f = 1;
    printf("\n++++++++++++Results of fit++++++++++++\n");
    printf(" t_1/2 = %.3f",f*ft.pars[0]*log(2.0));
    if (mode != 2) printf(" %s",(ft.freepars[0]) ? "+/- " : "FIXED\n");
    if (mode != 2 && ft.freepars[0]) printf("%.3f\n",ft.errs[0]*log(2.0));
    if (mode == 2 && ft.freepars[0]) printf(" (no fit)\n");
    
    printf(" FWHM = %.3f",ft.pars[1]/0.6005612);
    if (mode != 2) printf(" %s",(ft.freepars[1]) ? "+/- " : "FIXED\n");   
    if (mode != 2 && ft.freepars[1]) printf("%.3f\n",ft.errs[1]/0.6005612);
    if (mode == 2 && ft.freepars[1]) printf(" (no fit)\n");
    
    printf(" Centroid = %.2f",ft.pars[2]);
    if (mode != 2) printf(" %s",(ft.freepars[2]) ? "+/- " : "FIXED\n");   
    if (mode != 2 && ft.freepars[2]) printf("%.2f\n",ft.errs[2]);
    if (mode == 2 && ft.freepars[2]) printf(" (no fit)\n");
    
    printf(" Scaling factor = %.2f",ft.pars[3]);
    if (mode != 2) printf(" %s",(ft.freepars[3]) ? "+/- " : "FIXED\n");   
    if (mode != 2 && ft.freepars[3]) printf("%.2f\n",ft.errs[3]);
    if (mode == 2 && ft.freepars[3]) printf(" (no fit)\n");
    
    printf(" Background level = %.2f (not fitted)\n",bkgnd);
    printf("++++++++++++++++++++++++++++++++++++++\n");
    /*for mode == 1 just write out result of fit don't write to file*/
    if (mode == 1) return ;
    
    fprintf(fout, "# Input file: %s\n",inname);
    fprintf(fout, "# Output file: %s\n",outname);
    fprintf(fout, "# Initial guess t_1/2 = %.2f\n",init[0]);
    fprintf(fout, "# Initial guess FWHM = %.2f\n",init[1]);
    fprintf(fout, "# Initial guess centroid = %.2f\n",init[2]);
    fprintf(fout, "# Initial scaling factor = %.2f\n",init[3]);
    fprintf(fout, "#++++++++++++Results of fit++++++++++++\n");
    fprintf(fout, "# t_1/2 = %.3f %s",f*ft.pars[0]*log(2.0),
	    (ft.freepars[0]) ? "+/- " : "FIXED\n");
    if (mode != 2 && ft.freepars[0])
	fprintf(fout, "%.3f\n",ft.errs[0]*log(2.0));
    fprintf(fout, "# FWHM = %.3f %s",ft.pars[1]/0.6005612,
	    (ft.freepars[1]) ? "+/- " : "FIXED\n");
    if (mode != 2 && ft.freepars[1])
	fprintf(fout, "%.3f\n",ft.errs[1]/0.6005612);
    fprintf(fout, "# Centroid = %.2f %s",ft.pars[2],
	    (ft.freepars[2]) ? "+/- " : "FIXED\n");
    if (mode != 2 && ft.freepars[2]) fprintf(fout, "%.2f\n",ft.errs[2]);
    fprintf(fout, "# Scaling factor = %.2f %s",ft.pars[3],
	    (ft.freepars[3]) ? "+/- " : "FIXED\n");
    if (mode != 2 && ft.freepars[3]) fprintf(fout, "%.2f\n",ft.errs[3]);
    fprintf(fout, "# Background level = %.2f (not fitted)\n",bkgnd);
    fprintf(fout, "# Chisqr/D.O.F. = %.3f\n",chisqr);
    fprintf(fout, "#++++++++++++++++++++++++++++++++++++++\n");   
    fprintf(fout, "#   Chan \t    fit\n");
    
    /*write the fitted points to output file*/        
    for (i = 0; i < ft.ndp; i++)
    {
/*	printf("x[%d]=%f",i,x[i]);*/
    	eval(ft.pars, x[i], fit, derivs, 0);
    	fprintf(fout, "%8.1f\t%8.3f\n", x[i],*fit);    
    }
    fclose(fout);
    printf(" --> Output written to file: %s\n",outname);
} /*END write_output()*/

