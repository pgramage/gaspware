#include<stdio.h>
#include<math.h>
#define ESP 0.00000000001
#define RANGE 200





long double F(long double x)
{
return 2.159-13.589+(0.917948-0.837289)*x+(-0.756290E-05-0.169887E-03)*x*x + 0.120416E-06*x*x*x;
}



/*
double F(int deg, double *pol1, double *pol2)
  {
    double ans=0; 
    int i;
    for (i=0; i<deg; i++) ans+=(pol1[i]-pol2[i])
    return ans;
    
  }
*/  
  

int main(int argc, char **argv)
{
/*  
  if (argc < 2) {
    printf("ERROR:  .mcal file required as argument:\n\n");
    exit(0);
  }
  
    
  freopen(argv[1], "r", stdin);
*/  
  
  
  long double x1,x2,x3,f1,f2,t;
  printf("\nEnter the value of x1: ");
  scanf("%Lf",&x1);
  printf("\nEnter the value of x2: ");
  scanf("%Lf",&x2);
  printf("\n______________________________________________\n");
  printf("\n    x1\t  x2\t  x3\t     f(x1)\t f(x2)");
  printf("\n______________________________________________\n");
  do
  {
  f1=F(x1);
  f2=F(x2);
  x3=x2-((f2*(x2-x1))/(f2-f1));
  printf("\n%Lf   %Lf   %Lf   %Lf   %Lf",x1,x2,x3,f1,f2);
  x1=x2;
  x2=x3;
  if(f2<0)
    t=-f2;
  else
    t=f2;
  }while(t>ESP);
printf("\n______________________________________________\n");
printf("\n\nApp.root = %Lf\n\n",x3);


  

}