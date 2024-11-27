#include <stdio.h>
#include <math.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


// gcc -lm test_j.c -o test_j

double function_j(double f, double fp, double fptilde) {
 
//   double a  = 0.0081;
//   double b  = 0.6;
//   double g  = 9.807;
//   double pi = 4.*atan(1.);

//   fprintf(stdout,"pi = %lg\n",pi);

//   double fptildemin = (1.0/2.0/pi) * pow((4.0 * b / 5.0), (1.0/4.0));

//   fprintf(stdout,"fptildemin = %lg\n",fptildemin);

   const double fptildemin = 0.132474;
   const double aX = 0.675014;
   const double gX = 0.875572;
   const double cX = 0.000641624 * 96.1772;

//   fprintf(stderr,"cX = %lg\n",cX);

   double gC = 5.87;
   double aC = 0.0317;

   double saC = 0.0547;
   double saX = 0.32;

   double sbC = 0.0783;
   double sbX = 0.16;


   double fpt = MAX(fptilde, fptildemin);
   double logF = log(fpt);
   double alpha   = aC  * exp(aX*logF);
   double gamma   = gC  * exp(gX*logF);
   double sigma_a = saC * exp(saX*logF);
   double sigma_b = sbC * exp(sbX*logF);
   double sigma   = sigma_a;
   if ( f > fp ) sigma = sigma_b;   
   double tmp = fp/f;
   tmp = tmp * tmp;
   tmp = tmp * tmp;
   
   double exp1arg = -1.25 * tmp;

   tmp = (f-fp)/(sigma*fp);
   double exp2arg = -0.5 * tmp * tmp;
   
   double one_over_f = 1./f;
   tmp = one_over_f * one_over_f;
   tmp = tmp * tmp * one_over_f;
   double S = alpha * cX * tmp * exp(exp1arg) * pow(gamma, exp(exp2arg));

   return S;
}


int main( int argc, char **argv) {

  double S, f, fp, fptilde;
  int f_i,fp_i,fptilde_i;
#pragma omp parallel for private(f_i,fp_i,fptilde_i,f,fp,fptilde) 
  for (f_i=0;f_i<1001;++f_i)
  {
      f = 0.01 * f_i - 5.0;
      for (fp_i=0;fp_i<1001;++fp_i)
      {
          fp = fp_i * 0.01;
          for (fptilde_i=0;fptilde_i<1001;++fptilde_i)
          {
              fptilde = 0.01 * fptilde_i;
              S = function_j(f, fp, fptilde);
          }    
      }
  }
/*
  for (f = -5.; f <= 5.; f += 0.01) {
    for (fp = 0.; fp <= 10.; fp += 0.01) {
      for (fptilde = 0.; fptilde <= 10.; fptilde += 0.01) {
        S = function_j(f, fp, fptilde);
      }
    }
  }
*/
 f = 1.5;
 fp = 3.04;
 fptilde = 4.4;
 S = function_j(f,fp,fptilde);
 fprintf(stdout,"S = %lg\n",S);
}

