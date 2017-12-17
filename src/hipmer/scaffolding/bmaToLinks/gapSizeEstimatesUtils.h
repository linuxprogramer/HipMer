#include <math.h>
#define PI 3.14159265359

/* Util functions for gap size estimations */

double erf(double x)
{
   double absX = (x < 0.0) ? -x : x;
   double t = 1.0/(1.0+0.5*absX);
   
   double t2 = t*t;
   double t3 = t*t2;
   double t4 = t*t3;
   double t5 = t*t4;
   double t6 = t*t5;
   double t7 = t*t6;
   double t8 = t*t7;
   double t9 = t*t8;
   
   double a0 = -1.26551223;
   double a1 = 1.00002368;
   double a2 = 0.37409196;
   double a3 = 0.09678418;
   double a4 = -0.18628806;
   double a5 = 0.27886807;
   double a6 = -1.13520398;
   double a7 = 1.48851587;
   double a8 = -0.82215223;
   double a9 = 0.17087277;
   
   double tau = t*exp(-x*x + a0 + a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t5 + a6*t6 + a7*t7 + a8*t8 + a9*t9);
   
   if (x < 0.0) {
      return (tau-1.0);
   } else {
      return (1.0-tau);
   }
}

double G0(double a, double b, double mu, double sigma)
{
   double rt2sig = sqrt(2)*sigma;
   double erfa = erf((a-mu)/rt2sig);
   double erfb = erf((b-mu)/rt2sig);
   
   return ((sqrt(PI)/sqrt(2))*sigma*(erfb-erfa));
}

double G1(double a, double b, double mu, double sigma)
{
   
   double za = (a-mu)/(1.0*sigma);
   double zb = (b-mu)/(1.0*sigma);
   
   double expa = exp(-0.5*za*za);
   double expb = exp(-0.5*zb*zb);
   
   double g0 = G0(a,b,mu,sigma);
   
   return (sigma*sigma*(expa-expb) + mu*g0);
}

double G2(double a, double b, double mu, double sigma)
{
   double za = (a-mu)/(1.0*sigma);
   double zb = (b-mu)/(1.0*sigma);
   
   double expa = exp(-0.5*za*za);
   double expb = exp(-0.5*zb*zb);
   
   double g0 = G0(a,b,mu,sigma);
   double g1 = G1(a,b,mu,sigma);
   
   double sigma2 = sigma*sigma;
   
   return (sigma2*g0 + mu*g1 + sigma2*(a*expa - b*expb));
}

double meanSpanningClone(double g, int k, int l, int c1, int c2, double mu, double sigma)
{
   double x1 = g+2*k-1;
   double x2 = g+c1+l;
   double alpha = x2-x1;
   double x3 = g+c2+l;
   double x4 = x3+alpha;
   double num = 0.0;
   double den = 0.0;
   double N1 = G2(x1,x2,mu,sigma)-x1*G1(x1,x2,mu,sigma);
   double N2 = (x2-x1)*G1(x2,x3,mu,sigma);
   double N3 = x4*G1(x3,x4,mu,sigma)-G2(x3,x4,mu,sigma);
   double D1 = G1(x1,x2,mu,sigma)-x1*G0(x1,x2,mu,sigma);
   double D2 = (x2-x1)*G0(x2,x3,mu,sigma);
   double D3 = x4*G0(x3,x4,mu,sigma)-G1(x3,x4,mu,sigma);
   
   num = N1+N2+N3;
   den = D1+D2+D3;
   
   if (den > 0.0) {
      return (1.0*num/(1.0*den));
   } else {
      return 0.0;
   }
}

double estimateGapSize(double meanAnchor ,int k, int l, int c1, int c2, double mu, double sigma)
{
   double gMax = mu+3.0*sigma-2.0*k;
   double gMin = -1.0*(k-2);
   double gMid = 1.0*mu-meanAnchor;
   double aMax, aMin, aMid, deltaG, iterations;
   
   if (gMid < gMin) {
      gMid = gMin+1;
   } else if (gMid > gMax) {
      gMid = gMax-1;
   }
   
   aMax = meanSpanningClone(gMax,k,l,c1,c2,mu,sigma) - gMax;
   aMin = meanSpanningClone(gMin,k,l,c1,c2,mu,sigma) - gMin;
   aMid = meanSpanningClone(gMid,k,l,c1,c2,mu,sigma) - gMid;
   deltaG = gMax-gMin;
   iterations = 0;
   
   while (deltaG > 10.0) {
      iterations++;
      if (meanAnchor > aMid) {
         gMax = gMid;
         aMax = aMid;
         gMid = (gMid+gMin)/2;
         aMid = meanSpanningClone(gMid,k,l,c1,c2,mu,sigma) - gMid;
      } else if (meanAnchor < aMid) {
         gMin = gMid;
         aMin = aMid;
         gMid = (gMid+gMax)/2;
         aMid = meanSpanningClone(gMid,k,l,c1,c2,mu,sigma) - gMid;
      } else {
         break ;
      }
      deltaG = gMax-gMin;
   }
   
   return gMid;
}
