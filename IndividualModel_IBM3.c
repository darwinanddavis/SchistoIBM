/* file IndividualModel_IBM3.c */
#include <R.h>
#include <math.h>
static double parms[29];

#define iM parms[0]
#define k parms[1]
#define M parms[2]
#define EM parms[3]
#define Fh parms[4]
#define muD parms[5]
#define DR parms[6]
#define yRP parms[7]
#define ph parms[8]
#define yPE parms[9]
#define iPM parms[10]
#define eh parms[11]
#define mP parms[12]
#define alpha parms[13]
#define yEF parms[14]
#define LM parms[15]
#define kR parms[16]
#define d0 parms[17]
#define kk parms[18]
#define hb parms[19]
#define theta parms[20]
#define mR parms[21]
#define yVE parms[22]
#define ENV parms[23]
#define Lp parms[24]
#define SAtotal parms[25]
#define r parms [26]
#define K parms [27]
#define Det parms [28]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=29;
odeparms(&N, parms);
}

/* Derivatives and 2 output variables */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");

double Chi = M/(1 + EM);
double fH = y[0]/(y[0]+Fh);
double fP = y[2]/(y[2] + eh);

double L = y[1];
double LG = fmax(y[1], Lp);

double GVOL = pow(LG, 3);
double VOL = pow(L,3);
double SA = pow(L,2);

double Dens = y[5]/(Chi*GVOL);
double kstar = fmin(k + y[5]*alpha, 1);
double g = 1/(yVE*kstar*EM);
double aM = iM*yEF;
double mV = aM*k/(LM*Chi);
double mD = muD*mV;
double rp = Dens*Dens/(ph*ph + Dens*Dens);
double Jec = y[2]*g/(g + y[2])*(aM*SA + yVE*EM*(mV+mR*EM*y[7])*Chi*VOL);

ydot[0] = -(iM*SAtotal/ENV)*fH + r*y[0]*(1 - y[0]/K) + Det;
ydot[1] = yVE/(3*Chi*SA)*(kstar*Jec - (mV+mR*EM*y[7])*Chi*VOL);
ydot[2] = aM/(Chi*EM*L)*(fH - y[2]) - iPM*y[5]*fP/(EM*Chi*VOL);
ydot[5] = yPE*iPM*fP*(1 - rp)*y[5] - mP*y[5];
ydot[6] = fmax(yRP*yPE*iPM*fP*rp*y[5],0);
if(y[3] < DR){
  ydot[3] = (1 - kstar)*Jec - mD*y[3];
  ydot[4] = 0;}else{
  ydot[3] = fmin(0, (1 - kstar)*Jec - mD*DR);
  ydot[4] = fmax((1 - kstar)*Jec - mD*DR, 0);}
ydot[7] = theta/(Chi*VOL)*ydot[6] + kR*(1-y[2]) - kR*y[7] - 3*y[7]*ydot[1]/L;
ydot[8] = kk*fmax(y[7] - d0, 0) + hb;
if(y[2] <= 0){
  ydot[0] = 0;
  ydot[1] = 0;
  ydot[2] = 0;
  ydot[3] = 0;
  ydot[4] = 0;
  ydot[5] = 0;
  ydot[6] = 0;
  ydot[7] = 0;
  ydot[8] = 0;
  y[8] = 10;}

  yout[0] = exp(-y[8]);
  yout[1] = LG;

}

/* END file IndividualModel_IBM3.c */ 
