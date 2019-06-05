//  Toy model for the estimation of r_sp including nDGP gravity: 1806.04302

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

const double pi=4.0*atan(1.0);	//	pi number

double rta,rvir,rsp,avir,ata,acol,asp;
double Gama,di,Ri,a,rc;


// all masses in 10^11 Ms/h units
const double Om	= 3.0e-1;
const double Ol = 1.0e0 - Om;
const double h = 7.e-1;

//const double Mpc = 3.08567758149 *pow(10, 22);
//const double GMsc2 = 2953.25024 *0.5 /Mpc *pow(10, 11);	//	G Msun c^2 /2 e-11 [Mpc]
const double Mpc = 3.08567758149;			//	Mpc e-22 [m]
const double GMsc2 = 2953.25024e-11 *0.5 /Mpc;		//	G Msun c^2 /2 e-11 [Mpc]
const double c = 299792.458;
const double GMs = GMsc2 *c *c;				//	G Msun /2 e-11 [Mpc s^2/km^2]
const double H0 = 1.0e2 *h;
const double H02 = H0 *H0;
const double H0c = H0 /c;
const double rhocr = 2.77537 /h /h;			//	rho_cr e-11 [Msun Mpc^-3]

const double ai = 1.0e-5;
const double Mi = 1.0e3;

double H(double a){ return H0 *sqrt(Om *pow(a,-3) + Ol);}
double H2(double a){ return H02 *(Om *pow(a,-3) + Ol);}
double Hprime(double a){ return H0 *0.5 /sqrt(Om *pow(a,-3) + Ol) *(-3 *Om *pow(a,-4));}
double HprimeOverH(double a){ return -1.5 *Om /(Om *a + Ol *pow(a,4));}
double beta(double a){ return 1.0 + 2.0 *H(a) *rc /c *(1.0 + a *HprimeOverH(a) /3);}
double rs(double a){ return pow(16.0 *GMsc2 *Mi *rc *rc /(9.0 *beta(a) *beta(a)), 1.0/3);}
double g(double x){ double y3=abs(x*x*x); return sqrt(y3 *y3 + y3) - y3;}

double fNFW(double x){ return log(1. + x) - x/(1. + x);}
double Mass(double a, double r, double con){ return Mi *fNFW(con *abs(r) /rvir) /fNFW(con) *pow(a /avir, Gama);}

int MG;

double enhancement_x(double x1, double a){
  if(MG==0) return 1.;
  else return (1.0 + 2.0 *g(Ri *(x1 + a/ai) /rs(a)) /(3.0 *beta(a)));
}

double enhancement_r(double r1, double a){
  if(MG==0) return 1.;
  else return (1.0 + 2.0 *g(r1 /rs(a)) /(3.0 *beta(a)));
}

//	constant mass
double x1prime(double x2){return x2;}
double x2prime(double x1, double x2, double a){
  return -(1.0 +a *HprimeOverH(a)) *x2 /a +(1.0 + a *HprimeOverH(a)) *x1 /a /a - Om *H02 *pow(a,-5) /(2.0 *H2(a))
*enhancement_x(x1,a) *(x1 +a /ai) *((1. + di) *pow(ai *x1/a + 1.0, -3) - 1.);}

//	FRW mass
double r1prime(double r2){return r2;}
double r2prime(double r1, double r2, double a){
  return  -(1.0 +a *HprimeOverH(a)) *r2 /a + H02 *Ol *abs(r1) *pow(a *H(a),-2)
  -r1/abs(r1) *GMs *Mass(a, r1, 3.0) *pow(a *H(a) *r1,-2) *enhancement_r(r1,a);}


///////////////////////////////////////////////////
void iterate(int printing){
///////////////////////////////////////////////////
  double da,dda;
  double k1x1, k2x1, k3x1, k4x1;
  double k1x2, k2x2, k3x2, k4x2;
  double x1,x2;

  double r1,r2;

  FILE *SP;
  SP=fopen("sp_dgp.dat","w+");

  //	integrating with constant mass
  double steep=5.0e-1;
    a=ai;	da=a*steep;	dda=da*0.5;	x1=0.0;	x2=-di/(3.0*ai);
    rta=0.0;
    int initializing=1;
    fprintf(SP,"%le %le\n",a,Ri*(x1 +a /ai));
    while(initializing==1){
      k1x1=x1prime(x2);			k1x2=x2prime(x1,x2,a);				a+=dda;
      k2x1=x1prime(x2+k1x2*dda);	k2x2=x2prime(x1+k1x1*dda,x2+k1x2*dda,a);
      k3x1=x1prime(x2+k2x2*dda); 	k3x2=x2prime(x1+k2x1*dda,x2+k2x2*dda,a);	a+=dda;
      k4x1=x1prime(x2+k3x2*da);      	k4x2=x2prime(x1+k3x1*da,x2+k3x2*da,a);

      x1+=(k1x1+k2x1*2.0+k3x1*2.0+k4x1)/6*da;
      x2+=(k1x2+k2x2*2.0+k3x2*2.0+k4x2)/6*da;
//      da=a*steep;	dda=da*0.5;
      if(printing==1) fprintf(SP,"%le %le\n",a,Ri*(x1 +a /ai));
      if(rta <(Ri*(x1 +a /ai))){
        rta=Ri*(x1 +a /ai);
        ata=a;}
      else{
       if(rta*0.5 <Ri*(x1 +a /ai)){
          r1=Ri*(x1 +a /ai);
          r2=Ri*(x2 +1.0 /ai);
          avir=a;}
        else
          initializing=0;
      }
    }

  //	integrating with NFW * accretion
    rvir=r1;
    int collapsing=1;
    while(collapsing==1){
      k1x1=r1prime(r2);			k1x2=r2prime(r1,r2,a);				a+=dda;
      k2x1=r1prime(r2+k1x2*dda);	k2x2=r2prime(r1+k1x1*dda,r2+k1x2*dda,a);
      k3x1=r1prime(r2+k2x2*dda); 	k3x2=r2prime(r1+k2x1*dda,r2+k2x2*dda,a);	a+=dda;
      k4x1=r1prime(r2+k3x2*da);      	k4x2=r2prime(r1+k3x1*da,r2+k3x2*da,a);

      r1+=(k1x1+k2x1*2.0+k3x1*2.0+k4x1)/6*da;
      r2+=(k1x2+k2x2*2.0+k3x2*2.0+k4x2)/6*da;
//      da=a*steep;	dda=da*0.5;
      if(printing==1) fprintf(SP,"%le %le\n",a,r1);
      if((r1>0.0)&&(r2<0.0))
        acol=a;
      else{
        if(r2>0.0){
          asp=a;
          rsp=r1;
          collapsing=0;
        }
      }
    }
  fclose(SP);
}


///////////////////////////////////////////////////
///////////////////////////////////////////////////
int main(){
///////////////////////////////////////////////////
///////////////////////////////////////////////////

FILE * GA;
char FileName[200];
int rc_ind=-1;
//for(rc_ind=-1;rc_ind<10;rc_ind++){
if(rc_ind==-1){	rc=0.; 					MG=0;}
else{		rc = double(rc_ind+1)*0.03 *c / H0;	MG=1;}
sprintf(FileName,"sp_gama_rc%d.dat",rc_ind);
GA=fopen(FileName,"w+");
fprintf(GA,"# rc=%le\n",rc);
//for(Gama=1.;Gama<5.;Gama+=0.1){
Gama=3.;
  di = 2.81e-5;
  Ri = pow(3.0 *Mi /(4.0 *pi *Om *rhocr *pow(ai,-3) *(1.0 + di)), 1.0 /3);printf("Ri=%le\n",Ri);
  double ddi=di/10;
  double tolerance=1.0e-4;
  int sign=1;
  a=ai;
//  while(abs(a-1.0)>tolerance){
    iterate(1);
    printf("%le ",a);fflush(stdout);
    if((a-1.0)>0.){
      if(sign<0) ddi*=-0.5;
      else;
      di+=ddi;
      sign=1;}
    else{
      if(sign>0) ddi*=-0.5;
      else;
      di+=ddi;
      sign=-1;}
//  }	//a=1
  printf("\n\n");
  fprintf(GA,"%le %le\n",Gama,abs(rsp/rvir));fflush(GA);
//}	//gamma
fclose(GA);
//}	//rc


//  iterate(1);

    printf("\nata=%le\n",ata);
    printf("avir=%le\n",avir);
    printf("rta=%le\n",rta);
    printf("rvir=%le\n",rvir);
    printf("acol=%le\n",acol);
    printf("rsp=%le\n",rsp);
    printf("asp=%le\n",asp);
    printf("rsp/rvir=%le\n",abs(rsp/rvir));
  return 0;

}


