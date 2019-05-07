//  Based on the toy model for the estimation of r_sp including nDGP gravity: 1806.04302
//  in this case we have a diferential equation for matter without any dependence on the scalar field
//  SplashBack_nDGP.cpp
//  already made by 1805.09918

//  f(R) does not have an analitic expression for phi', then we are using the Fourier formulation used by 1805.09918
//  in this case we have to compute a couple of fourier transforms in the source term
//  SplashBack_fR.cpp

//  if we want to model the symmetron case there no other option but compute the phi profile
//  in this case we need to solve the scalr field equation and profile in order to source the matter equation
//  SplashBack_symmetron.cpp
//  already made for other people under the self-similar approximation (in EdS space)

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

const double pi=4.0*atan(1.0);	//	number pi

double Gama,di,dvir,Ri,a,rc;

// all masses in 10^11 Ms/h units
const double Om	= 3.0e-1;
const double Ol = 1.0e0 - Om;
const double h = 7.0e-1;

//const double Mpc = 3.08567758149 *pow(10, 22);
//const double GMsc2 = 2953.25024 *0.5 /Mpc *pow(10, 11);	//	G Msun c^2 /2 e-11 [Mpc]
const double Mpc = 3.08567758149;			//	Mpc e-22 [m]
const double GMsc2 = 2953.25024e-11 *0.5 /Mpc;		//	G Msun c^2 /2 e-11 [Mpc]
const double c = 2.99792458e5;
const double GMs = GMsc2 *c *c;				//	G Msun /2 e-11 [Mpc s^2/km^2]
const double H0 = 1.0e2 *h;
const double H02 = H0 *H0;
const double H0c = H0 /c;
const double rhocr = 2.77537e0 /h /h;			//	rho_cr e-11 [Msun Mpc^-3]

const double ai = 1.0e-5;
const double Mi = 1.0e3;

//	density profile needed to solve the fR case
const int n_r=200;				//	number of spherical shells to follow
const int n_k=n_r;				//	number of scales to sample ~delta
double d[n_r],dp[n_r],source[n_r],td[n_k];	//	delta(r) delta'(r) nabla2Psi(r) and ~delta(k)
double d1[n_r],d2[n_r];				//	delta and delta' aux for rungekutta steps
double r[n_r],k[n_k];				//	r ,k

double rta[n_r],rvir[n_r],rsp,ata[n_r],avir[n_r],acol,asp,rad,vel;

double r_min=8.5e-3, r_max=3.0e0;
double k_min=2.0*pi*8.5e-3, k_max=2.0*pi*3.0e0;

double H(double a){ return H0 *sqrt(Om *pow(a,-3) + Ol);}
double H2(double a){ return H02 *(Om *pow(a,-3) + Ol);}
double E2a5(double a){ double a2=a*a,a5=a2*a2*a;return (Om *a2 + Ol *a5);}
double E2a2(double a){ return (Om /a + Ol *a*a);}
double Hprime(double a){ return H0 *0.5 /sqrt(Om *pow(a,-3) + Ol) *(-3 *Om *pow(a,-4));}
double HprimeOverH(double a){ return -1.5 *Om /(Om *a + Ol *pow(a,4));}
double beta(double a){ return 1.0 + 2.0 *H(a) *rc /c *(1.0 + a *HprimeOverH(a) /3);}
double rs(double a){ return pow(16.0 *GMsc2 *Mi *rc *rc /(9.0 *beta(a) *beta(a)), 1.0/3);}
double g(double x){ double y3=abs(x*x*x); return sqrt(y3 *y3 + y3) - y3;}

double fR0;
int ind_ri;	//	the index of the edge of the profile at the initial time (r=1)
int ind_rta;	//	index of the latest shell that turned-around
int ind_rvir;	//	index of the more extern virialized shell

const double Rv=1.081945e+00*0.5;//pow(3.0 *Mi /(4.0 *pi *Om *rhocr *200.), 1.0 /3);	//	is r_200 at a=1 given Mi
double fNFW(double x){ return log(1.0 + x) - x/(1.0 + x);}
double DeltaVirial(double a, double r, double con){
//  return Mi *fNFW(con *abs(r/a) /Rv) /(4. *pi *rhocr *Om *pow(r,3) *fNFW(con)) -1.0;}
  return 200.0 *fNFW(con *abs(r) /Rv) /fNFW(con) *pow(Rv/r,3) -1.0;}
double Mass(double a, double r, double con){ return Mi *fNFW(con *abs(r) /Rv) /fNFW(con) *pow(a /avir[ind_ri], Gama);}

//	constant mass nDPG
/*double x1prime(double x2){return x2;}
double x2prime(double x1, double x2, double a){
  return -(1.0 +a *HprimeOverH(a)) *x2 /a +(1.0 + a *HprimeOverH(a)) *x1 /a /a - Om *H02 *pow(a,-5) /(2.0 *H2(a))
*(1.0 + 2.0 *g(Ri *(x1 + a/ai) /rs(a)) /(3.0 *beta(a))) *(x1 +a /ai) *((1. + di) *pow(ai *x1/a + 1.0, -3) - 1.);}

//	FRW mass nDPG
double r1prime(double r2){return r2;}
double r2prime(double r1, double r2, double a){
  return  -(1.0 +a *HprimeOverH(a)) *r2 /a + H02 *Ol *abs(r1) *pow(a *H(a),-2)
  -r1/abs(r1) *GMs *Mass(a, r1, 3.0) *pow(a *H(a) *r1,-2) *(1.0 + 2.0 *g(r1 /rs(a)) /(3.0 *beta(a)));}
*/

///////////////////////////////////////////////////////////////////////
//	f(R)
///////////////////////////////////////////////////////////////////////


void k_r_sample(void){
  double c_r=pow(r_max/r_min,1.0/(n_r-1));
  // you may want to multiply it by the current r in order to cover the entire halo
  for(int ind_r=0;ind_r<n_r;ind_r++){
    r[ind_r]=r_min*pow(c_r,ind_r);
    if(r[ind_r]<0.8)
      ind_ri=ind_r;}
  printf("index of the halos edge at ti =%d\n",ind_ri);
  // again you may want divide it by the curren calue of r in order to cover to not lose any information
  double c_k=pow(k_max/k_min,1.0/(n_k-1));
  for(int i=0;i<n_k;i++)
    k[i]=k_min*pow(c_k,i);
}

const double epsilon1=2.0*fR0/(H0c*H0c*Om);	//	Ri gives the right units to k
const double lm_ratio=4.0*Ol/Om;
const double epsilon2=pow(1.0+lm_ratio,-2);

//	fR correction to the source term
double epsilon(double k){
  double k2=k*k*epsilon1;
  return k2/3./(k2 + a *a *epsilon2 *pow(pow(a,-3) +lm_ratio,3) );
}

//	spherical Fourier kernel
double SFT(double x){
return sin(x)/x;
}

//	spherical Fourier kernel2
double SFT2(double x){
return sin(x)-x*cos(x);
}

//	Fourier transform of the matter profile
void calc_tilde_delta(void){
  int ind_k,ind_r;
  for(ind_k=0;ind_k<n_k;ind_k++){
    td[ind_k]=0.0;
    for(ind_r=0;ind_r<n_r-1;ind_r++)
      td[ind_k] += (r[ind_r+1] - r[ind_r])*(pow(r[ind_r],2)  *d1[ind_r]  *SFT(k[ind_k]*r[ind_r]) 
                                          + pow(r[ind_r+1],2)*d1[ind_r+1]*SFT(k[ind_k]*r[ind_r+1]));}
}

//	inverse Fourier transform of the gravitational potential
void calc_source(void){
  int ind_k,ind_r;
  for(ind_r=ind_rvir+1;ind_r<n_r;ind_r++){//	just collapsing shells
    source[ind_r]=0.0;
    for(ind_k=0;ind_k<n_k-1;ind_k++)
      source[ind_r]+=(k[ind_k+1]-k[ind_k])
      *(pow(k[ind_k],2)  *td[ind_k]  *SFT(k[ind_k]  *r[ind_r])*(1.0+epsilon(k[ind_k]/Ri))
       +pow(k[ind_k+1],2)*td[ind_k+1]*SFT(k[ind_k+1]*r[ind_r])*(1.0+epsilon(k[ind_k+1]/Ri)));
    source[ind_r]*=0.5/pi;
  }
  if(ind_rvir>=ind_ri){	//	orviting from the halo's edge
    source[ind_ri]=0.0;
    for(ind_k=0;ind_k<n_k-1;ind_k++)
      source[ind_r]+=(k[ind_k+1]-k[ind_k])
      *(pow(k[ind_k],-1)  *td[ind_k]  *SFT2(k[ind_k]  *abs(rad))*(1.0+epsilon(k[ind_k]/Ri))
       +pow(k[ind_k+1],-1)*td[ind_k+1]*SFT2(k[ind_k+1]*abs(rad))*(1.0+epsilon(k[ind_k+1]/Ri)));
      source[ind_r]*=0.5/pi;}
}


//	collapsing f(R)
double d1prime(double d2){return d2;}
double d2prime(double d1, double d2, double a, double source){
//  return -(3.0/a +HprimeOverH(a)) *d2 +4.0 *d2 *d2 /(3.0 *(1.0 +d1)) +3.0 *Om *(1.0 +d1) /(2.0 *E2a5(a)) *d1;}
  return -(3.0/a +HprimeOverH(a)) *d2 +4.0 *d2 *d2 /(3.0 *(1.0 +d1)) +3.0 *Om *(1.0 +d1) /(2.0 *E2a5(a)) *source;}


//	orviting f(R)
double r1prime(double r2){return r2;}
double r2prime(double r1, double r2, double a, double source){
  return  -(1.0 +a *HprimeOverH(a)) *r2 /a + H02 *Ol *abs(r1) *pow(a *H(a),-2)
  -r1/abs(r1)  *pow(r1,-2) *3.0 *Om /(2.0 *E2a2(a)) *source;}
//  -r1/abs(r1) *GMs *Mass(a, r1, 3.0) *pow(a *H(a) *r1,-2);}

//	charge the initial density profile and set all shperical shells to collapse
void set_initial_profile(int type){
  double s=0.3;
  for(int ind_r=0;ind_r<n_r;ind_r++){
    if(type==0){	//	top-hat
      if(r[ind_r]<=Rv)	d[ind_r]=di;	else	d[ind_r]=0.;}
    if(type==1)		//	tanh
      d[ind_r]=di*0.5*(1.0-tanh( (r[ind_r]-1.) /s ));
    dp[ind_r]=d[ind_r]/ai;				//	velocity
    rta[ind_r]=-1.0;
    rvir[ind_r]=-1.0;
  }
  ind_rta=-1;						//	no shell has already turned-around
  ind_rvir=-1;						//	no shell has already virialized
}


// returns the radius of the shell with index ind_r, computed by using the density inside that shell
double get_rad(int ind_r){
  return a *pow(3.0 *Mi /(4.0 *pi *Om *rhocr *(1.0 + d[ind_r])), 1.0 /3);
}

//	returns dr/da for the given shell index
double get_vel(int ind_r){
  return get_rad(ind_r) *(1.0 /a -dp[ind_r] /(3.0 *(1.0 +d[ind_r])));
}


///////////////////////////////////////////////////
void iterate(int printing){
///////////////////////////////////////////////////
  double da,dda;
  double k1d1[n_r], k2d1[n_r], k3d1[n_r], k4d1[n_r];
  double k1d2[n_r], k2d2[n_r], k3d2[n_r], k4d2[n_r];

  double r1,r2;
  double k1x1, k2x1, k3x1, k4x1;
  double k1x2, k2x2, k3x2, k4x2;

  FILE *SP;
  SP=fopen("sp.dat","w+");

  //	integrating with constant mass
  double steep=5.0e-0;
  a=ai;	da=a*steep;	dda=da*0.5;
  set_initial_profile(1);	//	initialize d1 and d2
  int running=1;

  //	printing profile
  /*FILE *INI;
  INI=fopen("data.dat","w+");
  for(int id_r=0;id_r<n_r;id_r++) d1[id_r]=d[id_r];
    calc_tilde_delta();			calc_source();
  for(int id_r=0;id_r<n_r;id_r++) fprintf(INI,"%le %le %le\n",r[id_r],d[id_r],source[id_r]);
  fclose(INI);*/

  for(int id_r=0;id_r<n_r;id_r++){
    d1[id_r]=d[id_r];			d2[id_r]=dp[id_r];}

  while(running==1){
    //	evolve only no virialized shells
    //	when the edge reaches virial state we start the orbit path
    calc_tilde_delta();			calc_source();
    for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
      k1d1[id_r]=d1prime(d2[id_r]);	k1d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
      d1[id_r]=d[id_r]+k1d1[id_r]*dda;	d2[id_r]=dp[id_r]+k1d2[id_r]*dda;}
    if(ind_rvir>=ind_ri){
      k1x1=r1prime(r2);			k1x2=r2prime(r1,r2,a,source[ind_ri]);
      rad=r1+k1x1*dda;			vel=r2+k1x2*dda;}

    a+=dda;
    calc_tilde_delta();			calc_source();
    for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
      k2d1[id_r]=d1prime(d2[id_r]);	k2d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
      d1[id_r]=d[id_r]+k2d1[id_r]*dda;	d2[id_r]=dp[id_r]+k2d2[id_r]*dda;}
    if(ind_rvir>=ind_ri){
      k2x1=r1prime(vel);		k2x2=r2prime(rad,vel,a,source[ind_ri]);
      rad=r1+k2x1*dda;			vel=r2+k2x2*dda;}

    calc_tilde_delta();			calc_source();
    for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
      k3d1[id_r]=d1prime(d2[id_r]); 	k3d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
      d1[id_r]=d[id_r]+k3d1[id_r]*da;	d2[id_r]=dp[id_r]+k3d2[id_r]*da;}
    if(ind_rvir>=ind_ri){
      k3x1=r1prime(vel);	 	k3x2=r2prime(rad,vel,a,source[ind_ri]);
      rad=r1+k3x1*da;			vel=r2+k3x2*da;}

    a+=dda;
    calc_tilde_delta();			calc_source();
    for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
      k4d1[id_r]=d1prime(d2[id_r]);     k4d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);

      d[id_r]+= (k1d1[id_r] +k2d1[id_r]*2.0 +k3d1[id_r]*2.0 +k4d1[id_r])/6*da;		d1[id_r]=d[id_r];
      dp[id_r]+=(k1d2[id_r] +k2d2[id_r]*2.0 +k3d2[id_r]*2.0 +k4d2[id_r])/6*da;		d2[id_r]=dp[id_r];}

    if(ind_rvir>=ind_ri){
      k4x1=r1prime(vel);      		k4x2=r2prime(rad,vel,a,source[ind_ri]);
      r1+=(k1x1+k2x1*2.0+k3x1*2.0+k4x1)/6*da;						rad=r1;
      r2+=(k1x2+k2x2*2.0+k3x2*2.0+k4x2)/6*da;						vel=r2;}

    if(printing==1) fprintf(SP,"%le %le\n",a,get_rad(ind_ri));
    fprintf(SP,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",a
      ,get_rad(ceil(double(ind_ri)/5)),get_rad(ceil(double(ind_ri)/2)),get_rad(ceil(double(ind_ri)/1.4)),get_rad(ceil(double(ind_ri)/1.2)),get_rad(ind_ri),get_rad(ceil(double(ind_ri)*1.2)),r1
      ,d[int(ceil(double(ind_ri)/5))] ,d[int(ceil(double(ind_ri)/2))] ,d[int(ceil(double(ind_ri)/1.4))] ,d[int(ceil(double(ind_ri)/1.2))] ,d[ind_ri]      ,d[int(ceil(double(ind_ri)*1.2))]);fflush(SP);


    //	stopping evolution when edge orbit reach second turn-around
/*    if(ind_rvir>=ind_ri){
      if((r1>0.0)&&(r2<0.0))
        acol=a;
      else{
        if(r2>0.0){
          asp=a;
          rsp=r1;
          running=0;
        }
      }
    }
*/

    //	labeling shells that already viarilized, also saving its avir and rvir
    /*if((ind_rvir<ind_rta)&&(rta[ind_rvir]*0.5 >get_rad(ind_rvir))){
      rvir[ind_rvir]=get_rad(ind_rvir);
      avir[ind_rvir]=a;
      printf("%d vir at %le when %le\n",ind_rvir,rvir[ind_rvir],avir[ind_rvir]);fflush(stdout);
      if(ind_rvir==ind_ri){	//	if edge just virilized then initialize orbit
        printf("%d is the halo's edge\n",ind_rvir);fflush(stdout);
        r1=get_rad(ind_rvir);
        r2=get_vel(ind_rvir);}
      ind_rvir++;}*/

  for(int id_r=0;id_r<n_r;id_r++){
    //	in last scenario we freeze the profile when it shell gets close to a virial top-hat
    //	instead, we want to freeze the profile for it to be a NFW: (completally indepndet of turn-around)
    if((rvir[id_r]<0.0)&&(rta[id_r]>0.0)&&( d[id_r] >= DeltaVirial(a,get_rad(id_r),3.0) )){// 3=conc, may we vary it with z?
      rvir[id_r]=get_rad(id_r);
      avir[id_r]=a;
      if(id_r==ind_ri){	//	if edge just virilized then initialize orbit
        r1=get_rad(id_r);	rad=r1;
        r2=get_vel(id_r);	vel=r2;
        dvir=d[id_r];
        printf("edge reaching virial condition, r1=%le r2=%le, d1=%le D=%le\n",r1,r2,d[id_r],DeltaVirial(a,get_rad(id_r),3.0));}
      d[id_r] = DeltaVirial(a,get_rad(id_r),3.0);
      ind_rvir=id_r;
      printf("%d vir at %le when %le reaching %le\n",ind_rvir,rvir[ind_rvir],avir[ind_rvir],d[id_r]);fflush(stdout);}

    //	labeling shells that have already turned-around, also saving its ata and rta
    if((rta[id_r]<0.0)&&(dp[id_r]>3.0 *(1.0 +d[id_r]) /a)){
      rta[id_r]=get_rad(id_r);
      ata[id_r]=a;
      ind_rta=id_r;
      if(id_r==ind_ri)printf("edge reaching turn-around condition, r1=%le a=%le\n", rta[id_r],a);
      printf("%d t-a at %le when %le\n",ind_rta,rta[ind_rta],ata[ind_rta]);fflush(stdout);}
  }
if(a>2.)running=0;
}
  fclose(SP);
}

/*
///////////////////////////////////////////////////
void integrate_nDGP(void){
///////////////////////////////////////////////////
FILE * GA;
char FileName[200];
int rc_ind;
for(rc_ind=0;rc_ind<10;rc_ind++){
rc = double(rc_ind+1)*0.03 *c / H0;
sprintf(FileName,"sp_gama_rc%d.dat",rc_ind);
GA=fopen(FileName,"w+");
fprintf(GA,"# rc=%le\n",rc);
for(Gama=1.;Gama<5.;Gama+=0.1){
  di = 1.8e-5;
  double ddi=di/10;
  double tolerance=1.0e-4;
  int sign=1;
  a=ai;
  while(abs(a-1.0)>tolerance){
    iterate(0);
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
  }
  printf("\n\n");
  fprintf(GA,"%le %le\n",Gama,abs(rsp/rvir));fflush(GA);
}
fclose(GA);
}
//  iterate(1);
}

*/




///////////////////////////////////////////////////
void integrate_fR(void){
///////////////////////////////////////////////////
  k_r_sample();

FILE * GA;
char FileName[200];
int fR0_ind=0;
double fR0min=1.0e-6;
double fR0max=1.0e-4;
double c_fR0=pow(fR0max/fR0min,1./(10-1));
//for(fR0_ind=0;fR0_ind<10;fR0_ind++){
//fR0 = 1.0e-6 *pow(c,fR0_ind);
fR0=0.;
sprintf(FileName,"sp_fR%d.dat",fR0_ind);
GA=fopen(FileName,"w+");
fprintf(GA,"# fR0=%le\n",fR0);
//for(Gama=1.;Gama<5.;Gama+=0.1){
  di = 3.7e-5;printf("di=%le\n",di);
  Ri = ai *pow(3.0 *Mi /(4.0 *pi *Om *rhocr *(1.0 + di)), 1.0 /3);
  double ddi=di/10;
  double tolerance=1.0e-4;
  int sign=1;
  a=ai;

//	just print original profile and FT^-1(FT profile)
/*    FILE *FT;
    char file[200];
    charge_initial_profile(0);
    sprintf(file,"FT%d.dat",0);
    FT=fopen(file,"w+");
    for(int id_r=0;id_r<n_r;id_r++)	d1[id_r]=d[id_r];
    calc_tilde_delta();			calc_source();
    for(int id_r=0;id_r<n_r;id_r++)	fprintf(FT,"%le %le %le\n",r[id_r],d[id_r],source[id_r]);
    fclose(FT);

    for(int index=1;index<10;index++){
      sprintf(file,"FT%d.dat",index);
      FT=fopen(file,"w+");
      for(int mas=0;mas<100;mas++){
        for(int id_r=0;id_r<n_r;id_r++)	d1[id_r]=source[id_r];
        calc_tilde_delta();			calc_source();}
      for(int id_r=0;id_r<n_r;id_r++)	fprintf(FT,"%le %le %le\n",r[id_r],d[id_r],source[id_r]);
      fclose(FT);
    }
*/


printf("ri=%le rv=%le\n",Ri,Rv);fflush(stdout);
printf("Delta=%le\n",DeltaVirial(1., Rv,3.));
printf("Delta=%le\n",DeltaVirial(1., 0.5*Rv,3.));
printf("Delta=%le\n",DeltaVirial(1., 0.25*Rv,3.));
printf("Delta=%le\n",DeltaVirial(1., 0.125*Rv,3.));
printf("Delta=%le\n",DeltaVirial(1., 0.0625*Rv,3.));

//	evolve as usual
    iterate(0);



/*  while(abs(a-1.0)>tolerance){
    iterate(0);
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
  }*/
  printf("\n\n");
  fprintf(GA,"%le %le\n",Gama,abs(rsp/rvir[ind_ri]));fflush(GA);
//}	//	gama
fclose(GA);
//}	//	fR0



}





///////////////////////////////////////////////////
///////////////////////////////////////////////////
int main(){
///////////////////////////////////////////////////
///////////////////////////////////////////////////

//  integrate_nDGP();
  integrate_fR();

for(int i=ind_ri;i<ind_ri+1;i++){
    printf("\nata=%le\n",ata[i]);
    printf("rta=%le\n\n",rta[i]);

    printf("avir=%le\n",avir[i]);
    printf("rvir=%le\n\n",rvir[i]);
}
    printf("acol=%le\n\n",acol);

    printf("asp=%le\n",asp);
    printf("rsp=%le\n",rsp);
    printf("rsp/rvir=%le\n",abs(rsp/rvir[ind_ri]));
  return 0;

}


