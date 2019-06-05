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


// last version used an static sampling for r and k while computing the direct and inverse FT
// now we update the sampling for every timestep

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

const double pi=4.0*atan(1.0);	//	number pi

double Gama,Di,dvir,Ri,a,rc,R0,norm_epsilon,norm_k2;	//R0 scaling units, Ri more physical edge

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
const double c_H0 = c /H0;
const double rhocr = 2.77537e0 /h /h;			//	rho_cr e-11 [Msun Mpc^-3]

const double ai = 1.0e-5;
const double Mi = 1.0e3;

//	density profile needed to solve the fR case
const int n_r=500;				//	number of spherical shells to follow
const int n_k=n_r;				//	number of scales to sample ~delta
double d[n_r],dp[n_r],source[n_r],td[n_k];	//	delta(r) delta'(r) nabla2Psi(r) and ~delta(k)
double d1[n_r],d2[n_r];				//	delta and delta' aux for rungekutta steps
double r[n_r],k[n_k];				//	r ,k
double ri[n_r],ki[n_k],di[n_r];			//	initial r, k and delta

double rta[n_r],rvir[n_r],rsp,ata[n_r],avir[n_r],acol,asp,rad,vel;

double r_min=1.0e-3, r_max=3.0e0;
double k_min=2.0*pi*1.0e-3, k_max=2.0*pi*3.0e0;

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
int ind_r0;	//	the index of the edge of the profile at the initial time (r=1.0), scaling units
int ind_ri;	//	the index of the edge of the profile at the initial time (r=0.8), physical edge
int ind_rta;	//	index of the latest shell that turned-around
int ind_rvir;	//	index of the more extern virialized shell

double Rv;//=1.081945e+00*0.5;//pow(3.0 *Mi /(4.0 *pi *Om *rhocr *200.), 1.0 /3);	//	is r_200 at a=1 given Mi
double fNFW(double x){ return log(1.0 + x) - x/(1.0 + x);}
double NFW(double r, double con){
  double rho_s = Mi *pow(con /Rv,3) /( 4.0 *pi *(log(1.0+ con)-con/(1.0 +con))) /(rhocr *Om);//	/(rhoc Om)=delta
  double x= con*r/Rv;
  double xp1= x +1.0;
  return rho_s/(x*xp1*xp1);}

double DeltaVirial(double a, double r, double con){
//  return Mi *fNFW(con *abs(r/a) /Rv) /(4. *pi *rhocr *Om *pow(r,3) *fNFW(con)) -1.0;}
  return 200.0 *fNFW(con *abs(r) /Rv) /fNFW(con) *pow(Rv/r,3) -1.0;}

double Mass(double a, double r, double con){ return Mi *fNFW(con *abs(r) /Rv) /fNFW(con) *pow(a /avir[ind_ri], Gama);}


void k_r_sample(void){
  double c_r=pow(r_max/r_min,1.0/(n_r-1));
  // you may want to multiply it by the current r in order to cover the entire halo
  for(int ind_r=0;ind_r<n_r;ind_r++){
    ri[ind_r]=r_min*pow(c_r,ind_r);
    if(ri[ind_r]<=0.5) ind_ri=ind_r;
    if(ri[ind_r]<=1.0) ind_r0=ind_r;}
  printf("index of the halos edge at ti =%d\n",ind_ri);
  // again you may want divide it by the curren calue of r in order to cover to not lose any information
  double c_k=pow(k_max/k_min,1.0/(n_k-1));
  for(int i=0;i<n_k;i++)
    ki[i]=k_min*pow(c_k,i);
}

// returns the radius of the shell with index ind_r, computed by using the density inside that shell
double get_rad(int ind_r){
//  return a *pow(3.0 *Mi /(4.0 *pi *Om *rhocr *(1.0 + d[ind_r])), 1.0 /3);
  if(ind_r>ind_rvir)
    return R0 *ri[ind_r0] *a /ai *pow((1.0 +di[ind_r])/(1.0 +d[ind_r]) ,1.0 /3);
  else
    return R0 *ri[ind_r0] *avir[ind_r] /ai *pow((1.0 +di[ind_r])/(1.0 +d[ind_r]) ,1.0 /3);
}

//	returns dr/da for the given shell index
double get_vel(int ind_r){
  return get_rad(ind_r) *(1.0 /a -dp[ind_r] /(3.0 *(1.0 +d[ind_r])));
}

//	fR correction to the source term
inline double epsilon(double k,double b){
  double k2=k*k*fR0*norm_k2*norm_epsilon;	//ne depends only in cosmo, nk2 evolves
  double inside= Om *pow(b,-2) +4.0 *Ol *b;
  return k2/3.0/(k2 + b *pow(inside,3));
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
    #pragma omp parallel
    {
    double td_priv=0.0;
    double k_priv=k[ind_k];
    #pragma omp for
    for(ind_r=0;ind_r<n_r-1;ind_r++)
      td_priv += (r[ind_r+1] - r[ind_r])*(pow(r[ind_r],2)  *d1[ind_r]  *SFT(k_priv*r[ind_r]) 
                                        + pow(r[ind_r+1],2)*d1[ind_r+1]*SFT(k_priv*r[ind_r+1]));
    #pragma omp critical
    {
    td[ind_k]+=td_priv;
    }
    }
    }
}

//	inverse Fourier transform of the gravitational potential
void calc_source(double a){
  int ind_k,ind_r;
  norm_k2 =pow(1.0 + d1[ind_r0] ,2.0 /3);
  if(ind_rvir>=ind_ri){	//	orviting from the halo's edge
    source[ind_ri]=0.0;
    #pragma omp parallel
    {
    double source_priv=0.0;
    double r_priv=abs(rad);
    double a_priv=a;
    #pragma omp for
    for(ind_k=0;ind_k<n_k-1;ind_k++)
      source_priv+=(k[ind_k+1]-k[ind_k])
        *(pow(k[ind_k],-1)  *td[ind_k]  *SFT2(k[ind_k]  *r_priv)*(1.0+epsilon(k[ind_k],a_priv))
         +pow(k[ind_k+1],-1)*td[ind_k+1]*SFT2(k[ind_k+1]*r_priv)*(1.0+epsilon(k[ind_k+1],a_priv)));
    #pragma omp critical
    {
    source[ind_ri]+=source_priv;
    }
    }
    source[ind_ri]*=0.5/pi;
  }
  else{
    for(ind_r=ind_rvir+1;ind_r<n_r;ind_r++){//	just collapsing shells
      source[ind_r]=0.0;
      #pragma omp parallel
      {
      double source_priv=0.0;
      double r_priv=r[ind_r];
      double a_priv=a;
      #pragma omp for
      for(ind_k=0;ind_k<n_k-1;ind_k++)
        source_priv+=(k[ind_k+1]-k[ind_k])
                     *(pow(k[ind_k],2)  *td[ind_k]  *SFT(k[ind_k]  *r_priv)*(1.0+epsilon(k[ind_k],a_priv))
                      +pow(k[ind_k+1],2)*td[ind_k+1]*SFT(k[ind_k+1]*r_priv)*(1.0+epsilon(k[ind_k+1],a_priv)));
      #pragma omp critical
      {
      source[ind_r]+=source_priv;
      }
      }
      source[ind_r]*=0.5/pi;
    }
  }
}


//	collapsing f(R)
double d1prime(double d2){return d2;}
double d2prime(double d1, double d2, double a, double source){
//  return -(3.0/a +HprimeOverH(a)) *d2 +4.0 *d2 *d2 /(3.0 *(1.0 +d1)) +3.0 *Om *(1.0 +d1) /(2.0 *E2a5(a)) *source;}
  return -(3.0/a +HprimeOverH(a)) *d2 +4.0 *d2 *d2 /(3.0 *(1.0 +d1)) +3.0 *Om *(1.0 +d1) /(2.0 *E2a5(a)) *d1;}


//	orviting f(R)
double r1prime(double r2){return r2;}
double r2prime(double r1, double r2, double a, double source){
  return  -(1.0 +a *HprimeOverH(a)) *r2 /a + H02 *Ol *abs(r1) *pow(a *H(a),-2)
//  -r1/abs(r1)  *pow(r1,-2) *3.0 *Om /(2.0 *E2a2(a)) *source *pow(a /avir[ind_ri], Gama);}
  -r1/abs(r1) *GMs *Mass(a, r1, 3.0) *pow(a *H(a) *r1,-2);}

//	charge the initial density profile and set all shperical shells to collapse
void set_initial_profile(int type){
  double s=0.3;
  for(int ind_r=0;ind_r<n_r;ind_r++){
    if(type==0){	//	top-hat
      if(r[ind_r]<=Rv)	di[ind_r]=Di;	else	di[ind_r]=0.;}
    if(type==1)		//	tanh
      di[ind_r]=Di*0.5*(1.0-tanh( (ri[ind_r]-1.) /s ));
  }

  R0 = ai *pow(3.0 *Mi /(4.0 *pi *Om *rhocr *(1.0 + Di   )), 1.0 /3);
  Ri = ai *pow(3.0 *Mi /(4.0 *pi *Om *rhocr *(1.0 + di[ind_ri])), 1.0 /3);printf("R0=%le Ri=%le\n",R0,Ri);
  norm_epsilon =2.0 *pow( c_H0 *(Om +4.0 *Ol) *ai /(R0 *pow( 1.0 +di[ind_r0] ,1.0 /3)) ,2);

  //printf("R0=%le Ri=%le d(Ri)=%le Di=%le\n",R0,Ri,di[ind_ri],Di);fflush(stdout);
  #pragma omp parallel for
  for(int ind_r=0;ind_r<n_r;ind_r++){
    d[ind_r]=di[ind_r];
    r[ind_r]=ri[ind_r];
    k[ind_r]=ki[ind_r];
    dp[ind_r]=di[ind_r]/ai;				//	velocity
    rta[ind_r]=-1.0;
    rvir[ind_r]=-1.0;
  }
  ind_rta=-1;						//	no shell has already turned-around
  ind_rvir=-1;						//	no shell has already virialized
}


void freeze_r(int ind_r){	//	when a shell reaches virial state we freeze its radius
  double norm_r =pow( (1.0 +d1[ind_r0]) /(1.0 +di[ind_r0]) ,1.0 /3) / ri[ind_r0];
  r[ind_r]= norm_r *ri[ind_r] *pow((1.0 +di[ind_r])/(1.0 +d1[ind_r]) ,1.0 /3);	// such that r[ind_r0]=1
}

//	imposing the NFW profile
void set_NFW(void){
  #pragma omp parallel for
  for(int id_r=0;id_r<n_r;id_r++){
    r[id_r]=ri[id_r];
    d1[id_r]=NFW(r[id_r]*Rv,3.0);		//	3.0 is the concentration
  }
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
  set_initial_profile(0);	//	initialize d1 and d2
  int running=1;

    if(printing==1) {fprintf(SP,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",a
      ,get_rad(ceil(double(ind_ri)/1.2)),get_rad(ceil(double(ind_ri)/1.15)),get_rad(ceil(double(ind_ri)/1.1)),get_rad(ceil(double(ind_ri)/1.05)),get_rad(ind_ri),get_rad(ceil(double(ind_ri)*1.05)),r1
      ,d[int(ceil(double(ind_ri)/1.2))] ,d[int(ceil(double(ind_ri)/1.15))] ,d[int(ceil(double(ind_ri)/1.1))] ,d[int(ceil(double(ind_ri)/1.05))] ,d[ind_ri]      ,d[int(ceil(double(ind_ri)*1.05))]);fflush(SP);}

  #pragma omp parallel for
  for(int id_r=0;id_r<n_r;id_r++){
    d1[id_r]=d[id_r];			d2[id_r]=dp[id_r];}

  while(running==1){
    //	evolve only no virialized shells
    //	when the edge reaches virial state we start the orbit path

    if(ind_rvir>=ind_ri){
//      calc_tilde_delta();		calc_source(a);
      k1x1=r1prime(r2);			k1x2=r2prime(r1,r2,a,source[ind_ri]);
      rad=r1+k1x1*dda;			vel=r2+k1x2*dda;}
    else{
      #pragma omp parallel for
      for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
        k1d1[id_r]=d1prime(d2[id_r]);	k1d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
        d1[id_r]=d[id_r]+k1d1[id_r]*dda;d2[id_r]=dp[id_r]+k1d2[id_r]*dda;}}

    a+=dda;
    if(ind_rvir>=ind_ri){
//      calc_tilde_delta();		calc_source(a);
      k2x1=r1prime(vel);		k2x2=r2prime(rad,vel,a,source[ind_ri]);
      rad=r1+k2x1*dda;			vel=r2+k2x2*dda;}
    else{
      #pragma omp parallel for
      for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
        k2d1[id_r]=d1prime(d2[id_r]);	k2d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
        d1[id_r]=d[id_r]+k2d1[id_r]*dda;d2[id_r]=dp[id_r]+k2d2[id_r]*dda;}}

    if(ind_rvir>=ind_ri){
//      calc_tilde_delta();		calc_source(a);
      k3x1=r1prime(vel);	 	k3x2=r2prime(rad,vel,a,source[ind_ri]);
      rad=r1+k3x1*da;			vel=r2+k3x2*da;}
    else{
      #pragma omp parallel for
      for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
        k3d1[id_r]=d1prime(d2[id_r]); 	k3d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
        d1[id_r]=d[id_r]+k3d1[id_r]*da;	d2[id_r]=dp[id_r]+k3d2[id_r]*da;}}

    a+=dda;
    if(ind_rvir>=ind_ri){
//      calc_tilde_delta();		calc_source(a);
      k4x1=r1prime(vel);      		k4x2=r2prime(rad,vel,a,source[ind_ri]);
      r1+=(k1x1+k2x1*2.0+k3x1*2.0+k4x1)/6*da;						rad=r1;
      r2+=(k1x2+k2x2*2.0+k3x2*2.0+k4x2)/6*da;						vel=r2;}
    else{
      #pragma omp parallel for
      for(int id_r=ind_rvir+1;id_r<n_r;id_r++){
        k4d1[id_r]=d1prime(d2[id_r]);     k4d2[id_r]=d2prime(d1[id_r],d2[id_r],a,source[id_r]);
        d[id_r]+= (k1d1[id_r] +k2d1[id_r]*2.0 +k3d1[id_r]*2.0 +k4d1[id_r])/6*da;		d1[id_r]=d[id_r];
        dp[id_r]+=(k1d2[id_r] +k2d2[id_r]*2.0 +k3d2[id_r]*2.0 +k4d2[id_r])/6*da;		d2[id_r]=dp[id_r];}}

//    if(printing==1) fprintf(SP,"%le %le\n",a,get_rad(ind_ri));
    if(printing==1){ fprintf(SP,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",a
      ,get_rad(ceil(double(ind_ri)/1.2)),get_rad(ceil(double(ind_ri)/1.15)),get_rad(ceil(double(ind_ri)/1.1)),get_rad(ceil(double(ind_ri)/1.05)),get_rad(ind_ri),get_rad(ceil(double(ind_ri)*1.05)),r1
      ,d[int(ceil(double(ind_ri)/1.2))] ,d[int(ceil(double(ind_ri)/1.15))] ,d[int(ceil(double(ind_ri)/1.1))] ,d[int(ceil(double(ind_ri)/1.05))] ,d[ind_ri]      ,d[int(ceil(double(ind_ri)*1.05))]);fflush(SP);}

    //	following the edge orbit starting from virial state until it reaches second turn-around
    if(ind_rvir>=ind_ri){//printf("ind_rvir=%d\n",ind_rvir);
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
    else{
      for(int id_r=0;id_r<n_r;id_r++){
        //	in last scenario we freeze the profile when it shell gets close to a virial top-hat
        //	instead, we want to freeze the profile for it to be a NFW: (completally indepndet of turn-around)
        if((ind_rvir<n_r)&&(rvir[id_r]<0.0)&&(rta[id_r]>0.0)&&( rta[id_r]*0.5 >get_rad(id_r) )){// 3=conc, may we vary it with z?
          rvir[id_r]=get_rad(id_r);
          avir[id_r]=a;
          ind_rvir=id_r;
          freeze_r(id_r);
          if(id_r==ind_ri){	//	if edge just virilized then initialize orbit
            for(int id=id_r;id<n_r;id++)	avir[id]=a;
            r1=get_rad(id_r);	rad=r1;
            r2=get_vel(id_r);	vel=r2;
            dvir=d[id_r];	ind_rvir=n_r;
            Rv=get_rad(id_r);	//	half of the turn around radius for the selected shell tat represents the edge of the halo
            set_NFW();		//	imposing the NFW profile and stopping the collapsing for every shell
          //printf("edge reaching virial condition, r1=%le r2=%le, d1=%le D=%le\n",r1,r2,d[id_r],DeltaVirial(a,get_rad(id_r),3.0));
          }
          //  d[id_r] = DeltaVirial(a,get_rad(id_r),3.0);

          //  printf("%d vir at %le when %le reaching %le\n",ind_rvir,rvir[ind_rvir],avir[ind_rvir],d[id_r]);fflush(stdout);
        }

        //	labeling shells that have already turned-around, also saving its ata and rta
        if((ind_rta<n_r)&&(rta[id_r]<0.0)&&(dp[id_r]>3.0 *(1.0 +d[id_r]) /a)){
          rta[id_r]=get_rad(id_r);
          ata[id_r]=a;
          ind_rta=id_r;
          //if(id_r==ind_ri)printf("edge reaching turn-around condition, r1=%le a=%le\n", rta[id_r],a);
          //printf("%d t-a at %le when %le\n",ind_rta,rta[ind_rta],ata[ind_rta]);fflush(stdout);
        }
      }
    }
    if(a>2.)running=0;
  }
  fclose(SP);
}



///////////////////////////////////////////////////
void integrate_fR(void){
///////////////////////////////////////////////////

FILE * GA;
char FileName[200];
int fR0_ind=-1;
double fR0min=1.0e-6;
double fR0max=1.0e-4;
double c_fR0=pow(fR0max/fR0min,1./(10-1));

  Di = 2.81e-5;printf("Di=%le\n",Di);
  k_r_sample();
  double tolerance=1.0e-4;
  double ddi;
  int Gama_ind=0;
//for(Gama_ind=0;Gama_ind<11;Gama_ind++){
//Gama=1.0+Gama_ind*2.0*0.1;
Gama=3.0;
sprintf(FileName,"sp_fR_%d.dat",Gama_ind);
GA=fopen(FileName,"w+");

//for(fR0_ind=-1;fR0_ind<10;fR0_ind++){
  if(fR0_ind==-1)	fR0=0.0;
  else		fR0 = 1.0e-6 *pow(c_fR0,fR0_ind);
  //fR0=0.0;

  if(fR0_ind==-1)	ddi=Di*0.05;
  else			ddi=Di*0.02;
  a=ai;
  //	iterate untill achieving sp at a=1
  int bigger=0, smaller=0, last=0;
//  while(abs(a-1.0)>tolerance){
    iterate(1);
    printf("fR0=%le Di=%le ata=%le rta=%le avir=%le rvir=%le acol=%le asp=%le rsp=%le rsp/rvir=%le\n",fR0,Di,ata[ind_ri],rta[ind_ri],
       avir[ind_ri],rvir[ind_ri],acol,asp,rsp,rsp/rvir[ind_ri]);fflush(stdout);
/*    if(a>1.0){
      if(last==-1)	bigger=1;
      last=1;}
    else{
      if(last==1)	smaller=1;
      last=-1;}
    if((bigger>0)&&(smaller>0))	ddi*=0.5;	//	when it already changed direction
    Di+=last*ddi;*/

    if(a>1.0){	bigger=1;	last=1;}
    else{	smaller=1;	last=-1;}
    if((bigger>0)&&(smaller>0))	ddi*=0.5;	//	when it already changed direction
    Di+=last*ddi;
//  }	//	while tolerance
  printf("\n");
  Di-=last*ddi;
//  iterate(1);
  fprintf(GA,"fR0=%le Di=%le ata=%le rta=%le avir=%le rvir=%le acol=%le asp=%le rsp=%le rsp/rvir=%le\n",fR0,Di,
    ata[ind_ri],rta[ind_ri],avir[ind_ri],rvir[ind_ri],acol,asp,rsp,rsp/rvir[ind_ri]);
//}	//	fR0
fclose(GA);
//}	//	gama



}





///////////////////////////////////////////////////
///////////////////////////////////////////////////
int main(){
///////////////////////////////////////////////////
///////////////////////////////////////////////////

  omp_set_num_threads(5);

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





















