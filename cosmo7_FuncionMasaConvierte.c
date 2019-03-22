//	dn=f/prefactor
//	entonces f=dn*prefactor
//	de esta forma voy a obtener f(sigma) a partir de las simulaciones
//	la idea es encontrar una función que lo ajuste, un f_fit(sigma)


//	entonces: la idea es calcular el prefactor
//	interolar la abundancia medida
//	y entonces calcular f para ajustarla con una función de sigma



//paquetes gnu
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
//NR3
#include "../l/nr3.h"
#include "../l/interp_1d.h"
#include "../l/parametros.h"


  Doub L=1.;//	1./1.717;			//	r=1.717rl    dn/dlnr = (f/v) dlns_/dlnrl |rl(r)
const double Og=0.0000475;			//	fracción de radiacion
	double Mv,Mc,E2,rho_c,o_m,Dc,r3v,rv,cc,rho_s;		//constantes NFW
const double Ok=0.;							//	fracción de curvatura hoy

const double dcc=1.686;
const double rhom=2.77536627*Om*1e11;

using namespace std;
const double pi=4.*atan(1.);		//	numero pi

inline double W_dW(double x, double &w, double &dw){//	funcion ventana al cuadrado
  double xc=x*cos(x);
  double s=sin(x);
  double x2=pow(x,-2);
  double x4=3.*x2*x2;
  w=x*x4*(s-xc);
  dw=x4*(3.*xc+(x*x-3.)*s);}


inline double T(double sigma){					//	función de masa de Tinker
//	return 0.186*(1.+pow(sigma/2.57,-1.47))*exp(-1.19*pow(sigma,-2.));} 
	return 0.482*(pow(sigma,-1.97)+pow(sigma,-0.51))*exp(-1.228*pow(sigma,-2.));}


int main(){
for(int t=0;t<=6;t++){

//////////////////////////////// Leyendo el espectro de potencias interpolado del CAMB
  Int nk=0,i,j;
  Doub aux1,aux2,kmin=1.*pow(10,8),kmax=-1.;

  char NomArch[400];
  FILE * D;
  sprintf(NomArch,"/home/cosmousp/Descargas/CAMB/my_folder/planck_2018_matterpower_0.dat");
  D=fopen(NomArch,"r");//1000x5.dat
  while(fscanf(D,"%le %le\n",&aux1, &aux2)!=EOF){
    nk++;
    if(kmin>aux1)	kmin=aux1;
    if(kmax<aux1)	kmax=aux1;}
  printf("kmin=%le, kmax=%le h/Mpc\n",kmin,kmax);fflush(stdout);
  VecDoub K(nk);
  VecDoub P(nk);
  rewind(D);//1000x5
  for(i=0;i<nk;i++){
    fscanf(D,"%le   %le\n",&aux1, &aux2);
    K[i]=aux1;
    P[i]=aux2;}
  fclose(D);


  cout<<"leyendo los datos de dnh_512x256_03.txt"<<endl;
  FILE * UN;
  sprintf(NomArch,"../%c/dnh_512x256_03.txt",NameTheories[t]);
  UN=fopen(NomArch,"r");
  double mas,ab,dab,rad;//	abundancia
  int Nr=0;
  while(fscanf(UN,"%le %le %le\n",&mas,&ab,&dab)!=EOF)
    Nr++;
  VecDoub masa(Nr);
  VecDoub radio(Nr);
  VecDoub sigma(Nr);
  VecDoub dsigma(Nr);
  VecDoub abundancia(Nr);
  VecDoub dabundancia(Nr);
  Doub prefac;
  rewind(UN);
  for(i=0;i<Nr;i++){
    fscanf(UN,"%le %le %le\n",&masa[i],&abundancia[i],&dabundancia[i]);
    radio[i]=pow(3*masa[i]*1e-11/(4*pi*2.77536627*Om),1./3);
  }
  fclose(UN);


//////////////////////////////// Integrando sigma^2 para algunos valores de R(M)


  FILE * S;
  double w,dw,wm,dwm,sigma2;
  sprintf(NomArch,"../%c/sigma.txt",NameTheories[t]);
  S=fopen(NomArch,"w+");
  for(j=0;j<Nr;j++){//	calcula sigma y derivada para todos los radios en que tenemos abundancia
    sigma2 = 0.0;
    dsigma[j] = 0.0;
    for(i=0;i<nk-1;i++){///	integral
      W_dW(radio[j]*K[i],w,dw);
      W_dW(radio[j]*K[i+1],wm,dwm);
       sigma2   += (K[i+1]-K[i])*(P[i]*pow(K[i],2)*pow(w,2) + P[i+1]*pow(K[i+1],2)*pow(wm,2));
      dsigma[j] += (K[i+1]-K[i])*(P[i]*pow(K[i],3)*w*dw     + P[i+1]*pow(K[i+1],3)*wm*dwm);
    }
    sigma[j]=sqrt(sigma2)/(2.0*pi);
    dsigma[j]/=(4.0*pi*pi);
    fprintf(S,"%le %le\n",radio[j],sigma[j]);
    }
  fclose(S);





  FILE * UF;
  sprintf(NomArch,"../%c/fh_medida.dat",NameTheories[t]);
  UF=fopen(NomArch,"w+");
  fprintf(UF,"# sigma f fmax fmin\n");
  for(i=0;i<Nr;i++){
    prefac=-pow(sigma[i],2)*masa[i]*3/(2.77536627*Om*1e11*radio[i]*dsigma[i]);
    fprintf(UF,"%le %le %le %le %le %le %le %le %le\n",
      masa[i],
      -T(sigma[i])*2.77536627*Om*1e11*radio[i]*dsigma[i]/(masa[i]*sigma[i]*sigma[i]*3.),
      T(sigma[i]),
      abundancia[i],
      dabundancia[i],
      1./sigma[i],                   //	sigma
      abundancia[i]*prefac,				//	f_medida
      (abundancia[i]+dabundancia[i])*prefac,				//	f_max
      (abundancia[i]-dabundancia[i])*prefac);				//	f_min
    }
  fclose(UF);
  
}
cout<<"acaba"<<endl;
	return 0;

}


