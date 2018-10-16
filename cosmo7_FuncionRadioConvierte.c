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


const double dcc=1.686;

using namespace std;
const double pi=4.*atan(1.);		//	numero pi

inline double W_dW(double x, double &w, double &dw){//	funcion ventana al cuadrado
  double xc=x*cos(x);
  double s=sin(x);
  double x2=pow(x,-2);
  double x4=3.*x2*x2;
  w=x*x4*(s-xc);
  dw=x4*(3.*xc+(x*x-3.)*s);}




int main(){
for(int t=0;t<=6;t++){

//////////////////////////////// Leyendo el espectro de potencias interpolado del CAMB
  Int nk=0,i,j;
  Doub aux1,aux2,kmin=1.*pow(10,8),kmax=-1.;

  char NomArch[400];
  FILE * D;
  sprintf(NomArch,"../%c/Matter_Power_MG.dat",NameTheories[t]);
  D=fopen(NomArch,"r");//1000x5.dat
  while(fscanf(D,"%lf %lf\n",&aux1, &aux2)!=EOF){
    nk++;
    if(kmin>aux1)	kmin=aux1;
    if(kmax<aux1)	kmax=aux1;}
  printf("kmin=%le, kmax=%le h/Mpc\n",kmin,kmax);fflush(stdout);
  VecDoub K(nk);
  VecDoub P(nk);
  rewind(D);//1000x5
  for(i=0;i<nk;i++){
    fscanf(D,"%lf   %lf\n",&aux1, &aux2);
    K[i]=aux1;
    P[i]=aux2;}
  fclose(D);


  cout<<"leyendo los datos de dnv_esferas_512x256_03.txt"<<endl;
  FILE * UN;
  sprintf(NomArch,"../%c/dnv_esferas_512x256_03.txt",NameTheories[t]);
  UN=fopen(NomArch,"r");
  float rad;
  double ab,dab;//	abundancia
  int Nr=0;
  while(fscanf(UN,"%f %le %le\n",&rad,&ab,&dab)!=EOF)
    Nr++;
  VecDoub radio(Nr);
  VecDoub sigma(Nr);
  VecDoub dsigma(Nr);
  VecDoub abundancia(Nr);
  VecDoub dabundancia(Nr);
  Doub prefac;
  rewind(UN);
  Doub radi;
  for(i=0;i<Nr;i++){
    fscanf(UN,"%f %le %le\n",&rad,&abundancia[i],&dabundancia[i]);
    radi=double(rad);
    radio[i]=radi;
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
  sprintf(NomArch,"../%c/f_medida.dat",NameTheories[t]);
  UF=fopen(NomArch,"w+");
  fprintf(UF,"# sigma f fmax fmin\n");
  for(i=0;i<Nr;i++){
    prefac=-4.*pi*pow(sigma[i]*radio[i],2)/(3.*dsigma[i]);
    fprintf(UF,"%le %le %le %le\n",
      sigma[i],                   //	sigma
      abundancia[i]*prefac,				//	f_medida
      (abundancia[i]+dabundancia[i])*prefac,				//	f_max
      (abundancia[i]-dabundancia[i])*prefac);				//	f_min
    }
  fclose(UF);
  
}
cout<<"acaba"<<endl;
	return 0;

}


