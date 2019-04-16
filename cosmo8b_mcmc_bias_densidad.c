#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "cosmo8b_parametros.h"
#include "../l/cosmo8_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR
#include "../l/cosmo8_densidad_bias.h"	//	subrutinas densidad

  using namespace std;

  const int segunda_vez=0;

  const int solo_analiza=1;

  const int NP=100000;					//	numero de puntos de la cadena
  const int BI=1000;					//	burn-in




  int n_salta;				//	pasos que desecha antes de colocar un nuevo eslabón


  double param[Np],sigma[Np];


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void carga(int t,int s){
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
  param[0]=0.55-0.014*pow(radio_stack[t][s],2);	//	b
  param[1]=5.;					//	A
  param[2]=2.;					//	B
  sigma[0]=0.07*pow(radio_stack[t][s],0.5)*0.3;


  FILE * NOM;
  char nomAr[410];
  sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",first_theory,last_theory,1,1);
  NOM=fopen(nomAr,"r");
  fscanf(NOM,"%d %le\n",&n_salta,&co);
  for(int i=0;i<Npbd;i++)
    fscanf(NOM,"%le ",&param_bd[i]);
  fclose(NOM);

  param[1]=param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+0];
  param[2]=param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+1];
  param[0]=param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+2];

  sigma[1]=0.1;
  sigma[2]=0.1;

  if(segunda_vez==1){
    FILE * NOM;
    char nomAr[410];
    sprintf(nomAr,"lista_%c_%d.txt",NameTheories[t],s);
    NOM=fopen(nomAr,"r");
    while(fscanf(NOM,"%d %le %le %le %le\n",&n_salta,&co,&param[0],&param[1],&param[2])!=EOF);
    fclose(NOM);
  }
  printf("co=%le\n",co);
}









//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void montecarlo(Crandom & ran){
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
  int i;
  for(i=0;i<Np;i++){
    delta[i]=ran.gauss(0.,sigma[i]);
    param[i]+=delta[i];
  }
  while(param[1]>20.){
    param[1]-=delta[1];    
    delta[1]=ran.gauss(0.,sigma[1]);
    param[1]+=delta[1];
    }
  while(param[2]>30.){
    param[2]-=delta[2];    
    delta[2]=ran.gauss(0.,sigma[2]);
    param[2]+=delta[2];
    }
}






// no es necesario, el MontePithon los encuentra en la cadena
//guarda os valores de los parametros para los cuales el chi2 es maximo y minimo 
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void MaxChi(void){
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
    int i;
    if(ChiNew<ChiMin){
      for(i=0;i<Np;i++)
      pMin[i]=param[i];			//	parametros
      ChiMin=ChiNew;			//	chi_minimo
    }
    else if(ChiNew>ChiMax){
      ChiMax=ChiNew;
      for(i=0;i<Np;i++)
	pMax[i]=param[i];}
    else;
}





//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void retrocede(void){
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
  for(int i=0;i<Np;i++)
    param[i]-=delta[i];
  ChiNew=ChiOld;
//cout<<ChiNew<<"   ";
}










//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
void media_error(int t, int s){
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

  int i;

  printf("mean_error...\n");fflush(stdout);

  FILE * ES;
  char nomAr[500];
  double leed[Np],c2;
  int peso,j;

  for(i=0;i<Np;i++){			//	analisis
    mean[i]=0.;
    stdev[i]=0.;
  }
  contador=0;

  sprintf(nomAr,"lista_%c_%d.txt",NameTheories[t],s);
  ES=fopen(nomAr,"r");

  for(j=0;j<NP;j++){
    fscanf(ES,"%d %le",&peso,&ChiNew);
    contador++;//=peso;
    for(i=0;i<Np;i++)
      fscanf(ES," %le",&leed[i]);
    fscanf(ES,"\n");


    for(i=0;i<Np;i++){
      mean[i]+=leed[i];//*peso;
    }

    if(j==0) ChiMin=ChiNew;

    if(ChiNew<ChiMin){
      ChiMin=ChiNew;				//	chi_minimo
      for(i=0;i<Np;i++){
        pMin[i]=leed[i];			//	parametros
       }

    }

  }

  for(i=0;i<Np;i++)
    mean[i]/=contador;

  rewind(ES);
  for(j=0;j<NP;j++){
    fscanf(ES,"%d %le",&peso,&c2);
    for(i=0;i<Np;i++)
      fscanf(ES," %le",&leed[i]);
    fscanf(ES,"\n");

    for(i=0;i<Np;i++){
      stdev[i]+=pow(leed[i]-mean[i],2);//*peso;
    }

  }
  for(i=0;i<Np;i++)
    stdev[i]/=contador;

  fclose(ES);

  //	sirve para graficar los parametros
  FILE * SA;
  sprintf(nomAr,"ajuste_BIAS_%c.txt",NameTheories[t]);
  if(s==first_stack)
    SA=fopen(nomAr,"w+");
  else
    SA=fopen(nomAr,"a+");

  //	mean best-fit, solo bias
  for(i=0;i<1;i++)
    fprintf(SA,"%le %le %le %le\n",radio_stack[t][s],pMin[i],mean[i],sqrt(stdev[i]));
  fclose(SA);

  sprintf(nomAr,"ajuste_S1_%c.txt",NameTheories[t]);
  if(s==first_stack)
    SA=fopen(nomAr,"w+");
  else
    SA=fopen(nomAr,"a+");

  //	mean best-fit, s1
  for(i=1;i<2;i++)
    fprintf(SA,"%le %le %le %le\n",radio_stack[t][s],pMin[i],mean[i],sqrt(stdev[i]));
  fclose(SA);

  sprintf(nomAr,"ajuste_S2_%c.txt",NameTheories[t]);
  if(s==first_stack)
    SA=fopen(nomAr,"w+");
  else
    SA=fopen(nomAr,"a+");

  //	mean best-fit, s2
  for(i=2;i<3;i++)
    fprintf(SA,"%le %le %le %le\n",radio_stack[t][s],pMin[i],mean[i],sqrt(stdev[i]));
  fclose(SA);

 


  printf("mean_error...hecho\n");fflush(stdout);


}












double Suppress_(double y){
return exp(-param[1]*pow(y,-param[2]));
}



//	calculados usando cosmo8b_desiende_bias_densidad.c
double alpha=7.593649e+00;
double beta=9.381333e+00;
double ca=1.415600e+00;
double cb=1.106479e+00;



double PerfilUniversal_(double r){
return -dc*(1.-pow(r/ca,alpha))/(1.+pow(r/cb,beta));
}





//Programa principal

int main(){
  Crandom ran(764);
  int i,l;
  double cn;

  FILE * NOM;
  FILE * OML;
  char nomAr[210];

  LeePerfiles();	//	read data
  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100

  calcula_correlation_fuction=1;	//calcula la funcion de correlacion cuando entre a main_rodrigo
  imprime_correlation_fuction=1;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory<3)	{param_cosmo[0]=log10(parametro[id_theory]);	MG=1;}
    else		{param_cosmo[0]=parametro[id_theory];		MG=2;}
    main_rodrigo(id_theory);	//	calcula P(k,z=0), los sigma de perfiles y la funcion de correlacion materia, usa param_cosmo[0]
  }


for(int t=0;t<=6;t++)
for(int s=first_stack;s<=last_stack;s++){

  carga(t,s);		//	carga: parametros iniciales mcmc bias+densidad
  CargaCorrelacion(t);//	pone en xi_actual la correlacion correpondiente a id_theory
  Spline_interp correl(radio_xi,xi_actual);

  //	empieza calculo del chi2
  double chi2=0.0;
  for(i=0;i<bines_radio_densidad[s];i++)
    if((i!=19)&&(i!=20))
      chi2+=pow( (  PerfilUniversal_(radio_densidad[i]) + param[0] * correl.interp(radio_densidad[i]*radio_stack[t][s]) *  Suppress_(radio_densidad[i]) -  densidad[t][s][i] ) / error_densidad[t][s][i] ,2.0);



  ChiOld=ChiNew;
  ChiNew=chi2;
  n_salta++;
  //	termina calculo del chi2
  ChiOld=ChiNew;
  ChiMax=ChiNew;
  ChiMin=ChiNew;
  //	termina funcion inicie

printf("chi2=%le\n",chi2);



  if(solo_analiza==0){
  sprintf(nomAr,"lista_%c_%d.txt",NameTheories[t],s);
  NOM=fopen(nomAr,"w+");
  cout<<"comienza la cadena "<<endl;
  l=0;
  while(l<NP+BI){// numero de puntos a calcular en el espacio de fase

  montecarlo(ran);

  //	empieza calculo del chi2
  chi2=0.0;

  for(i=0;i<bines_radio_densidad[s];i++)
    if((i!=19)&&(i!=20))
      chi2+=pow( (  PerfilUniversal_(radio_densidad[i]) + param[0] * correl.interp(radio_densidad[i]*radio_stack[t][s]) *  Suppress_(radio_densidad[i]) -  densidad[t][s][i] ) / error_densidad[t][s][i] ,2.0);

if((t==0)&&(s==6)&&(param[2]<1.)) chi2*=2.;

if((t==3)&&(s==4)&&(param[1]<10.)) chi2*=2.;

  ChiOld=ChiNew;
  ChiNew=chi2;
  n_salta++;
  //	termina calculo del chi2
  //	termina montecarlo

    // Metropolis
    cn=ChiNew;
    co=ChiOld;

    if(cn<=co) {
      if(l>BI){
        fprintf(NOM,"%d	%le	",n_salta+1,cn);
        for(int k=0;k<Np;k++)
          fprintf(NOM,"%le	",param[k]);
        fprintf(NOM,"\n");
      }
     // MaxChi();
      l++;
      n_salta=0;
    }

    else {
      if((exp(0.5*(co-cn)))>=ran.r()) {																			//sin priors
        if(l>BI){
          fprintf(NOM,"%d	%le	",n_salta+1,cn);
          for(int k=0;k<Np;k++)
            fprintf(NOM,"%le	",param[k]);
        fprintf(NOM,"\n");
      }
      //  MaxChi();
        l++;
        n_salta=0;
      }
      else{
        retrocede();
      }
    }

  fflush(NOM);
//  printf("\ni=%d co=%f cn=%f\n\n",i,co,cn);

  }



fclose(NOM);
cout<<endl<<" numero de puntos ="<<l<<endl;

  }		//	if solo analiza

cout<<endl<<" calcula la media y el error en los parametros "<<endl;
media_error(t,s);


}



cout<<"terminó"<<endl;
  return 0;
}
















