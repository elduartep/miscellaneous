#include <iostream>
#include <fstream>
#include <cmath>
#include "../l/RANDOM.h"
//#include "param.h"
//NR3
#include "../l/nr3.h"
#include "../l/interp_1d.h"

  using namespace std;
  #define seed 2176

  const int segunda_vez=0;
  const int solo_analiza=0;
  const int norma_error=0;


  #define N 9
  const int without_id_theory=40;

  const int theories=7;
  int first_theory=3;
  int  last_theory=6;

  char NameTheories[theories]={'4','5','6','l','A','B','D'};
  const double parametro[theories]={1.000000e-04,1.000000e-05,1.000000e-06,1.000000e-08,1.000000e-00,2.000000e-00,3.000000e-00};

  const int bin=15;
  const int first_bin=0;
  const int last_bin=12;

  const int NP=500000;				//	numero de puntos de un parametro
  const int BI=100;				//	numero de puntos de un parametro
  const int ciclo_lineal=100;
  const int tercio_ciclo_lineal=70;
  const double disminuye=0.869;
  const double aumenta=1.15;




  double valor_param[NP];
  double medial;				//	la media del parametro que estoy variando
  double stdevl;				//	variacion del parametro
  double fac=0.3;
  double ChiMin1;
  double pMin1;

  const double pi=4.*atan(1.);			//	numero pi

  double ChiMax, ChiMin, param[N],delta[N],sigma[N]; 	//	variables basicas
  double pMax[N], pMin[N], ChiOld, ChiNew;		//	variables auxiliares

  double delta_param[N],param_old[N];	//	para evolucion lineal

  double alpha;
  double cn, co;

  double Sigma[theories][bin];
  double medida[theories][bin];
  double error[theories][bin];

  int n_salta;				//	pasos que desecha antes de colocar un nuevo eslabón

  double mean[N],stdev[N]; int contador;	//	analisis










void carga(void){
  //	alpha: potencia
  param[0]=1.45e-00;	sigma[0]=1.e-01;	//	4
  param[1]=1.30e-00;	sigma[1]=1.e-01;	//	5
  param[2]=1.15e-00;	sigma[2]=1.e-01;	//	6
  param[3]=1.4e-00;	sigma[3]=1.e-01;	//	l
  param[4]=1.2e-00;	sigma[4]=1.e-01;	//	A
  param[5]=1.6e-00;	sigma[5]=1.e-01;	//	B
  param[6]=1.8e-00;	sigma[6]=1.e-01;	//	D

  param[7]=4.7e-01;	sigma[7]=1.0e-02;	//	A: factor global
  param[8]=1.9e-01;	sigma[8]=1.0e-02;	//	c: factor dentro de la exponencial

  //	A sigma^-alpha ( 1 + b nu^beta) exp( -c nu^gamma )

if(segunda_vez==1){
  FILE * NOM;
  char nomAr[210];
  sprintf(nomAr,"lista_abundancia_%d_%d.txt",first_theory,last_theory);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n"
   ,&n_salta,&co,&param[0],&param[1],&param[2],&param[3],&param[4],&param[5],&param[6],&param[7],&param[8],&param[9],&param[10],&param[11],&param[12],&param[13],&param[14],&param[15])!=EOF);
  fclose(NOM);
  cn=co;
}
printf("co=%le\n",co);
}






void MonteCarlo(Crandom & ran){
for(int i=0;i<N;i++){
  delta[i]=ran.gauss(0.,sigma[i])*fac;
  param[i]+=delta[i];//printf("%d %d,,, ",i,j);
}
}











void retrocede(void){
for(int i=0;i<N;i++){
  param[i]-=delta[i];
  ChiNew=ChiOld;
}
}








//	A nu^alpha ( 1 + b nu^beta) exp( -c nu^gamma )
inline double FuncionRadio(double sigma){
//double nu=1.686/sigma;
return param[7]*pow(sigma,-alpha)*exp(-param[8]*pow(sigma,-4.5));
}








// lee el archivo que contiene la función de radio que fué calculada por FuncionRadioPredice.c
void LeeFuncionRadio(void){
  printf("Leyendo funcion radio...");fflush(stdout);
  int i;
  char NomArch[300];
  char line[450];
  FILE * Ud;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    sprintf(NomArch,"../%c/f_medida.dat",NameTheories[id_theory]);

    Ud=fopen(NomArch,"r");
    //	primer renglon: nombres de las columnas
    fgets(line, 450, Ud);

    //	los demás renglones
    double f_max,f_min;
    for(i=0;i<bin;i++){
      fscanf(Ud,"%le %le %le %le\n",&Sigma[id_theory][i], &medida[id_theory][i], &f_max, &f_min);
      error[id_theory][i]=f_max-f_min;
      if((error[id_theory][i]<0.004)&&(norma_error==1))
        error[id_theory][i]=0.004;
    }
    fclose(Ud);
  }
  printf("echo\n");fflush(stdout);
}




















void media_error(void){
  printf("mean_error...\n");fflush(stdout);

  FILE * ES;
  char nomAr[500];
  double a[N],c2;
  int peso,i,j;

  for(i=0;i<N;i++){			//	analisis
    mean[i]=0.;
    stdev[i]=0.;
  }
  contador=0;

  sprintf(nomAr,"lista_abundancia2_%d_%d.txt",first_theory,last_theory);
  ES=fopen(nomAr,"r");

  for(j=0;j<NP;j++){
    fscanf(ES,"%d %le",&peso,&ChiNew);
    contador++;//=peso;

    for(i=0;i<N;i++){
      fscanf(ES," %le",&a[i]);
      mean[i]+=a[i];//*peso;
    }
    fscanf(ES,"\n");

    if(j==0) ChiMin=ChiNew;

    if(ChiNew<ChiMin){
      ChiMin=ChiNew;				//	chi_minimo
      for(i=0;i<N;i++){
        pMin[i]=a[i];			//	parametros
       }

    }

  }

  for(i=0;i<N;i++)
    mean[i]/=contador;

  rewind(ES);
  for(j=0;j<NP;j++){
    fscanf(ES,"%d %le",&peso,&c2);
    for(i=0;i<N;i++){
      fscanf(ES," %le",&a[i]);
      stdev[i]+=pow(a[i]-mean[i],2);//*peso;
    }
    fscanf(ES,"\n");

  }
  for(i=0;i<N;i++)
    stdev[i]/=contador;

  fclose(ES);




  //	sirve para graficar los parametros
  FILE * SA;
  sprintf(nomAr,"ajuste_abundancia2f_%d_%d.txt",first_theory,last_theory);
  SA=fopen(nomAr,"w+");
//	media
  for(i=first_theory;i<=3;i++)
    fprintf(SA,"%le %le %le %le %le %le %le\n",parametro[i]
     ,mean[i],sqrt(stdev[i]),mean[7],sqrt(stdev[7]),mean[8],sqrt(stdev[8]));
//	best fit
  for(i=first_theory;i<=3;i++)
    fprintf(SA,"%le %le %le %le %le %le %le\n",parametro[i]
     ,pMin[i],sqrt(stdev[i]),pMin[7],sqrt(stdev[7]),pMin[8],sqrt(stdev[8]));
  fclose(SA);
  sprintf(nomAr,"ajuste_abundancia2s_%d_%d.txt",first_theory,last_theory);
  SA=fopen(nomAr,"w+");
//	media
  for(i=3;i<=last_theory;i++)
    fprintf(SA,"%le %le %le %le %le %le %le\n",parametro[i]
     ,mean[i],sqrt(stdev[i]),mean[7],sqrt(stdev[7]),mean[8],sqrt(stdev[8]));
//	best fit
  for(i=3;i<=last_theory;i++)
    fprintf(SA,"%le %le %le %le %le %le %le\n",parametro[i]
     ,pMin[i],sqrt(stdev[i]),pMin[7],sqrt(stdev[7]),pMin[8],sqrt(stdev[8]));
  fclose(SA);

  //	sirve para graficar los perfiles
  sprintf(nomAr,"ajuste_abundancia2_%d_%d.dat",first_theory,last_theory);

  SA=fopen(nomAr,"w+");

  fprintf(SA,"parametro=%le\n",parametro[first_theory]);
for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
  fprintf(SA,"alpha%c=%le\n d_alpha%c=%le\n"
	,NameTheories[id_theory],pMin[id_theory],NameTheories[id_theory],sqrt(stdev[id_theory]));
}
  fprintf(SA,"A=%le\n d_A=%le\n",pMin[7],sqrt(stdev[7]));
  fprintf(SA,"c=%le\n d_c=%le\n",pMin[8],sqrt(stdev[8]));
  fprintf(SA,"\n");
  fclose(SA);

  printf("chi2_min=%le\n",ChiMin);

  printf("mean_error...hecho\n");fflush(stdout);


}






















//Programa principal

//Programa principal
int main(int argc,char **argv){


  Crandom ran(seed);
  int i,j,l;

  FILE * NOM;
  FILE * OML;
  char nomAr[210];

  LeeFuncionRadio();
  carga();			//	varianza de parametros mcmc



  //	empieza calculo del chi2
  double chi2=0.0;
  double valor;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory!=without_id_theory){
    alpha=param[id_theory];
    for(i=first_bin;i<=last_bin;i++)
      chi2+=pow((FuncionRadio(Sigma[id_theory][i])-medida[id_theory][i])/error[id_theory][i],2.0);
  }}
  ChiNew=chi2;
  n_salta++;
  ChiOld=ChiNew;
  ChiMax=ChiNew;
  ChiMin=ChiNew;
  //	termina funcion inicie

printf("chi2=%le\n",chi2);





  if(solo_analiza==0){
  sprintf(nomAr,"lista_abundancia2_%d_%d.txt",first_theory,last_theory);
  if(segunda_vez==1)	NOM=fopen(nomAr,"a+");
  else			NOM=fopen(nomAr,"w+");

  cout<<"comienza la cadena "<<endl;

  l=0;
  n_salta=0;
  while(l<NP+BI){// numero de puntos a calcular en el espacio de fase

  MonteCarlo(ran);

  //	empieza calculo del chi2
  chi2=0.0;  
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory!=without_id_theory){
    alpha=param[id_theory];
    for(i=first_bin;i<=last_bin;i++)
      chi2+=pow((FuncionRadio(Sigma[id_theory][i])-medida[id_theory][i])/error[id_theory][i],2.0);
  }}
  ChiOld=ChiNew;
  ChiNew=chi2;
  n_salta++;


  // Metropolis
  cn=ChiNew;
  co=ChiOld;

  if(ChiMin>ChiNew)
    ChiMin=ChiNew;
  co=ChiMin;	//	comparo conChiMin para la exponencial en vez del chi anterior


    if(cn<=co) {
      if(l>BI){
        fprintf(NOM,"%d	%le ",n_salta,cn);
        for(int k=0;k<N;k++)
          fprintf(NOM,"%le ",param[k]);
        fprintf(NOM,"\n");
      }
      l++;
      n_salta=0;
    }

    else {
      if((exp(0.5*(co-cn)))>=ran.r()) {
        if(l>BI){
          fprintf(NOM,"%d %le ",n_salta,cn);
          for(int k=0;k<N;k++)
            fprintf(NOM,"%le	",param[k]);
        fprintf(NOM,"\n");
      }
        l++;
        n_salta=0;
      }
      else{
        retrocede();
      }
    }

  fflush(NOM);

 // printf("l=%d lr=%d co=%f cn=%f\n\n",l,lr,co,cn);

  }	//	while l<NP




fclose(NOM);
cout<<endl<<" numero de puntos ="<<l<<endl;

  }		//	if solo analiza

cout<<endl<<" calcula la media y el error en los parametros "<<endl;
media_error();






cout<<"terminó"<<endl;
  return 0;
}







