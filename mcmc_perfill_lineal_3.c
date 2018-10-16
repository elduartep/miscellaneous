// mcmc que ajusta los parametros del perfil de densidad de los voids
// construido a partir de ~/con_masa1/H2_class.cpp





// noviembre 3 de 2016

// la idea es que ajuste los parametros libres para los bines_stack
// debo guardar las cadenas en un archivo
// debe incluir el valor de chi^2 y del número de pasos que desechó
// la idea es analizar las cadenas con MontePython, porque me dá la media y el sigma encima de las graficas

// la idea es que pueda re-correr a partir de los datos anteriores, como lo hace MontePython
// como también, usar la información del analisis para colocar un nuevo sigma_param en la segunda corrida...

// depués de haber analizado las cadenas: parametros y errores
// debo graficarlas en función del radio medio de cada stack
// dependiendo de la forma de la grafica, en la que debo confiar solo para valores grandes de R_stack
// debo ajustar una función f(R_stack) para cada parametro

// la idea, al final, es hacer lo mismo para gravedad modificada
// ver la diferencia entre los graficos de los parametros f(R_stack)
// e intentar generalizar esta función para que sea fr(R_stack,fof0)






// para las cadenas: tendré 31 bin_stack, podría usar 12 cores
// y mandar a cada uno a hacer una cadena diferente






//	Julio 6 de 2017
//	Perfiles_double.c fue alterado para imprimir en el mismo archivo de salida
//	tanto el perfil de densidad como la densidad integrada
//	entonces hay que alterar este programa para que pueda extraer los datos necesarios del
//	arquivo de salida de Perfiles_double.c



//	septiembre 18
//	ya que el calculo de los perfiles ahora estahecho es escala lineal para
//	valorar cada parte del perfil con el mismo peso
//	entonces debemos adaptar este para que se adapte a la nueva salida
//	cual es la verdadera diferencia: que cada bin_stack tiene un diferente alcance en r
//	200 bines para los voids pequeños, 150 para los medianos, y 100 para los dos grandes





//      febrero 27
//	voy a estuadiar a partir de ahora solo los voids esfericos
//      voy a explorar tanto el perfil como el perfil integrado
//      voy a ver solo lcdm para ver bien que relacion hay entre los parametros del perfil
//      voy a usar varias parametrizaciones de la misma función

//      voy a usar también una potencia adicional


//	marzo 2 2018
//	pequeño ajuste para que calcule el delta_c medio

#include <iostream>
#include <fstream>
#include <cmath>
#include "../l/RANDOM.h"
#include "param.h"
//NR3
#include "../l/nr3.h"
#include "../l/interp_1d.h"

  using namespace std;

  #define N 15 						//	número de parámetros
  #define M 5						//	parámetros independientes


  //	solo uno de los tres debe ser igual a 1, los otros deben ser iguales a 0
  #define perfil_completo 1
  #define perfil_interno 0
  #define perfil_externo 0

  //	1 si delta_c es constante
  #define deltac_constante 1
  #define uniformiza_error 1


  #define beta_propt_alpha 1
  #define beta_propto_alpha 0


  //	primero debo correr interno para saber cual es el delta_c
  //	luego, al correr completo, este re-escribe los resultados de interno incluyendo los valores de concentarción

  const	int esferas=1;
  const int solo_analiza=0;
  const int densidad_integrada=0;

  const int bin_stack=10;
  const int first_stack=1;
  const int last_stack=7;

  const int NP=10000;					//	numero de puntos de la cadena
  const int BI=1000;					//	burn-in

  const double pi=4.*atan(1.);				//	numero pi
  const int bin_radio=200;			//	maximo número de bines
  const int bines_radio[bin_stack]={200,200,200,200,150,150,150,150,100,100};

  double ChiMax, ChiMin, param[M],delta[M],sigma[M]; 	//	variables basicas
  double pMax[N], pMin[N], ChiOld, ChiNew;		//	variables auxiliares
  double radio_stack[bin_stack];
  double radio[bin_radio];
  double densidad[bin_stack][bin_radio];
  double error[bin_stack][bin_radio];

  int n_salta;				//	pasos que desecha antes de colocar un nuevo eslabón
  int id_stack;				//	indice del stack

  double mean[N],stdev[N]; int contador;	//	analisis

  double deltac[bin_stack],sigma_deltac[bin_stack];

inline double perfil(double r){

double beta=param[4];

#if beta_propt_alpha>0
    beta=param[4]*param[3];
#endif

#if beta_menos_alpha>0
    beta=param[3]+param[4];
#endif

#if perfil_interno>0
   return -param[0]*(1.-param[1]*pow(r,param[3]));
#else
   return -param[0]*(1.-param[1]*pow(r,param[3]))/(1.+pow(r/param[2],beta));
#endif
}




// lee el archivo que contiene los perfiles
void LeePerfiles(void){
  int i,j;
  double aux;
  cout<<"Leyendo perfiles...";
  FILE * Ud;
	if (esferas==1)
	  Ud=fopen("densidad_n_esferas_512x256_03.txt","r");
	else
	  Ud=fopen("densidad_n_union_512x256_03.txt","r");

  //	primer renglon: los radios stack en Mpc/h
	char caracter;
	fscanf(Ud,"%s ", &caracter);
  for(j=0;j<bin_stack-1;j++)
    fscanf(Ud,"%lf ", &radio_stack[j]);
  fscanf(Ud,"%lf\n", &radio_stack[j]);

  //	los demás renglones
  for(i=0;i<bin_radio;i++){
    //	radio medio del bin radial i
    fscanf(Ud,"%le ",&radio[i]);
    //	las demás columnas columnas de la densidad
    for(j=0;j<bin_stack;j++)
      fscanf(Ud,"%le %le ",&densidad[j][i], &error[j][i]);

    //	las demás columnas columnas de la densidad integrada
    if(densidad_integrada==0){
    //	borde externo del bin radial i
    fscanf(Ud,"%le ",&aux);
      for(j=0;j<bin_stack-1;j++)
        fscanf(Ud,"%le %le ",&aux, &aux);
      fscanf(Ud,"%le %le\n",&aux, &aux);}
    else{
    //	borde externo del bin radial i
    fscanf(Ud,"%le ",&radio[i]);
    //	las demás columnas columnas de la densidad integrada
    for(j=0;j<bin_stack-1;j++)
      fscanf(Ud,"%le %le ",&densidad[j][i], &error[j][i]);
    fscanf(Ud,"%le %le\n",&densidad[j][i], &error[j][i]);}
  }

#if uniformiza_error>0
  for(j=0;j<bin_stack-1;j++){
    for(i=0;i<bin_radio;i++){
      if(error[j][i]<0.025){
        if(error[j][i]<0.025){
          error[j][i]=0.025;}
      }
    }
//	fijando el valor de un punto
densidad[j][0+j]=-0.867554;		///////////////////
error[j][0+j]=0.005;		///////////////////
  }
#endif

  fclose(Ud);
  cout<<"echo"<<endl;
}





// calcula el chi2
void Chi(void){
  int i;
  double chi2=0.0;
  //	el numero de bines radiales con información ahora depende del valor de bin_stack
  for(i=0;i<bines_radio[id_stack];i++){
//    if ((id_stack>1)||(((radio[i]>1.1)||(radio[i]<0.7))&&(id_stack<2)))
//    if (((radio[i]>0.5)&&(id_stack>3))||((radio[i]>0.8)&&((id_stack==3)||(id_stack==2)||(id_stack==1)))||((radio[i]>1.1)&&(id_stack==0))||(todo_el_perfil==1))
if((perfil_completo==1) || ((perfil_externo==1)&&(radio[i]>1.)) || ((perfil_interno==1)&&(radio[i]<1.)) )
      chi2+=pow((perfil(radio[i])-densidad[id_stack][i])/error[id_stack][i],2.0);
  }
  ChiOld=ChiNew;
  ChiNew=chi2;
  n_salta++;
}




#if perfil_completo>0
#include "mcmc_perfil_3_completo.h"
#endif

#if perfil_interno>0
#include "mcmc_perfil_3_interno.h"
#endif

#if perfil_externo>0
#include "mcmc_perfil_3_externo.h"
#endif

// no es necesario, el MontePithon los encuentra en la cadena
//guarda os valores de los parametros para los cuales el chi2 es maximo y minimo 
void MaxChi(void){
    int i;
    if(ChiNew<ChiMin){
      for(i=0;i<M;i++)
      pMin[i]=param[i];			//	parametros
      ChiMin=ChiNew;			//	chi_minimo
    }
    else if(ChiNew>ChiMax){
      ChiMax=ChiNew;
      for(i=0;i<M;i++)
	pMax[i]=param[i];}
    else;
}




void Montecarlo(Crandom & ran){
  int i;
  for(i=0;i<M;i++){
    delta[i]=ran.gauss(0.,sigma[i]);
    param[i]+=delta[i];
  }
  Chi();					//	compara con observaciones
}




void retrocede(void){
  for(int i=0;i<M;i++)
    param[i]-=delta[i];
  ChiNew=ChiOld;
//cout<<ChiNew<<"   ";
}











void media_error(void){
  printf("mean_error para el id_stack=%d...\n",id_stack);fflush(stdout);
  FILE * SIG;
  SIG=fopen("Sigma_MG.dat","r");
  Int cont=0,i,j;
  Doub masa,sigma,rad;
  while(fscanf(SIG,"%le %le\n",&rad,&sigma)!=EOF)
    cont++;
  VecDoub Radio(cont);
  VecDoub Sigma(cont);
  rewind(SIG);
  for(i=0;i<cont;i++){
    fscanf(SIG,"%le %le\n",&rad,&sigma);
    Radio[i]=rad;
    Sigma[i]=sigma;}
  fclose(SIG);

  Spline_interp sig(Radio,Sigma);		//	interpola el sigma         sig.interp(0.002)

//for(i=0;i<cont;i++)
//printf("%le %le %le\n",Radio[i],Sigma[i],sig.interp(Radio[i]));
//for(i=0;i<bin_stack;i++)
//printf("%le\n",radio_stack[i]);


  FILE * ES;
  char nomAr[500];
  double a[N],c2;
  int peso;

  for(i=0;i<N;i++){			//	analisis
    mean[i]=0.;
    stdev[i]=0.;
  }
  contador=0;

  sprintf(nomAr,"lista_%d.txt",id_stack);
  ES=fopen(nomAr,"r");

  for(j=0;j<NP;j++){
    fscanf(ES,"%d %lf",&peso,&ChiNew);
    contador++;//=peso;
    for(i=0;i<M;i++)
      fscanf(ES," %le",&a[i]);
    fscanf(ES,"\n");

#if beta_propt_alpha>0
    a[4]*=a[3];
#endif

#if beta_menos_alpha>0
    a[4]+=a[3];
#endif

a[i]=(a[4]-a[3]);			i++;		//	beta-alpha
a[i]=pow(abs(a[1]),-1./a[3]);		i++;		//	rs
a[i]=pow(a[2],-a[4]);			i++;		//	rv^-beta
a[i]=a[0]*a[1];				i++;		//	dc*A
a[i]=a[0]/a[1]*pow(a[2],a[4]);		i++;		//	dc rs^-alpha rv^beta
a[i]=pow(a[2],-a[4])/a[0];		i++;		//	dc^-1 rv^-beta
a[i]=pow(abs(a[1]),-1./a[3])/a[2];	i++;		//	rs/rv
a[i]=1./a[2];				i++;		//	1/rv
a[i]=a[1]*pow(a[2],a[3]);		i++;		//	(rv/rs)^alpha
a[i]=a[4]/a[3];						//	beta/alpha

    for(i=0;i<N;i++){
      mean[i]+=a[i];//*peso;
    }

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
    fscanf(ES,"%d %lf",&peso,&c2);
    for(i=0;i<M;i++)
      fscanf(ES," %le",&a[i]);
    fscanf(ES,"\n");

#if beta_propt_alpha>0
    a[4]*=a[3];
#endif

#if beta_menos_alpha>0
    a[4]+=a[3];
#endif

a[i]=(a[4]-a[3]);			i++;		//	beta-alpha
a[i]=pow(abs(a[1]),-1./a[3]);		i++;		//	rs
a[i]=pow(a[2],-a[4]);			i++;		//	rv^-beta
a[i]=a[0]*a[1];				i++;		//	dc*A
a[i]=a[0]/a[1]*pow(a[2],a[4]);		i++;		//	dc rs^-alpha rv^beta
a[i]=pow(a[2],-a[4])/a[0];		i++;		//	dc^-1 rv^-beta
a[i]=pow(abs(a[1]),-1./a[3])/a[2];	i++;		//	rs/rv
a[i]=1./a[2];				i++;		//	1/rv
a[i]=a[1]*pow(a[2],a[3]);		i++;		//	(rv/rs)^alpha
a[i]=a[4]/a[3];						//	beta/alpha

    for(i=0;i<N;i++){
      stdev[i]+=pow(a[i]-mean[i],2);//*peso;
    }

  }
  for(i=0;i<N;i++)
    stdev[i]/=contador;

  fclose(ES);


  //	correccion de rs para
if(id_stack==6){
mean[6]=pow(abs(mean[1]),-1./mean[3]);
mean[11]=mean[6]/mean[2];
stdev[6]=(stdev[1]*pow(abs(mean[1]),-2)+stdev[3]*pow(log(abs(mean[1]))/mean[3],2))*pow(mean[6]/mean[3],2);
stdev[11]=mean[11]*mean[11]*(stdev[6]*pow(mean[6],-2)+stdev[2]*pow(mean[2],-2));
}




  //	sirve para graficar los parametros
#if perfil_interno>0
  if (esferas==1)
    sprintf(nomAr,"ajuste_esferas_%d_int.txt",id_stack);
  else
    sprintf(nomAr,"ajuste_union_%d_int.txt",id_stack);
#else
  if (esferas==1)
    sprintf(nomAr,"ajuste_esferas_%d.txt",id_stack);
  else
    sprintf(nomAr,"ajuste_union_%d.txt",id_stack);
#endif
  FILE * SA;
  SA=fopen(nomAr,"w+");


//	media
  fprintf(SA,"%le %le ",radio_stack[id_stack],parametro);
  for(i=0;i<N;i++)
    fprintf(SA,"%le %le ",mean[i],sqrt(stdev[i]));
  fprintf(SA,"%le %le\n",sig.interp(radio_stack[id_stack]),sig.interp(radio_stack[id_stack]*mean[2]));

//	best fit
  fprintf(SA,"%le %le ",radio_stack[id_stack],parametro);
  for(i=0;i<N;i++)
    fprintf(SA,"%le %le ",pMin[i],sqrt(stdev[i]));
  fprintf(SA,"%le %le\n",sig.interp(radio_stack[id_stack]),sig.interp(radio_stack[id_stack]*mean[2]));
  fclose(SA);

  //	sirve para graficar los perfiles
#if perfil_interno>0
  if(esferas==1)
    sprintf(nomAr,"ajuste_esferas_int.txt");
  else
    sprintf(nomAr,"ajuste_union_int.txt");
#else
  if(esferas==1)
    sprintf(nomAr,"ajuste_esferas.txt");
  else
    sprintf(nomAr,"ajuste_union.txt");

#endif

  if((id_stack==0)||(first_stack==last_stack))
    SA=fopen(nomAr,"w+");
  else
    SA=fopen(nomAr,"a");
  fprintf(SA,"\nradio=%le\n parametro=%le\n",radio_stack[id_stack],parametro);
  fprintf(SA,"a%c%d=%le\n da%c%d=%le\n",nombre,id_stack,mean[0],nombre,id_stack,sqrt(stdev[0]));
  fprintf(SA,"b%c%d=%le\n db%c%d=%le\n",nombre,id_stack,mean[1],nombre,id_stack,sqrt(stdev[1]));
  fprintf(SA,"c%c%d=%le\n dc%c%d=%le\n",nombre,id_stack,mean[2],nombre,id_stack,sqrt(stdev[2]));
  fprintf(SA,"d%c%d=%le\n dd%c%d=%le\n",nombre,id_stack,mean[3],nombre,id_stack,sqrt(stdev[3]));
  fprintf(SA,"e%c%d=%le\n de%c%d=%le\n",nombre,id_stack,mean[4],nombre,id_stack,sqrt(stdev[4]));
/*  fprintf(SA,"\nradio=%le\n parametro=%le\n\n",radio_stack[id_stack],parametro);
  fprintf(SA,"a%c%d=%le\n da%c%d=%le\n",nombre,id_stack,pMin[0],nombre,id_stack,sqrt(stdev[0]));
  fprintf(SA,"b%c%d=%le\n db%c%d=%le\n",nombre,id_stack,pMin[1],nombre,id_stack,sqrt(stdev[1]));
  fprintf(SA,"c%c%d=%le\n dc%c%d=%le\n",nombre,id_stack,pMin[2],nombre,id_stack,sqrt(stdev[2]));
  fprintf(SA,"d%c%d=%le\n dd%c%d=%le\n",nombre,id_stack,pMin[3],nombre,id_stack,sqrt(stdev[3]));
  fprintf(SA,"e%c%d=%le\n de%c%d=%le\n",nombre,id_stack,pMin[4],nombre,id_stack,sqrt(stdev[4]));*/

  fprintf(SA,"\n");
  fclose(SA);



//	sirve para calcular deltac medio de todos los bines_stack
deltac[id_stack]=mean[0];
sigma_deltac[id_stack]=sqrt(stdev[0]);







  printf("mean_error para el id_stack=%d...hecho\n",id_stack);fflush(stdout);


}









void calcula_media_deltac(){
  double media_deltac=0.;
  double pesos_deltac=0.;  
  for (id_stack=first_stack;id_stack<=last_stack;id_stack++){
    media_deltac+=deltac[id_stack]/sigma_deltac[id_stack];
    pesos_deltac+=1./sigma_deltac[id_stack];
  }
  media_deltac/=pesos_deltac;
  printf("media de los deltac de los stack [%d,%d]=%le\n",first_stack,last_stack,media_deltac);
}










//Programa principal

int main(){
  Crandom ran2(764);
  int i,l;
  double cn, co;

  FILE * NOM;
  FILE * OML;
  char nomAr[21];

  LeePerfiles();


for (id_stack=first_stack;id_stack<=last_stack;id_stack++){


  if(solo_analiza==0){
  sprintf(nomAr,"lista_%d.txt",id_stack);
  NOM=fopen(nomAr,"w+");

  Inicie();
  cout<<"comienza la cadena del stack "<<id_stack<<endl;
  i=0;
  while(i<NP+BI){// numero de puntos a calcular en el espacio de fase

    Montecarlo(ran2);

    // Metropolis
    cn=ChiNew;
    co=ChiOld;

    if(cn<=co) {
      if(i>BI){
        fprintf(NOM,"%d	%.2f	",n_salta+1,cn);
        for(int k=0;k<M;k++)
          fprintf(NOM,"%le	",param[k]);
        fprintf(NOM,"\n");
      }
      MaxChi();
      i++;
      n_salta=0;
    }

    else {
      if((exp(0.5*(co-cn)))>=ran2.r()) {																			//sin priors
        if(i>BI){
          fprintf(NOM,"%d	%.2f	",n_salta+1,cn);
          for(int k=0;k<M;k++)
            fprintf(NOM,"%le	",param[k]);
        fprintf(NOM,"\n");
      }
        MaxChi();
        i++;
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
cout<<endl<<" numero de puntos ="<<i<<endl;

  }		//	if solo analiza

cout<<endl<<" calcula la media y el error en los parametros "<<endl;
media_error();


}


calcula_media_deltac();





cout<<"terminó"<<endl;
  return 0;
}







