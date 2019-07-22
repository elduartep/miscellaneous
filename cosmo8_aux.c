//	tomado de cosmo6.c
//	imprime todo lo necesario para graficar los mejores ajustes

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"
#include "../l/cosmo8_parametros.h"


#include "../l/cosmo8_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_abundancia.h"//	subrutinas abundancia
#include "../l/cosmo8_densidad_bias.h"	//	subrutinas bias y densidad



int NP=500;

double chi2;
char NameCase[]={'f','s'};




//Programa principal

int main(int argc,char **argv){
  if (argc != 4) {
    printf("Wrong number of arguments.\n");
    printf("si_densidad (0 o 1)\n");
    printf("si_abundancia (0 o 1)\n");
    printf("si_bias (0 o 1)\n");
  exit(0);
  }


  if (sscanf(argv[1],"%d",&si_densidad) == 0) {
    printf("1 si quiere usar datos de densidad, 0 si no\n");
    exit(0);
  }
  if (sscanf(argv[2],"%d",&si_abundancia) == 0) {
    printf("1 si quiere usar datos de abundancia, 0 si no\n");
    exit(0);
  }
  if (sscanf(argv[3],"%d",&si_bias) == 0) {
    printf("1 si quiere usar datos del bias lineal, 0 si no\n");
    exit(0);
  }

  char nomAr[210];
  FILE *CHI;
//	normal
  sprintf(nomAr,"chi2_%d%d%d.dat",si_densidad,si_abundancia,si_bias);
//	cruzado
//  sprintf(nomAr,"chi2c_%d%d%d.dat",si_densidad,si_abundancia,si_bias);
  CHI=fopen(nomAr,"w+");



  LeePerfiles();
  carga_espectro_99();
  LeeAbundancia();
  LeeBias();

  for(MG=1;MG<3;MG++){

  CargaAjustesBiasDensidad();
  CargaAjustesAbundancia();
  CargaAjustesBias();

  FILE * SA;
  if(MG==1){
    first_theory=0;
    sprintf(nomAr,"ajuste_cosmo5_fin_f.dat");}
  else{
    first_theory=3;
    sprintf(nomAr,"ajuste_cosmo5_fin_s.dat");}
  SA=fopen(nomAr,"w+");

  for(int id_theory=(MG-1)*3;id_theory<=(MG-1)*3+3;id_theory++){

  int i,l;

  FILE * NOM;
  FILE * OML;
  char nomAr[210];




// modelos cruzados
//  sprintf(nomAr,"lista2_%d_%d_%d_%d_%d.txt",id_theory,MG,si_densidad,si_abundancia,si_bias);
// modelos correctos
  sprintf(nomAr,"lista8_%d_%d_%d_%d_%d.txt",id_theory,MG,si_densidad,si_abundancia,si_bias);
  NOM=fopen(nomAr,"r");

  double chi_min,param_min;
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2,&param_cosmo[0]);
    if(i==0){
      chi_min=chi2;
      param_min=param_cosmo[0];
    }
    if(chi_min>chi2){
      chi_min=chi2;
      param_min=param_cosmo[0];
    }
  }
  fclose(NOM);

// modelos cruzados
//  fprintf(CHI,"chi2min2_%c_%c=%le\n",NameCase[MG-1],NameTheories[id_theory],chi_min);
// modelos correctos
  fprintf(CHI,"chi2min_%c_%c=%le\n",NameCase[MG-1],NameTheories[id_theory],chi_min);

  param_cosmo[0]=param_min;

  main_rodrigo(id_theory);


  if(si_densidad==1){
  for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
    ParametrosPerfilCosmo(id_theory,id_stack);
    fprintf(SA,"d=%le\n",dc);
    fprintf(SA,"a%c%d=%le\n",NameTheories[id_theory],id_stack,alpha);
    fprintf(SA,"b%c%d=%le\n",NameTheories[id_theory],id_stack,beta);
    fprintf(SA,"v%c%d=%le\n",NameTheories[id_theory],id_stack,B);
    fprintf(SA,"s%c%d=%le\n\n",NameTheories[id_theory],id_stack,A);
  }
  }
  if(si_abundancia==1){
  //	sirve para graficar la abundancia
  ParametrosAbundancia(id_theory);
  fprintf(SA,"A=%le\n",param_abundancia[3]);
  fprintf(SA,"a%c=%le\n",NameTheories[id_theory],gama);
  fprintf(SA,"b=%le\n",param_abundancia[4]);
  fprintf(SA,"c=%le\n\n",param_abundancia[5]);

  imprime_mejor_ajuste_abundancia(id_theory,MG);


  }
  if(si_bias==1){
  double aux=bias_cosmo(1.,id_theory);
  fprintf(SA,"b1=%le\n",b0);
  fprintf(SA,"b2%c=%le\n",NameTheories[id_theory],b3);
  fprintf(SA,"b3%c=%le\n",NameTheories[id_theory],b4);
  }


}

fclose(SA);
}

fclose(CHI);






cout<<"terminó"<<endl;
  return 0;
}









