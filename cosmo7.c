//	intenta recuperar fR0 o ZSSB a partir de la densidad y/o la abundancia
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

#include "../l/cosmo7_parametros.h"




  const int NP=1000;			//	numero de puntos




#include "../l/cosmo7_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_abundancia.h"//	subrutinas abundancia
#include "../l/cosmo7_densidad.h"	//	subrutinas densidad








  double chi2;

// calcula el chi2
void Chi(int id_theory){
  int i;
  main_rodrigo(id_theory);
  if(si_densidad==1){
    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
      ParametrosPerfil(id_theory,id_stack);
      for(i=0;i<bines_radio_densidad[id_stack];i++)
        chi2+=pow((perfil(radio_densidad[i])-densidad[id_theory][id_stack][i])
               /error_densidad[id_theory][id_stack][i],2.0);
    }
  }
  if(si_abundancia==1){
    ParametrosAbundancia(id_theory);
    calcula_sigma_derivadas(id_theory);	//	calcula P(z=0)
    calcula_abundancia_teorica(id_theory);//	calcula sigma y derivadas
    for(i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      chi2+=pow((teoria_abundancia[id_theory][i]-abundancia[id_theory][i])/error_abundancia[id_theory][i],2.0);
    }
}

























//Programa principal

int main(int argc,char **argv){
  if (argc != 3) {
    printf("Wrong number of arguments.\n");
    printf("si_densidad (0 o 1)\n");
    printf("si_abundancia (0 o 1)\n");
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



  for(MG=1;MG<3;MG++){
  for(int id_theory=(MG-1)*3;id_theory<=(MG-1)*3+3;id_theory++){
//  for(MG=2;MG<3;MG++){
//  for(int id_theory=(MG-1)*3;id_theory<=(MG-1)*3+3;id_theory++){


  LeePerfiles();
  CargaAjustesDensidad();
  carga_espectro_99();	//	carga el espectro lcdm en 99 y 100

  LeeAbundancia();
  CargaAjustesAbundancia();


  int i,l;

  FILE * NOM;
  FILE * OML;
  char nomAr[210];


  sprintf(nomAr,"lista_%d_%d_%d_%d.txt",id_theory,MG,si_densidad,si_abundancia);
  NOM=fopen(nomAr,"w+");


  double sup,inf;
  if(MG==1){
    if(id_theory==0){      inf=-4.5;   sup=-3.5;}
    if(id_theory==1){      inf=-5.5;   sup=-4.5;}
    if(id_theory==2){      inf=-6.5;   sup=-5.5;}
    if(id_theory==3){      inf=-10.;   sup=-6.;}
  }
  else{
    if(id_theory==3){      inf=1.e-5;  sup=0.5;}
    if(id_theory==4){      inf=0.5;  sup=1.5;}
    if(id_theory==5){      inf=1.5;  sup=2.5;}
    if(id_theory==6){      inf=2.5;  sup=3.5;}
  }

  for(param_cosmo[0]=inf;param_cosmo[0]<=sup;param_cosmo[0]+=(sup-inf)/double(NP)){
    chi2=0.0;
    Chi(id_theory);				//	compara con observaciones
    fprintf(NOM,"%le %le\n",chi2,param_cosmo[0]);//printf(",");fflush(stdout);
    fflush(NOM);
  }



fclose(NOM);
}}






cout<<"terminó"<<endl;
  return 0;
}





