//	calcula el espectro usando el codigo de Rodrigo, lo imprime en z=0 para compararlo con CAMB

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


#include "../l/cosmo7_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR




double first_theory,last_theory;	
char NameCase[]={'f','s'};




//Programa principal

int main(int argc,char **argv){


  for(MG=1;MG<3;MG++){

  for(int id_theory=(MG-1)*3;id_theory<=(MG-1)*3+3;id_theory++){

  int i,l;


  carga_espectro_99();

  if(MG==1)
    param_cosmo[0]=log10(parametro[id_theory]);
  else
    param_cosmo[0]=parametro[id_theory];

  si_densidad=0;
  si_bias=0;
  main_rodrigo(id_theory);


  FILE * IM;
  char Ar[200];
  sprintf(Ar,"espectro_%c_%d.dat",NameTheories[id_theory],MG);
  IM=fopen(Ar,"w+");

  for(i=0;i<nk-1;i++){///	integral
    fprintf(IM,"%le %le \n",K[i],P[i]);
  }
  fclose(IM);


}


}







cout<<"terminó"<<endl;
  return 0;
}









