//	intenta recuperar fR0 o ZSSB a partir de la densidad y/o la abundancia
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "/home/cosmousp/Documentos/sao/l/RANDOM.h"

#include "/home/cosmousp/Documentos/sao/l/nr3.h"
#include "/home/cosmousp/Documentos/sao/l/interp_1d.h"

#include "../l/cosmo8_parametros.h"




  const int NP=500;			//	numero de puntos


  #define bias_teoria 1
  #define densidad_teoria 1



#include "../l/cosmo8_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_abundancia.h"//	subrutinas abundancia
#include "../l/cosmo8_densidad_bias.h"//	subrutinas densidad y bias (estan juntas porque se puede ajustar las dos al tiempo, pero no da buen resultado)






  double chi2;

// calcula el chi2
double Chi(int MG,int id_theory){
  double chi2=0.;
  double bias_aux,densidad_aux;

    if(si_bias==1){
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        bias_aux=bias_cosmo(sigma_bias[id_theory][id_stack],id_theory);
        chi2+=pow((bias_aux - medida_bias[id_theory][id_stack])/error_bias[id_theory][id_stack],2.0);
//printf("%d %le %le %le %le\n",id_stack,radio_bias[id_theory][id_stack],bias_aux,medida_bias[id_theory][id_stack],sigma_bias[id_theory][id_stack]);
      }
    }

    if(si_densidad==1){
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        ParametrosPerfilCosmo(id_theory,id_stack);
        for(int i=0;i<bines_radio_densidad[id_stack];i++)
          if((i!=19)&&(i!=20))
            chi2+=pow( (perfil(radio_densidad[i])-densidad[id_theory][id_stack][i]) /error_densidad[id_theory][id_stack][i],2.0);
      }
    }


    if(si_abundancia==1){
      ParametrosAbundancia(id_theory);
      calcula_sigma_derivadas(id_theory);
      calcula_abundancia_teorica(id_theory);	//	calcula la abundancia teorica, necesita gama
      for(int i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      chi2+=pow((teoria_abundancia[id_theory][i]-abundancia[id_theory][i])
                 /error_abundancia[id_theory][i],2.0);
    }

  
  return chi2;
}

























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
    printf("1 si quiere usar datos de bias, 0 si no\n");
    exit(0);
  }



  if(si_densidad==1)	calcula_correlation_fuction=1;



  LeePerfiles();

  LeeBias();

  LeeAbundancia();


  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100





  for(MG=1;MG<3;MG++){		//	clase de teoria
  CargaAjustesAbundancia();
  CargaAjustesBiasDensidad();	//	densidad + bias cuando se usa el fit
  CargaAjustesBias();


  for(int id_theory=(MG-1)*3;id_theory<=(MG-1)*3+3;id_theory++){//	cual teoria




  int i,l;

  FILE * NOM;
  FILE * OML;
  char nomAr[210];


  sprintf(nomAr,"lista8_%d_%d_%d_%d_%d.txt",id_theory,MG,si_densidad,si_abundancia,si_bias);
  NOM=fopen(nomAr,"w+");


  double sup,inf;
  if(MG==1){
  //	bias teoria, errores normales
    if(id_theory==0){      inf=-4.4;   sup=-3.;}
    if(id_theory==1){      inf=-5.2;   sup=-3.9;}
    if(id_theory==2){      inf=-6.6;   sup=-4.8;}
    if(id_theory==3){      inf=-8.;   sup=-5.4;}
  //	bias medido, errores internos x 3
    if(id_theory==0){      inf=-4.4;   sup=-3.;}
    if(id_theory==1){      inf=-5.5;   sup=-4.4;}
    if(id_theory==2){      inf=-7.5;   sup=-5.;}
    if(id_theory==3){      inf=-10.;   sup=-7.;}

if((si_bias==1)&&(si_abundancia==0)&&(si_densidad==0)){
    if(id_theory==0){      inf=-5.0;   sup=-2.5;}
    if(id_theory==1){      inf=-6.5;   sup=-3.5;}
    if(id_theory==2){      inf=-10.0;   sup=-4.5;}
    if(id_theory==3){      inf=-10.;   sup=-5.;}
}

  }
  else{

//	bias teoria, errores normales
    if(id_theory==3){      inf=0.;     sup=1.8;}
    if(id_theory==4){      inf=0.;    sup=1.8;}
    if(id_theory==5){      inf=1.6;    sup=2.8;}
    if(id_theory==6){      inf=1.9;    sup=4.;}

//	bias medido, errores internos x 3
    if(id_theory==3){      inf=0.;    sup=1.4;}
    if(id_theory==4){      inf=0.;    sup=1.7;}
    if(id_theory==5){      inf=1.3;    sup=2.6;}
    if(id_theory==6){      inf=2.;    sup=3.5;}

if((si_bias==1)&&(si_abundancia==0)&&(si_densidad==0)){
    if(id_theory==3){      inf=1.e-5;  sup=1.5;}
    if(id_theory==4){      inf=0.0;    sup=2.5;}
    if(id_theory==5){      inf=1.0;    sup=3.5;}
    if(id_theory==6){      inf=2.0;    sup=5.0;}
}

  }

  for(param_cosmo[0]=inf;param_cosmo[0]<=sup;param_cosmo[0]+=(sup-inf)/double(NP)){
    main_rodrigo(id_theory);	//	calcula P(k,z=0), los sigma de perfiles y la funcion de correlacion materia, usa param_cosmo[0]

    double chi2=Chi(MG,id_theory);
    fprintf(NOM,"%le %le\n",chi2,param_cosmo[0]);//printf(",");fflush(stdout);
    fflush(NOM);
  }



fclose(NOM);
}}






cout<<"terminó"<<endl;
  return 0;
}





