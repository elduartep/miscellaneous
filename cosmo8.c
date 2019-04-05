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

#include "../l/cosmo8_parametros.h"




  const int NP=100;			//	numero de puntos


  #define bias_teoria 1
  #define densidad_teoria 1



#include "../l/cosmo8_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_abundancia.h"//	subrutinas abundancia
#include "../l/cosmo8_densidad_bias.h"//	subrutinas densidad y bias (estan juntas porque se puede ajustar las dos al tiempo, pero no da buen resultado)






  double chi2;

// calcula el chi2
double Chi(int id_theory){
  double chi2=0.;
  double bias_aux,densidad_aux;
  if(si_abundancia==1){
    ParametrosAbundancia(id_theory);
    calcula_sigma_derivadas(id_theory);	//	calcula P(z=0)
    calcula_abundancia_teorica(id_theory);//	calcula sigma y derivadas
    for(int i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      chi2+=pow((teoria_abundancia[id_theory][i]-abundancia[id_theory][i])/error_abundancia[id_theory][i],2.0);
    }

    if(si_bias==1){
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        bias_aux=bias(sigma_bias[id_theory][id_stack]);
        chi2+=pow((bias_aux - medida_bias[id_theory][id_stack])/error_bias[id_theory][id_stack],2.0);
//        chi2+=pow((bias_aux - medida_bias2[id_theory][id_stack])/error_bias2[id_theory][id_stack],2.0);
//        chi2+= pow((bias_aux - (medida_bias[id_theory][id_stack] +medida_bias2[id_theory][id_stack])*0.5),2.0)/(pow(error_bias[id_theory][id_stack],2)*0.25 + pow(error_bias2[id_theory][id_stack],2)*0.25);
      }
    }

    if(si_densidad==1){
      CargaCorrelacion(id_theory);//	pone en xi_actual la correlacion correpondiente a id_theory
      Spline_interp correl(radio_xi,xi_actual);
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        if(bias_teoria==1)	bias_aux = bias(sigma_bias[id_theory][id_stack]);
        else 			bias_aux = medida_bias[id_theory][id_stack];
        ParametrosSuppress(id_theory,id_stack);
        for(int i=0;i<bines_radio_densidad[id_stack];i++){
          if(densidad_teoria==1)	densidad_aux = PerfilUniversal(radio_densidad[i]);
          else				densidad_aux = (densidad[1][5][i]+densidad[5][5][i])*0.5;
          if((i!=19)&&(i!=20))
            chi2+=pow( (densidad_aux + bias_aux * correl.interp(radio_densidad[i]*radio_stack[id_theory][id_stack]) *  Suppress(radio_densidad[i]) - densidad[id_theory][id_stack][i]) / error_densidad[id_theory][id_stack][i],2.0);
        }
      }
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
  CargaAjustesBiasDensidad();	//	densidad + bias cuando se usa el fit

  LeeBias();
  CargaAjustesBias();

  LeeAbundancia();

  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100


  for(int id_theory=0;id_theory<=1;id_theory++){

  if(id_theory<3)	MG=1;
  else		MG=2;

  CargaAjustesAbundancia();	//	todavia depende de la teoria


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
    if(id_theory==3){      inf=-10.;   sup=-7.;}
  //	bias medido, errores internos x 3
    if(id_theory==0){      inf=-4.5;   sup=-3.2;}
    if(id_theory==1){      inf=-5.7;   sup=-4.2;}
    if(id_theory==2){      inf=-7.;   sup=-5.5;}
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
    if(id_theory==3){      inf=-0.5;  sup=1.;}
    if(id_theory==4){      inf=0.3;    sup=1.8;}
    if(id_theory==5){      inf=2.6;    sup=3.8;}
    if(id_theory==6){      inf=3.9;    sup=5.;}

//	bias medido, errores internos x 3
    if(id_theory==3){      inf=-0.5;  sup=1.;}
    if(id_theory==4){      inf=0.;    sup=2.;}
    if(id_theory==5){      inf=2.5;    sup=4.;}
    if(id_theory==6){      inf=4.2;    sup=5.4;}

if((si_bias==1)&&(si_abundancia==0)&&(si_densidad==0)){
    if(id_theory==3){      inf=1.e-5;  sup=1.5;}
    if(id_theory==4){      inf=0.0;    sup=2.5;}
    if(id_theory==5){      inf=1.0;    sup=3.5;}
    if(id_theory==6){      inf=2.0;    sup=5.0;}
}

  }

  for(param_cosmo[0]=inf;param_cosmo[0]<=sup;param_cosmo[0]+=(sup-inf)/double(NP)){
    main_rodrigo(id_theory);	//	calcula P(k,z=0), los sigma de perfiles y la funcion de correlacion materia, usa param_cosmo[0]

    double chi2=Chi(id_theory);
    fprintf(NOM,"%le %le\n",chi2,param_cosmo[0]);//printf(",");fflush(stdout);
    fflush(NOM);
  }



fclose(NOM);
}






cout<<"terminó"<<endl;
  return 0;
}





