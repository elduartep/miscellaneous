// imprime algumas predicciones linealess
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





//Programa principal

int main(int argc,char **argv){



  calcula_correlation_fuction=1;



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




  double sup,inf;


  param_cosmo[0]=parametro[id_theory];
    main_rodrigo(id_theory);	//	calcula P(k,z=0), los sigma de perfiles y la funcion de correlacion materia, usa param_cosmo[0]

 
 

}






cout<<"terminó"<<endl;
  return 0;
}




