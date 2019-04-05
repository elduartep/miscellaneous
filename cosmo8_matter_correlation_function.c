//	calcula correlacion y sigma usando los valores fiduciales para cada caso f(R) y Symmetron

#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "../l/cosmo8_parametros.h"

  using namespace std;
 


 


#include "../l/cosmo8_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR


//Programa principal

int main(){

  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100
  imprime_correlation_fuction=1;

  for (int aux_first_theory=0;aux_first_theory<4;aux_first_theory+=3){
  int  aux_last_theory=aux_first_theory+3;

if(aux_first_theory==0)MG=1;
if(aux_first_theory==3)MG=2;

  for(int id_theory=aux_first_theory;id_theory<=aux_last_theory;id_theory++){
    if(aux_first_theory==0)param_cosmo[0]=log10(parametro[id_theory]);
    if(aux_first_theory==3)param_cosmo[0]=parametro[id_theory];

    main_rodrigo(id_theory);	//	calcula P(k,z=0) y los sigma de perfilesy la funcion de correlacion, usa param_cosmo[0]
  }

}




cout<<"terminó"<<endl;
  return 0;
}







