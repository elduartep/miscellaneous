//	calcula correlacion y sigma usando los valores fiduciales para cada caso f(R) y Symmetron

#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "../l/cosmo7_parametros.h"

  using namespace std;
 


 


#include "../l/cosmo7_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR


//Programa principal

int main(){

  for (int first_theory=0;first_theory<4;first_theory+=3){
  int  last_theory=first_theory+3;

if(first_theory==0)MG=1;
if(first_theory==3)MG=2;



  correlation_fuction=1;

  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(first_theory==0)param_cosmo[0]=log10(parametro[id_theory]);
    if(first_theory==3)param_cosmo[0]=parametro[id_theory];
    carga_espectro_99();	//	carga: espectro lcdm en 99 y 100
    main_rodrigo(id_theory);	//	calcula P(k,z=0) y los sigma de perfilesy la funcion de correlacion, usa param_cosmo[0]
  }

}




cout<<"terminó"<<endl;
  return 0;
}







