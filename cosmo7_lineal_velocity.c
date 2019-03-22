
#include <iostream>
#include <fstream>
#include <cmath>
#include "../l/RANDOM.h"
//#include "param.h"
//NR3
#include "../l/nr3.h"
//#include "../l/interp_1d.h"

#include "../l/cosmo7_parametros.h"

  using namespace std;








//	tomado de		mcmc_todo_libre5nuevofast_2.c
// lee el archivo que contiene los perfiles
void LeePerfiles(void){
  int i,j;
  double aux;
  printf("Leyendo perfiles...");fflush(stdout);
  char NomArch[300];
  FILE * Ud;
for(int id_theory=0;id_theory<=6;id_theory++){
  sprintf(NomArch,"../%c/densidad_n_esferas_512x256_03.txt",NameTheories[id_theory]);

  Ud=fopen(NomArch,"r");
  //	primer renglon: los radios stack en Mpc/h
  char caracter;
  fscanf(Ud,"%s ", &caracter);
  for(j=0;j<bin_stack-1;j++)
    fscanf(Ud,"%lf ", &radio_stack[id_theory][j]);
  fscanf(Ud,"%lf\n", &radio_stack[id_theory][j]);

  //	los demás renglones
  for(i=0;i<bin_radio_densidad;i++){
    //	radio medio del bin radial i
    fscanf(Ud,"%le ",&radio_densidad[i]);
    //	las demás columnas columnas de la densidad
    for(j=0;j<bin_stack;j++)
      fscanf(Ud,"%le %le ",&densidad[id_theory][j][i], &error_densidad[id_theory][j][i]);

    //	las demás columnas columnas de la densidad integrada
    if(densidad_integrada==0){
    //	borde externo del bin radial i
    fscanf(Ud,"%le ",&aux);
      for(j=0;j<bin_stack-1;j++)
        fscanf(Ud,"%le %le ",&aux, &aux);
      fscanf(Ud,"%le %le\n",&aux, &aux);}
    else{
    //	borde externo del bin radial i
    fscanf(Ud,"%le ",&radio_densidad[i]);
    //	las demás columnas columnas de la densidad integrada
    for(j=0;j<bin_stack-1;j++)
      fscanf(Ud,"%le %le ",&densidad[id_theory][j][i], &error_densidad[id_theory][j][i]);
    fscanf(Ud,"%le %le\n",&densidad[id_theory][j][i], &error_densidad[id_theory][j][i]);}
  }

#if uniformiza_densidad>0
  for(j=0;j<bin_stack-1;j++){
    for(i=0;i<bin_radio_densidad;i++){
//printf("%d %le %le %le\n",j,radio_densidad[i],densidad[id_theory][j][i],error_densidad[id_theory][j][i]);fflush(stdout);
      if(error_densidad[id_theory][j][i]<0.025){
         error_densidad[id_theory][j][i]=0.025;
      }
    }
//	fijando el valor de un punto
//densidad[id_theory][j][0+j]=-dc;		///////////////////
//error_densidad[id_theory][j][0+j]=0.005;		///////////////////
  }
#endif
  fclose(Ud);
}
  printf("echo\n");fflush(stdout);
}








void imprime_teoria_radial(){
  int i,s,t=3;
  double     integral[bin_stack]={0.};
  double integral_inf[bin_stack]={0.};
  double integral_sup[bin_stack]={0.};
  double drb=10./bin_radio_densidad;		//	delta radio bin (o sea, en unidades de r.2)
  double prefactor=-pow(Om,0.55)*H0*drb;

  FILE * AL;
  AL=fopen("teoria_radial.dat","w+");

  //	primera aproximación de v, usando dot_delta lineal
  //	itero los siguientes dos pasos
  //	con dicha velocidad estimo dot_delta
  //	con dicho dot_delta estimo velocidad
  for(i=0;i<bin_radio_densidad;i++){
    fprintf(AL,"%le ",drb*(0.5+1.*i));	// radio medio del bin radial i
    for(s=0;s<bin_stack;s++){
          integral[s]+=pow(0.5+1.*i,2)*(densidad[t][s][i]);
      integral_inf[s]+=pow(0.5+1.*i,2)*(densidad[t][s][i]-error_densidad[t][s][i]);
      integral_sup[s]+=pow(0.5+1.*i,2)*(densidad[t][s][i]+error_densidad[t][s][i]);
      if(i<bines_radio_densidad[s]){
        fprintf(AL,"%le ",prefactor*radio_stack[t][s]*pow(0.5+1.*i,-2)*integral[s]);
        fprintf(AL,"%le ",prefactor*radio_stack[t][s]*pow(0.5+1.*i,-2)*integral_inf[s]);
        fprintf(AL,"%le ",prefactor*radio_stack[t][s]*pow(0.5+1.*i,-2)*integral_sup[s]);}
    }
    fprintf(AL,"\n");
  }
  fclose(AL);

}








//Programa principal

int main(){

LeePerfiles();

//alloca_vol_bines();
 
imprime_teoria_radial();


//printf("bin_radio_densidad=%d \n",bin_radio_densidad);





  return 0;
}







