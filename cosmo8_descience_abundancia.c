#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "../l/cosmo8_parametros.h"

  using namespace std;
  #define seed 2176
  #define sigma_propto_param 1

  const int segunda_vez=0;

  const int without_id_theory=40;

  const int first_theory=0;
  const int  last_theory=3;



  const int ciclo_lineal=100;
  const int tercio_ciclo_lineal=70;
  const double disminuye=0.869;
  const double aumenta=1.15;




  double lateral_param[3];	//	valores laterales
  double lateral_chi[3];
  int n_avanza;				//	pasos que desecha antes de colocar un nuevo eslabón





  double chi_actual; 	//	variables basicas
  double delta_param[Npca],param_old[Npca];	//	para evolucion lineal
  double sigma_param_abundancia[Npca];



void carga(void){
if(first_theory==0){
  param_abundancia[0]=1.13e-01;	sigma_param_abundancia[0]=1.e-02;	//	alpha:0
  param_abundancia[1]=1.25e-01; sigma_param_abundancia[1]=1.e-02;	//	alpha:1
  param_abundancia[2]=3.00e-08;	sigma_param_abundancia[2]=1.e-09;	//	alpha:log
  param_abundancia[3]=2.01e-01;	sigma_param_abundancia[3]=1.e-02;	//	A
  param_abundancia[4]=3.17e-00;	sigma_param_abundancia[4]=1.e-01;	//	beta
  param_abundancia[5]=1.13e-00;	sigma_param_abundancia[5]=1.e-01;	//	c: factor dentro de la exponencial
}
else if(first_theory==3){
  param_abundancia[0]=3.0e-01;	sigma_param_abundancia[0]=1.e-02;	//	alpha:0
  param_abundancia[1]=1.9e-01;	sigma_param_abundancia[1]=1.e-02;	//	alpha:1
  param_abundancia[2]=3.0e-00;	sigma_param_abundancia[2]=1.e-01;	//	alpha:exp
  param_abundancia[3]=2.21e-01;	sigma_param_abundancia[3]=1.e-02;	//	A
  param_abundancia[4]=3.07e-00;	sigma_param_abundancia[4]=1.e-01;	//	beta
  param_abundancia[5]=1.17e-00;	sigma_param_abundancia[5]=1.e-01;	//	c: factor dentro de la exponencial
}


if(segunda_vez==1){
  FILE * NOM;
  char nomAr[210];int n_salta;
  sprintf(nomAr,"list_abundancia_%d_%d.txt",first_theory,last_theory);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le\n"
   ,&n_salta,&co,&param_abundancia[0],&param_abundancia[1],
&param_abundancia[2],&param_abundancia[3],&param_abundancia[4],&param_abundancia[5])!=EOF);
  fclose(NOM);
}
printf("co=%le\n",co);


#if sigma_propto_param > 0
  for(int aux_param=0;aux_param<Npca;aux_param++){
    sigma_param_abundancia[aux_param]=0.1*abs(param_abundancia[aux_param]);}
#endif


}









void coeficientes(int id_theory){
double x;
  if(first_theory==0){
    x=log10(param_abundancia[2]+parametro[id_theory]);
    gama=param_abundancia[0]+param_abundancia[1]*x;}
  else if(first_theory==3){
    x=parametro[id_theory];
    gama=param_abundancia[0]+param_abundancia[1]*x/(1.+param_abundancia[2]*exp(-x*x));}
}








#include "../l/cosmo7_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_abundancia.h"//	subrutinas abundancia
#include "../l/cosmo7_densidad.h"	//	subrutinas densidad





//Programa principal

int main(){

if(first_theory==0)MG=1;
if(first_theory==3)MG=2;

int global=1;
//while(global>0){
global=0;
int no_tercio=1;

  Crandom ran(seed);
  int i,j,max_n_salta;

  FILE * NOM;
  FILE * SOM;
  FILE * OML;
  char nomAr[210];

  LeeAbundancia();	//	carga: abundancia,error_,radio_
  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100
  carga();		//	carga: parametros iniciales mcmc abundancia

  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(first_theory==0)param_cosmo[0]=log10(parametro[id_theory]);
    if(first_theory==3)param_cosmo[0]=parametro[id_theory];
    main_rodrigo(id_theory);	//	calcula P(k,z=0) y los sigma de perfiles, usa param_cosmo[0]
    calcula_sigma_derivadas(id_theory);
  }


  //	empieza calculo del chi2
  double chi2=0.0;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory!=without_id_theory){
      coeficientes(id_theory);
      calcula_abundancia_teorica(id_theory);	//	calcula la abundancia teorica, necesita gama
      for(i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      chi2+=pow((teoria_abundancia[id_theory][i]-abundancia[id_theory][i])
                 /error_abundancia[id_theory][i],2.0);
    }
  }
  lateral_chi[0]=chi2;
  co=chi2;
  cn=chi2;
printf("chi2=%le\n",chi2);

  for(int aux_param=0;aux_param<Npca;aux_param++){
    param_old[aux_param]=param_abundancia[aux_param];}










  sprintf(nomAr,"list_abundancia_%d_%d.txt",first_theory,last_theory);
  if(segunda_vez==1)	NOM=fopen(nomAr,"a+");
  else			NOM=fopen(nomAr,"w+");

  cout<<"comienza la cadena "<<endl;


n_avanza=1;
double delta,delta_min;
int lineal;
int id_param;

while(n_avanza>0){
  n_avanza=0;

int primer_lineal_global=1;
int lineal_individual=0;
int lin_param=0;
int hizo_almenos_uno=0;


  for(lineal=0;lineal<ciclo_lineal;lineal++){

  for(id_param=0;id_param<Npca;id_param++){





  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	ciclo lineal	////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  if((lineal==ciclo_lineal-1)&&(id_param==Npca-1)){

  //	calcula delta lineal
  if(primer_lineal_global==1){		// defino delta lineal
    fprintf(NOM,"%le ",param_abundancia[id_param]);
    primer_lineal_global=0;
    for(int aux_param=0;aux_param<Npca;aux_param++){
        delta_param[aux_param]=0.1*(param_abundancia[aux_param]-param_old[aux_param]);
  //      printf("%le %le\n",param_abundancia[aux_param],delta_param[aux_param]);
}}

  //	doy un paso lineal
  if(lineal_individual==1)	//	da un paso lineal individual
    param_abundancia[lin_param]+=delta_param[lin_param];
  else				//	da un paso lineal global
    for(int aux_param=0;aux_param<Npca;aux_param++)
        param_abundancia[aux_param]+=delta_param[aux_param];

  chi2=0.0;  
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory!=without_id_theory){
      coeficientes(id_theory);
      calcula_abundancia_teorica(id_theory);	//	calcula la abundancia teorica, necesita gama
      for(i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      chi2+=pow((teoria_abundancia[id_theory][i]-abundancia[id_theory][i])
                 /error_abundancia[id_theory][i],2.0);
    }
  }
  co = cn;
  cn = chi2;




    if(lineal_individual==1){		//	lineal individual
      id_param=Npca-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        hizo_almenos_uno++;		//	me aseguro hacer otro ciclo lineal individual
//printf(";");fflush(stdout);
      }
      else{
//printf(",");fflush(stdout);
	//	retrocede
        cn=co;
        param_abundancia[lin_param]-=delta_param[lin_param];
        lin_param++;
        if(lin_param>=Npca)
          lin_param-=Npca;

        if((lin_param==0)){	//	cliclo de intentos individuales completo
          if(hizo_almenos_uno>0){	//	intento otro global
            printf("%d ",hizo_almenos_uno);fflush(stdout);
            lineal_individual=0;
           }
          else{		//	sale
            id_param=Npca;			//	salgo de los ciclos lineales
            primer_lineal_global=1;
            lineal_individual=0;	//	reset
            lin_param=0;		//	reset
            hizo_almenos_uno=0;		//	reset
            printf(".\n");fflush(stdout);
          }
        }
      }
    }
    else{			//	lineal global
      id_param=Npca-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        printf("G");fflush(stdout);}		//	continuo, no cambio nada
      else{
        printf("g");fflush(stdout);
  //      printf(" %le %le ",co,cn);fflush(stdout);
	//	retrocede
        cn=co;
        for(int aux_param=0;aux_param<Npca;aux_param++)
            param_abundancia[aux_param]-=delta_param[aux_param];
        lineal_individual=1;
        hizo_almenos_uno=0;
        lin_param=0;
      }
    }



}
  //////////////////////////////////////////////////////////////////////////
  ////////////////////////////	ciclo normal(lateral)
  //////////////////////////////////////////////////////////////////////////
  else{
  lateral_param[0]=param_abundancia[id_param];
  lateral_chi[0]=cn;/////////////////////////////////////////////////////
  delta=sigma_param_abundancia[id_param]*0.02;
  delta_min=0.0001*delta;
  lateral_param[1]=lateral_param[0]-delta;
  lateral_param[2]=lateral_param[0]+delta;

  while(delta>delta_min){// numero de puntos a calcular en el espacio de fase

  for(int vecino=1;vecino<=2;vecino++){
  param_abundancia[id_param]=lateral_param[vecino];

  chi2=0.;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory!=without_id_theory){
      coeficientes(id_theory);
      calcula_abundancia_teorica(id_theory);	//	calcula la abundancia teorica, necesita gama
      for(i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      chi2+=pow((teoria_abundancia[id_theory][i]-abundancia[id_theory][i])
                 /error_abundancia[id_theory][i],2.0);
    }
  }
  lateral_chi[vecino]=chi2;


  }  //	for vecinos

//printf("%le %le %le\n",lateral_chi[0],lateral_chi[1],lateral_chi[2]);

int aux_avanza=0;
for(int vecino=1;vecino<=2;vecino++)
  if(lateral_chi[0]>lateral_chi[vecino]){
    lateral_chi[0]=lateral_chi[vecino];
    lateral_param[0]=lateral_param[vecino];
    aux_avanza++;
  }
if(aux_avanza==0)
  delta*=0.5;
else
  n_avanza++;

  lateral_param[1]=lateral_param[0]-delta;
  lateral_param[2]=lateral_param[0]+delta;
//printf("%d %d %d %le %le %le %le %le %le\n",id_param,id_coeff,n_avanza,lateral_chi[0],lateral_chi[1],lateral_chi[2], param[id_param][id_coeff],delta,delta_min);
}  //	while delta

  param_abundancia[id_param]=lateral_param[0];
  cn=lateral_chi[0];

  if((id_param==0)){
    fflush(NOM);
    fprintf(NOM,"\n%d %le ",n_avanza,lateral_chi[0]);}
  fprintf(NOM,"%le ",param_abundancia[id_param]);



}  //	else:ciclo normal (lateral)



  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	a la mitad del ciclo lineal defino param_old
  //////////////////////////////////////////////////////////////////////////
  if((lineal==tercio_ciclo_lineal-1)&&(id_param==Npca-1)){
  if(no_tercio==0){
    printf("\ntercio lineal\n");fflush(stdout);
  for(int aux_param=0;aux_param<Npca;aux_param++){
    param_old[aux_param]=param_abundancia[aux_param];}}
  else
  no_tercio=0;}




}  //	for id_param

}  //	for lineal
global+=n_avanza;
printf("na %d\n",n_avanza);fflush(stdout);
}  //	while hace algun paso


printf("n_avanza=%d lineal=%d id_param=%d\n",n_avanza,lineal,id_param);
printf("delta=%le delta_min=%le\n",delta,delta_min);

fprintf(NOM,"\n");
fclose(NOM);









//}




cout<<"terminó"<<endl;
  return 0;
}







