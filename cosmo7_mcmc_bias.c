//	necesita Npbb=21 en cosmo7_parametros.h

//	primero intente con todo libre
//	eso sugiere una relacion lineal entre los coeficientes b1,b2,b3
//	parece que es mejor dejar fijo a b1

//	cuando lo hice todavia se nota una relacion lineal entre b2 y b3

//	aplico esta relacion en cosmo7_desciende_bias.c

#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "../l/cosmo7_parametros.h"

  using namespace std;

  #define sigma_propto_param 1


  const int segunda_vez=0;

  char NameCase[]={'f','s'};

  const int first_theory=0;
  const int  last_theory=6;



  const int ciclo_lineal=100;
  const int tercio_ciclo_lineal=70;
  const double disminuye=0.869;
  const double aumenta=1.15;




  double lateral_param[3];	//	valores laterales
  double lateral_chi[3];
  int n_avanza;				//	pasos que desecha antes de colocar un nuevo eslabón





  double delta_param[Npbb],param_old[Npbb];	//	para evolucion lineal
  double sigma_param_bias[Npbb];


void carga(void){
  for (int t=0;t<7;t++){
    param_bias[t]=0.6;		// corte
    param_bias[t+7]=-0.4;	// 4
    param_bias[t+14]=-0.17;	// coeficiente potencia 4
}


if(segunda_vez==1){
  FILE * NOM;
  char nomAr[410];int n_salta;
  sprintf(nomAr,"lista_bias_%d_%d.txt",first_theory,last_theory);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le\n"
,&n_salta,&co,&param_bias[0],&param_bias[1],&param_bias[2],&param_bias[3],&param_bias[4]
,&param_bias[5],&param_bias[6],&param_bias[7],&param_bias[8])!=EOF);
  fclose(NOM);
  lateral_chi[0]=co;
}
printf("co=%le\n",co);



#if sigma_propto_param > 0
  for(int aux_param=0;aux_param<Npbb;aux_param++){
    sigma_param_bias[aux_param]=0.1*abs(param_bias[aux_param]);}
#endif
}








void coeficientes(int t){
double x,x2;
  if(first_theory==0)
      x=log10(parametro[t]);
  else if(first_theory==3)
      x=parametro[t];
  b1=param_bias[t];
  b2=param_bias[t+7];
  b3=param_bias[t+14];
}







#include "../l/cosmo7_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_bias.h"//	subrutinas abundancia












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
  char nomAr[410];

  LeeBias();
  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100
  carga();		//	carga: parametros iniciales mcmc abundancia

  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(first_theory==0)param_cosmo[0]=log10(parametro[id_theory]);
    if(first_theory==3)param_cosmo[0]=parametro[id_theory];
    si_bias=1;
    main_rodrigo(id_theory);	//	calcula P(k,z=0) y los sigma de perfiles, usa param_cosmo[0]
  }

  //	empieza calculo del chi2
  double chi2=0.0;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    coeficientes(id_theory);
    for(int id_bin=first_bin_bias;id_bin<=last_bin_bias;id_bin++){
      chi2+=pow((bias(sigma_bias[id_theory][id_bin])-medida_bias[id_theory][id_bin])
               /error_bias[id_theory][id_bin],2.0);
    }}
  lateral_chi[0]=chi2;
  co=chi2;
  cn=chi2;
printf("chi2=%le\n",chi2);

  for(int aux_param=0;aux_param<Npbb;aux_param++){
    param_old[aux_param]=param_bias[aux_param];}










  sprintf(nomAr,"lista_bias_%d_%d.txt",first_theory,last_theory);
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

  for(id_param=0;id_param<Npbb;id_param++){





  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	ciclo lineal	////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  if((lineal==ciclo_lineal-1)&&(id_param==Npbb-1)){

  //	calcula delta lineal
  if(primer_lineal_global==1){		// defino delta lineal
    fprintf(NOM,"%le ",param_bias[id_param]);
    primer_lineal_global=0;
    for(int aux_param=0;aux_param<Npbb;aux_param++){
        delta_param[aux_param]=0.1*(param_bias[aux_param]-param_old[aux_param]);
        if(delta_param[aux_param]==0.)
          delta_param[aux_param]=param_bias[aux_param]*0.0001;
//        printf("%le %le\n",param_bias[aux_param],delta_param_bias[aux_param]);
}}

  //	doy un paso lineal
  if(lineal_individual==1)	//	da un paso lineal individual
    param_bias[lin_param]+=delta_param[lin_param];
  else				//	da un paso lineal global
    for(int aux_param=0;aux_param<Npbb;aux_param++)
        param_bias[aux_param]+=delta_param[aux_param];

  //	empieza calculo del chi2
  chi2=0.0;
  //	el numero de bines radiales con información ahora depende del valor de bin_stack
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    coeficientes(id_theory);
    for(int id_bin=first_bin_bias;id_bin<=last_bin_bias;id_bin++){
        chi2+=pow((bias(sigma_bias[id_theory][id_bin])-medida_bias[id_theory][id_bin])
               /error_bias[id_theory][id_bin],2.0);
    }}
  co = cn;
  cn = chi2;




    if(lineal_individual==1){		//	lineal individual
      id_param=Npbb-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        hizo_almenos_uno++;		//	me aseguro hacer otro ciclo lineal individual
//printf(";");fflush(stdout);
      }
      else{
//printf(",");fflush(stdout);
	//	retrocede
        cn=co;
        param_bias[lin_param]-=delta_param[lin_param];
        lin_param++;
        if(lin_param>=Npbb)
          lin_param-=Npbb;

        if(lin_param==0){	//	cliclo de intentos individuales completo
          if(hizo_almenos_uno>0){	//	intento otro global
            printf("%d ",hizo_almenos_uno);fflush(stdout);
            lineal_individual=0;
           }
          else{		//	sale
            id_param=Npbb;			//	salgo de los ciclos lineales
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
      id_param=Npbb-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        printf("G");fflush(stdout);}		//	continuo, no cambio nada
      else{
        printf("g");fflush(stdout);
  //      printf(" %le %le ",co,cn);fflush(stdout);
	//	retrocede
        cn=co;
        for(int aux_param=0;aux_param<Npbb;aux_param++)
            param_bias[aux_param]-=delta_param[aux_param];
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
  lateral_param[0]=param_bias[id_param];
  lateral_chi[0]=cn;
  delta=sigma_param_bias[id_param]*0.02;
  delta_min=0.0001*delta;
  lateral_param[1]=lateral_param[0]-delta;
  lateral_param[2]=lateral_param[0]+delta;

  while(delta>delta_min){// numero de puntos a calcular en el espacio de fase

  for(int vecino=1;vecino<=2;vecino++){
  param_bias[id_param]=lateral_param[vecino];

  chi2=0.0;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    coeficientes(id_theory);
    for(int id_bin=first_bin_bias;id_bin<=last_bin_bias;id_bin++){
        chi2+=pow((bias(sigma_bias[id_theory][id_bin])-medida_bias[id_theory][id_bin])
               /error_bias[id_theory][id_bin],2.0);
    }}
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
}  //	while delta

  param_bias[id_param]=lateral_param[0];
  cn=lateral_chi[0];

  if(id_param==0){
    fflush(NOM);
    fprintf(NOM,"\n%d %le ",n_avanza,lateral_chi[0]);}
  fprintf(NOM,"%le ",param_bias[id_param]);


}  //	else:ciclo normal (lateral)



  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	a la mitad del ciclo lineal defino param_old
  //////////////////////////////////////////////////////////////////////////
  if((lineal==tercio_ciclo_lineal-1)&&(id_param==Npbb-1)){
  if(no_tercio==0){
    printf("\ntercio lineal\n");fflush(stdout);
  for(int aux_param=0;aux_param<Npbb;aux_param++){
    param_old[aux_param]=param_bias[aux_param];
#if sigma_propto_param > 0
    sigma_param_bias[aux_param]=0.1*abs(param_bias[aux_param]);
#endif
}}
  else
  no_tercio=0;}



}  //	for id_param

}  //	for lineal
global+=n_avanza;
}  //	while hace algun paso


printf("n_avanza=%d lineal=%d id_param=%d\n",n_avanza,lineal,id_param);
printf("delta=%le delta_min=%le\n",delta,delta_min);

fprintf(NOM,"\n");
fclose(NOM);


  //	imprime 
  for(int t =first_theory;t<=last_theory;t++)
    printf("%le %le %le %le\n",parametro[t],param_bias[t],param_bias[t+7],param_bias[t+14]);
  fflush(stdout);

  //	imprime
  FILE *PLO;
  sprintf(nomAr,"bias_%d_%d.dat",first_theory,last_theory);
  PLO=fopen(nomAr,"w+");
  for(int t =first_theory;t<=last_theory;t++){
    fprintf(PLO,"a%c=%le\n",NameTheories[t],param_bias[t]);
    fprintf(PLO,"b%c=%le\n",NameTheories[t],param_bias[t+7]);
    fprintf(PLO,"c%c=%le\n",NameTheories[t],param_bias[t+14]);}
  fclose(PLO);




//}











cout<<"terminó"<<endl;
  return 0;
}







