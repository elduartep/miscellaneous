
#include <iostream>
#include <fstream>
#include <cmath>
#include "../l/RANDOM.h"
//#include "param.h"
//NR3
#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "../l/cosmo7_parametros.h"

  using namespace std;


  const	int esferas=1;
  const int solo_analiza=0;
  const int segunda_vez=0;

  const int first_theory=3;
  const int  last_theory=6;

  #define sigma_propto_param 1


  char NameCase[]={'f','s'};

  const int xyerrorbars=0;


  const int BI=0;				//	numero de puntos de un parametro
  const int ciclo_lineal=400;
  const int tercio_ciclo_lineal=350;
  const double disminuye=0.869;
  const double aumenta=1.15;



  double lateral_param_densidad[3];	//	valores laterales
  double lateral_chi[3];
  int n_avanza;				//	pasos que desecha antes de colocar un nuevo eslabón




  double sigma_param_densidad[Npd][Ncd],chi_actual; 	//	variables basicas
  double delta_param_densidad[Npd][Ncd],param_old[Npd][Ncd];	//	para evolucion lineal


  double sigma_inv[theories],sigma_inv_max[theories],sigma_inv_min[theories];
  int cont_sigma[theories];
  double maximo,minimo;		//	limites globales de la interplación de sigma^-1



void carga(void){
if(first_theory==0){
  param_densidad[0][0]=13.9235;	// alpha
  param_densidad[0][1]=0.36234;
  param_densidad[0][2]=1.e-8;
  param_densidad[0][3]=-3.0869;
  param_densidad[0][4]=1.20586;	// beta

  param_densidad[1][0]=0.09273;	// rv
  param_densidad[1][1]=0.25986;
  param_densidad[1][2]=1.634;

  param_densidad[2][0]=-0.9631;	// G
  param_densidad[2][1]=3.4017e-02;
  param_densidad[2][2]=1.e-8;
  param_densidad[2][3]=1.60655;
  param_densidad[2][4]=-0.25995;
}
if(first_theory==3){
  param_densidad[0][0]=10.71;	// alpha
  param_densidad[0][1]=0.685;
  param_densidad[0][2]=0.1;
  param_densidad[0][3]=-2.9;
  param_densidad[0][4]=1.205;	// beta

  param_densidad[1][0]=0.054;	// rv
  param_densidad[1][1]=0.319;
  param_densidad[1][2]=1.375;

  param_densidad[2][0]=-1.248;	// G
  param_densidad[2][1]=8.23e-02;
  param_densidad[2][2]=-5.12e-03;
  param_densidad[2][3]=1.6;
  param_densidad[2][4]=-0.255;
}



int n_salta;
if(segunda_vez==1){
  FILE * NOM;
  char nomAr[410];
  sprintf(nomAr,"lista_todo_libre_%d_%d.txt",first_theory,last_theory);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n"
,&n_salta,&co,&param_densidad[0][0],&param_densidad[0][1],&param_densidad[0][2],&param_densidad[0][3],&param_densidad[0][4],&param_densidad[1][0],&param_densidad[1][1],&param_densidad[1][2],&param_densidad[1][3],&param_densidad[1][4],&param_densidad[2][0],&param_densidad[2][1],&param_densidad[2][2],&param_densidad[2][3],&param_densidad[2][4])!=EOF);
  fclose(NOM);
  lateral_chi[0]=co;
}
printf("co=%le\n",co);



#if sigma_propto_param > 0
  for(int aux_param=0;aux_param<Npd;aux_param++){
  for(int aux_coeff=0;aux_coeff<Ncd;aux_coeff++){
    sigma_param_densidad[aux_param][aux_coeff]=0.1*abs(param_densidad[aux_param][aux_coeff]);}}
#endif
}











void ParametrosPerfilDesc(int t,int s){
double x,x2;
  if(first_theory==0){
      x=log10(param_densidad[0][2]+parametro[t]);
  alpha=param_densidad[0][0] +param_densidad[0][1]*x
        +param_densidad[0][3]*sigma_stack[t][s];
   beta=param_densidad[0][4]*alpha;
      B=param_densidad[1][0]+param_densidad[1][1]*pow(sigma_stack[t][s],-param_densidad[1][2]);
        A=(param_densidad[2][0] +param_densidad[2][1]*x
         +param_densidad[2][3]*sigma_stack[t][s] 
         +param_densidad[2][4]*pow(sigma_stack[t][s],2)) *pow(B,1./param_densidad[0][4]);}
  else if(first_theory==3){
      x=parametro[t];
     x2=x*x;
  alpha=param_densidad[0][0] +param_densidad[0][1]*x +param_densidad[0][2]*x2
        +param_densidad[0][3]*sigma_stack[t][s];
   beta=param_densidad[0][4]*alpha;
      B=param_densidad[1][0]+param_densidad[1][1]*pow(sigma_stack[t][s],-param_densidad[1][2]);
      A=(param_densidad[2][0] +param_densidad[2][1]*x +param_densidad[2][2]*x2
         +param_densidad[2][3]*sigma_stack[t][s] 
         +param_densidad[2][4]*pow(sigma_stack[t][s],2)) *pow(B,1./param_densidad[0][4]);}
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
  char nomAr[410];

  LeePerfiles();
  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100
  carga();		//	carga: parametros iniciales mcmc abundancia

  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(first_theory==0)param_cosmo[0]=log10(parametro[id_theory]);
    if(first_theory==3)param_cosmo[0]=parametro[id_theory];
    si_densidad=1;
    main_rodrigo(id_theory);	//	calcula P(k,z=0) y los sigma de perfiles, usa param_cosmo[0]
  }

  //	empieza calculo del chi2
  double chi2=0.0;
  //	el numero de bines radiales con información ahora depende del valor de bin_stack
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++)
    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
      ParametrosPerfilDesc(id_theory,id_stack);
      for(i=0;i<bines_radio_densidad[id_stack];i++)
        chi2+=pow((perfil(radio_densidad[i])-densidad[id_theory][id_stack][i])
               /error_densidad[id_theory][id_stack][i],2.0);
    }
  lateral_chi[0]=chi2;
  co=chi2;
  cn=chi2;
printf("chi2=%le\n",chi2);

  for(int aux_param=0;aux_param<Npd;aux_param++){
  for(int aux_coeff=0;aux_coeff<Ncd;aux_coeff++){
    param_old[aux_param][aux_coeff]=param_densidad[aux_param][aux_coeff];}}










  sprintf(nomAr,"lista_todo_libre_%d_%d.txt",first_theory,last_theory);
  if(segunda_vez==1)	NOM=fopen(nomAr,"a+");
  else			NOM=fopen(nomAr,"w+");

  cout<<"comienza la cadena "<<endl;


n_avanza=1;
double delta,delta_min;
int lineal;
int id_param;
int id_coeff;

while(n_avanza>0){
  n_avanza=0;

int primer_lineal_global=1;
int lineal_individual=0;
int lin_param=0;
int lin_coeff=0;
int hizo_almenos_uno=0;


  for(lineal=0;lineal<ciclo_lineal;lineal++){

  for(id_param=0;id_param<Npd;id_param++){
  for(id_coeff=0;id_coeff<Ncd;id_coeff++){





  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	ciclo lineal	////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  if((lineal==ciclo_lineal-1)&&(id_param==Npd-1)&&(id_coeff==Ncd-1)){

  //	calcula delta lineal
  if(primer_lineal_global==1){		// defino delta lineal
    fprintf(NOM,"%le ",param_densidad[id_param][id_coeff]);
    primer_lineal_global=0;
    for(int aux_param=0;aux_param<Npd;aux_param++)
      for(int aux_coeff=0;aux_coeff<Ncd;aux_coeff++){
        delta_param_densidad[aux_param][aux_coeff]=0.1*(param_densidad[aux_param][aux_coeff]-param_old[aux_param][aux_coeff]);
        if(delta_param_densidad[aux_param][aux_coeff]==0.)
          delta_param_densidad[aux_param][aux_coeff]=param_densidad[aux_param][aux_coeff]*0.0001;
//        printf("%le %le\n",param_densidad[aux_param][aux_coeff],delta_param_densidad[aux_param][aux_coeff]);
}}

  //	doy un paso lineal
  if(lineal_individual==1)	//	da un paso lineal individual
    param_densidad[lin_param][lin_coeff]+=delta_param_densidad[lin_param][lin_coeff];
  else				//	da un paso lineal global
    for(int aux_param=0;aux_param<Npd;aux_param++)
      for(int aux_coeff=0;aux_coeff<Ncd;aux_coeff++)
        param_densidad[aux_param][aux_coeff]+=delta_param_densidad[aux_param][aux_coeff];

  //	empieza calculo del chi2
  double chi2=0.0;
  //	el numero de bines radiales con información ahora depende del valor de bin_stack
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++)
    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
      ParametrosPerfilDesc(id_theory,id_stack);
      for(i=0;i<bines_radio_densidad[id_stack];i++)
        chi2+=pow((perfil(radio_densidad[i])-densidad[id_theory][id_stack][i])
               /error_densidad[id_theory][id_stack][i],2.0);
    }
  co = cn;
  cn = chi2;




    if(lineal_individual==1){		//	lineal individual
      id_coeff=Ncd-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        hizo_almenos_uno++;		//	me aseguro hacer otro ciclo lineal individual
//printf(";");fflush(stdout);
      }
      else{
//printf(",");fflush(stdout);
	//	retrocede
        cn=co;
        param_densidad[lin_param][lin_coeff]-=delta_param_densidad[lin_param][lin_coeff];
        lin_coeff++;
        if(lin_coeff>=Ncd){
          lin_coeff=0;
          lin_param++;}
        if(lin_param>=Npd)
          lin_param=0;

        if((lin_param==0)&&(lin_coeff==0)){	//	cliclo de intentos individuales completo
          if(hizo_almenos_uno>0){	//	intento otro global
            printf("%d ",hizo_almenos_uno);fflush(stdout);
            lineal_individual=0;
           }
          else{		//	sale
            id_coeff=Ncd;			//	salgo de los ciclos lineales
            primer_lineal_global=1;
            lineal_individual=0;	//	reset
            lin_param=0;		//	reset
            lin_coeff=0;		//	reset
            hizo_almenos_uno=0;		//	reset
            //	si en todo el ciclo anterior no azanzó nada
 /*           int todos=0;
            for(int aux_param=0;aux_param<M;aux_param++)
              for(int aux_coeff=0;aux_coeff<N;aux_coeff++)
                if(delta_param_densidad[aux_param][aux_coeff]!=0.)
                  todos=1;
            if(todos!=0)
            n_avanza++;*/
            printf(".\n");fflush(stdout);
          }
        }
      }
    }
    else{			//	lineal global
      id_coeff=Ncd-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        printf("G");fflush(stdout);}		//	continuo, no cambio nada
      else{
        printf("g");fflush(stdout);
  //      printf(" %le %le ",co,cn);fflush(stdout);
	//	retrocede
        cn=co;
        for(int aux_param=0;aux_param<Npd;aux_param++)
          for(int aux_coeff=0;aux_coeff<Ncd;aux_coeff++)
            param_densidad[aux_param][aux_coeff]-=delta_param_densidad[aux_param][aux_coeff];
        lineal_individual=1;
        hizo_almenos_uno=0;
        lin_param=0;
        lin_coeff=0;
      }
    }


}
  //////////////////////////////////////////////////////////////////////////
  ////////////////////////////	ciclo normal(lateral)
  //////////////////////////////////////////////////////////////////////////
  else{
  lateral_param_densidad[0]=param_densidad[id_param][id_coeff];
  lateral_chi[0]=cn;
  delta=sigma_param_densidad[id_param][id_coeff]*0.02;
  delta_min=0.0001*delta;
  lateral_param_densidad[1]=lateral_param_densidad[0]-delta;
  lateral_param_densidad[2]=lateral_param_densidad[0]+delta;

  while(delta>delta_min){// numero de puntos a calcular en el espacio de fase

  for(int vecino=1;vecino<=2;vecino++){
  param_densidad[id_param][id_coeff]=lateral_param_densidad[vecino];


  //	empieza calculo del chi2
  double chi2=0.0;
  //	el numero de bines radiales con información ahora depende del valor de bin_stack
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++)
    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
      ParametrosPerfilDesc(id_theory,id_stack);
      for(i=0;i<bines_radio_densidad[id_stack];i++)
        chi2+=pow((perfil(radio_densidad[i])-densidad[id_theory][id_stack][i])
               /error_densidad[id_theory][id_stack][i],2.0);
    }
  lateral_chi[vecino]=chi2;


  }  //	for vecinos

//printf("%le %le %le\n",lateral_chi[0],lateral_chi[1],lateral_chi[2]);

int aux_avanza=0;
for(int vecino=1;vecino<=2;vecino++)
  if(lateral_chi[0]>lateral_chi[vecino]){
    lateral_chi[0]=lateral_chi[vecino];
    lateral_param_densidad[0]=lateral_param_densidad[vecino];
    aux_avanza++;
  }
if(aux_avanza==0)
  delta*=0.5;
else
  n_avanza++;

  lateral_param_densidad[1]=lateral_param_densidad[0]-delta;
  lateral_param_densidad[2]=lateral_param_densidad[0]+delta;
}  //	while delta

  param_densidad[id_param][id_coeff]=lateral_param_densidad[0];
  cn=lateral_chi[0];

  if((id_param==0)&&(id_coeff==0)){
    fflush(NOM);
    fprintf(NOM,"\n%d %le ",n_avanza,lateral_chi[0]);}
  fprintf(NOM,"%le ",param_densidad[id_param][id_coeff]);


}  //	else:ciclo normal (lateral)



  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	a la mitad del ciclo lineal defino param_old
  //////////////////////////////////////////////////////////////////////////
  if((lineal==tercio_ciclo_lineal-1)&&(id_param==Npd-1)&&(id_coeff==Ncd-1)){
  if(no_tercio==0){
    printf("\ntercio lineal\n");fflush(stdout);
  for(int aux_param=0;aux_param<Npd;aux_param++){
  for(int aux_coeff=0;aux_coeff<Ncd;aux_coeff++){
    param_old[aux_param][aux_coeff]=param_densidad[aux_param][aux_coeff];
#if sigma_propto_param > 0
    sigma_param_densidad[aux_param][aux_coeff]=0.1*abs(param_densidad[aux_param][aux_coeff]);
#endif
}}}
  else
  no_tercio=0;}



}  //	for id_coeff

}  //	for id_param

}  //	for lineal
global+=n_avanza;
}  //	while hace algun paso


printf("n_avanza=%d lineal=%d id_param=%d id_coeff=%d\n",n_avanza,lineal,id_param,id_coeff);
printf("delta=%le delta_min=%le\n",delta,delta_min);

fprintf(NOM,"\n");
fclose(NOM);









//}












cout<<"terminó"<<endl;
  return 0;
}







