#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "cosmo8b_parametros.h"

  using namespace std;


int  yes_bias=1;
int  yes_densidad=1;

  const	int esferas=1;

  const int segunda_vez=0;


  #define sigma_propto_param 1

  #define densidad_teoria 1


  const int ciclo_lineal=100;
  const int tercio_ciclo_lineal=70;
  const double disminuye=0.869;
  const double aumenta=1.15;


  double lateral_param_bd[3];	//	valores laterales
  double lateral_chi[3];
  int n_avanza;				//	pasos que desecha antes de colocar un nuevo eslabón


  double sigma_param_bd[Npbd],chi_actual; 	//	variables basicas
  double delta_param_bd[Npbd],param_old[Npbd];	//	para evolucion lineal



//////////////////////////////////////////////////////
void carga(void){
//////////////////////////////////////////////////////


  param_bd[0]=8.965;	// alpha	3
  param_bd[1]=10.583;	// beta		4
  param_bd[2]=1.38;	// c		5
  param_bd[3]=1.0802;	// G		6

for(int s=first_stack;s<=last_stack;s++){
for(int t=first_theory;t<=last_theory;t++){
param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+0]=5.;	//A
param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+1]=radio_stack[t][s]*1.25;	//B
param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+2]=medida_bias[t][s];	//b	0.6-0.013*pow(radio_stack[t][s],2)
}}

  int i;
  for(i=0;i<4;i++)
    printf("Pbd%d=%le %d\n",i,param_bd[i],i);
  for(int s=first_stack;s<=last_stack;s++){
    for(int t=first_theory;t<=last_theory;t++){
      printf("A%c%d=%le %d\n",NameTheories[t],s,param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+0],i);i++;
      printf("B%c%d=%le %d\n",NameTheories[t],s,param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+1],i);i++;
      printf("b%c%d=%le %d\n",NameTheories[t],s,param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+2],i);i++;
    }
  }


int n_salta;
if(segunda_vez==1){
  FILE * NOM;
  char nomAr[410];
  sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",first_theory,last_theory,si_bias,si_densidad);
  NOM=fopen(nomAr,"r");
  fscanf(NOM,"%d %le\n",&n_salta,&co);
  for(int i=0;i<Npbd;i++)
    fscanf(NOM,"%le ",&param_bd[i]);
  fclose(NOM);
  lateral_chi[0]=co;
}
printf("co=%le\n",co);



#if sigma_propto_param > 0
  for(int aux_param=0;aux_param<Npbd;aux_param++){
    sigma_param_bd[aux_param]=0.1*abs(param_bd[aux_param]);}
#endif
}





















#include "../l/cosmo8_MG.h"//	corre camb hasta z=100, integra con MG hasta z=0, calcula sigma y dlns/dlnR

#include "../l/cosmo7_abundancia.h"//	subrutinas abundancia
#include "../l/cosmo8_densidad_bias.h"	//	subrutinas densidad




//////////////////////////////////////////////////////
double Chi2_bd(void){
//////////////////////////////////////////////////////
double chi2=0.;
  double bias_aux,densidad_aux;
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){

    if(si_bias==1){
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        bias_aux=param_bd[4+3*(last_theory-first_theory+1)*(id_stack-first_stack)+3*(id_theory-first_theory)+2];
        chi2+=pow((bias_aux - medida_bias[id_theory][id_stack])/error_bias[id_theory][id_stack],2.0);
      }
    }

    if(si_densidad==1){
      CargaCorrelacion(id_theory);//	pone en xi_actual la correlacion correpondiente a id_theory
      Spline_interp correl(radio_xi,xi_actual);
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        bias_aux = param_bd[4+3*(last_theory-first_theory+1)*(id_stack-first_stack)+3*(id_theory-first_theory)+2];
        ParametrosSuppress(id_theory,id_stack);
        for(int i=0;i<bines_radio_densidad[id_stack];i++){
          if(densidad_teoria==1)	densidad_aux = PerfilUniversal(radio_densidad[i]);
          else				densidad_aux = (densidad[1][5][i]+densidad[5][5][i])*0.5;
          if((i!=19)&&(i!=20))
            chi2+=pow( (densidad_aux + bias_aux * correl.interp(radio_densidad[i]*radio_stack[id_theory][id_stack]) *  Suppress(radio_densidad[i]) - densidad[id_theory][id_stack][i]) / error_densidad[id_theory][id_stack][i],2.0);
        }
      }
    }

  }

  for(int s=first_stack;s<=last_stack;s++){
    for(int t=first_theory;t<=last_theory;t++){
      if(param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+0]>20.)  chi2*=20.;
      if(param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+1]>30.)  chi2*=30.;
  }}

if(param_bd[0]>param_bd[1])
  chi2*=20.;

return chi2;
}












//Programa principal
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
int main(){
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

si_bias=yes_bias;
si_densidad=yes_densidad;


if(si_bias==0)
  CargaAjustesBias();

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

  LeePerfiles();	//	read data
  LeeBias();		//	read data
//  LeeBias2();		//	read data
  carga();		//	carga: parametros iniciales mcmc bias+densidad
  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100

  calcula_correlation_fuction=1;	//calcula la funcion de correlacion cuando entre a main_rodrigo
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(id_theory<3)	{param_cosmo[0]=log10(parametro[id_theory]);	MG=1;}
    else		{param_cosmo[0]=parametro[id_theory];		MG=2;}
    main_rodrigo(id_theory);	//	calcula P(k,z=0), los sigma de perfiles y la funcion de correlacion materia, usa param_cosmo[0]
  }

  //	empieza calculo del chi2
  double chi2=0.0;
  //	el numero de bines radiales con información ahora depende del valor de bin_stack
  chi2=Chi2_bd();
  lateral_chi[0]=chi2;
  co=chi2;
  cn=chi2;
printf("chi2=%le\n",chi2);

  for(int aux_param=0;aux_param<Npbd;aux_param++){
    param_old[aux_param]=param_bd[aux_param];}









  sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",first_theory,last_theory,si_bias,si_densidad);
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

if ((10*(lineal+1))%ciclo_lineal==0){
  printf("%d%% ",100*(lineal+1)/ciclo_lineal);fflush(stdout);}
//printf("lineal=%d\n",lineal);

  for(id_param=0;id_param<Npbd;id_param++){

//printf("id_param=%d\n",id_param);



  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	ciclo lineal	////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  if((lineal==ciclo_lineal-1)&&(id_param==Npbd-1)){
//printf("ciclo lineal\n");
  //	calcula delta lineal
  if(primer_lineal_global==1){		// defino delta lineal
    printf("\n");
    fprintf(NOM,"%le ",param_bd[id_param]);
    primer_lineal_global=0;
    for(int aux_param=0;aux_param<Npbd;aux_param++){
        delta_param_bd[aux_param]=0.1*(param_bd[aux_param]-param_old[aux_param]);
        if(delta_param_bd[aux_param]==0.)
          delta_param_bd[aux_param]=param_bd[aux_param]*0.0001;
//        printf("%le %le\n",param_densidad[aux_param][aux_coeff],delta_param_densidad[aux_param][aux_coeff]);
}}

  //	doy un paso lineal
  if(lineal_individual==1)	//	da un paso lineal individual
    param_bd[lin_param]+=delta_param_bd[lin_param];
  else				//	da un paso lineal global
    for(int aux_param=0;aux_param<Npbd;aux_param++)
        param_bd[aux_param]+=delta_param_bd[aux_param];

  //	empieza calculo del chi2
  chi2=Chi2_bd();
  co = cn;
  cn = chi2;
//printf("chi2=%le\n",chi2);fflush(stdout);



    if(lineal_individual==1){		//	lineal individual
      id_param=Npbd-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        hizo_almenos_uno++;		//	me aseguro hacer otro ciclo lineal individual
//printf(";");fflush(stdout);
      }
      else{
//printf(",");fflush(stdout);
	//	retrocede
        cn=co;
        param_bd[lin_param]-=delta_param_bd[lin_param];
        lin_param++;
        if(lin_param>=Npbd)
          lin_param=0;

        if(lin_param==0){	//	cliclo de intentos individuales completo
          if(hizo_almenos_uno>0){	//	intento otro global
            printf("%d ",hizo_almenos_uno);fflush(stdout);
            lineal_individual=0;
           }
          else{		//	sale
            id_param=Npbd;			//	salgo de los ciclos lineales
            primer_lineal_global=1;
            lineal_individual=0;	//	reset
            lin_param=0;		//	reset
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
      id_param=Npbd-2;	//	no salgo del ciclo lineal
      if(cn<co){
        n_avanza++;
        printf("G");fflush(stdout);}		//	continuo, no cambio nada
      else{
        printf("g");fflush(stdout);
  //      printf(" %le %le ",co,cn);fflush(stdout);
	//	retrocede
        cn=co;
        for(int aux_param=0;aux_param<Npbd;aux_param++)
            param_bd[aux_param]-=delta_param_bd[aux_param];
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
//printf("aaa id_param=%d\n",id_param);fflush(stdout);
  lateral_param_bd[0]=param_bd[id_param];
  lateral_chi[0]=cn;
  delta=sigma_param_bd[id_param]*0.02;
  delta_min=0.0001*delta;
  lateral_param_bd[1]=lateral_param_bd[0]-delta;
  lateral_param_bd[2]=lateral_param_bd[0]+delta;

  int control=0;
  while(delta>delta_min){// numero de puntos a calcular en el espacio de fase
  control++;
  if(control%1000==0){
    printf(" %dk->%d %le ",control/1000,id_param,lateral_param_bd[0]);fflush(stdout);}
  for(int vecino=1;vecino<=2;vecino++){
  param_bd[id_param]=lateral_param_bd[vecino];


  //	empieza calculo del chi2
  chi2=Chi2_bd();
  lateral_chi[vecino]=chi2;
//printf("chi2=%le\n",chi2);fflush(stdout);
//printf("/");fflush(stdout);

  }  //	for vecinos

//printf("%le %le %le\n",lateral_chi[0],lateral_chi[1],lateral_chi[2]);

int aux_avanza=0;
for(int vecino=1;vecino<=2;vecino++)
  if(lateral_chi[0]>lateral_chi[vecino]){
    lateral_chi[0]=lateral_chi[vecino];
    lateral_param_bd[0]=lateral_param_bd[vecino];
    aux_avanza++;
  }
if(aux_avanza==0)
  delta*=0.5;
else
  n_avanza++;

  lateral_param_bd[1]=lateral_param_bd[0]-delta;
  lateral_param_bd[2]=lateral_param_bd[0]+delta;
//printf("%d %d %le %le %le %le %le %le\n",id_param,n_avanza,lateral_chi[0],lateral_chi[1],lateral_chi[2], param_bd[id_param],delta,delta_min);
}  //	while delta

  param_bd[id_param]=lateral_param_bd[0];
  cn=lateral_chi[0];

  if(id_param==0){
    fflush(NOM);
    fprintf(NOM,"\n%d %le ",n_avanza,lateral_chi[0]);}
  fprintf(NOM,"%le ",param_bd[id_param]);


}  //	else:ciclo normal (lateral)



  //////////////////////////////////////////////////////////////////////////
  /////////////////////////////	a la mitad del ciclo lineal defino param_old
  //////////////////////////////////////////////////////////////////////////
  if((lineal==tercio_ciclo_lineal-1)&&(id_param==Npbd-1)){
  if(no_tercio==0){
    printf(" tercio lineal ");fflush(stdout);
  for(int aux_param=0;aux_param<Npbd;aux_param++){
    param_old[aux_param]=param_bd[aux_param];
#if sigma_propto_param > 0
    sigma_param_bd[aux_param]=0.1*abs(param_bd[aux_param]);
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
  FILE *PLO;
  sprintf(nomAr,"params_best-fit_bd_%d_%d.dat",first_theory,last_theory);
  PLO=fopen(nomAr,"w+");
  for(i=0;i<4;i++)
    fprintf(PLO,"Pbd%d=%le\n",i,param_bd[i]);
  for(int s=first_stack;s<=last_stack;s++){
    for(int t=first_theory;t<=last_theory;t++){
      fprintf(PLO,"A%c%d=%le\n",NameTheories[t],s,param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+0]);
      fprintf(PLO,"B%c%d=%le\n",NameTheories[t],s,param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+1]);
      fprintf(PLO,"b%c%d=%le\n",NameTheories[t],s,param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+2]);
    }
  }
  fclose(PLO);

  for(int t=first_theory;t<=last_theory;t++){
    sprintf(nomAr,"../%c/best-fit_bias_%d_%d.dat",NameTheories[t],first_theory,last_theory);
    PLO=fopen(nomAr,"w+");
    for(int s=first_stack;s<=last_stack;s++)
      fprintf(PLO,"%le %le\n",radio_stack[t][s],param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+2]);
    fclose(PLO);
  }






cout<<"terminó"<<endl;
  return 0;
}







