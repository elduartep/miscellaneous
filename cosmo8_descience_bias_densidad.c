#include <iostream>
#include <fstream>
#include <cmath>

#include "../l/RANDOM.h"

#include "../l/nr3.h"
#include "../l/interp_1d.h"

#include "cosmo8_parametros.h"

  using namespace std;


int  yes_bias=0;
int  yes_densidad=1;

  const	int esferas=1;

  const int segunda_vez=0;


  #define sigma_propto_param 1

  #define bias_teoria 1
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


if(si_bias==1){
  if(MG==1){
    param_bd[0]=-0.284;	// c
    param_bd[1]=3.365e-3;	// c
    param_bd[2]=0.594;	// c
    param_bd[3]=-0.211153;}	// c
  else{
    param_bd[0]=-0.347;	// c
    param_bd[1]=2.38e-2;	// c
    param_bd[3]=0.607;	// c
    param_bd[2]=-0.2198;}	// c
}

if(si_densidad==1){
if(first_theory==0){
  param_bd[0]=1.227669e+01;	// alpha
  param_bd[1]=3.168121e-01;
  param_bd[2]=1.596831e-08;
  param_bd[3]=-2.636782e+00;
  param_bd[4]=1.240896e+00;	// beta

  param_bd[5]=0.;	// rv
  param_bd[6]=3.959814e-01;
  param_bd[7]=1.347607e+00;

  param_bd[8]=0.4;//-1.105609e+01;	// G
  param_bd[9]=3.612350e-02;
  param_bd[10]=0.;
  param_bd[11]=3.;//1.148487e+01;
  param_bd[12]=0.;//1.079775e-01;
}
if(first_theory==3){
  param_bd[0]=9.790662e+00;	// alpha
  param_bd[1]=4.664299e-01;
  param_bd[2]=1.;
  param_bd[3]=-2.521648e+00;
  param_bd[4]=1.228912e+00;	// beta

  param_bd[5]=0.;	// rv
  param_bd[6]=4.212162e-01;
  param_bd[7]=1.389330e+00;

  param_bd[8]=0.16;//-1.530806e+01;	// G
  param_bd[9]=8.624882e-02;
  param_bd[10]=-7.881742e-03;
  param_bd[11]=2.9;//1.543693e+01;
  param_bd[12]=0.;//8.092481e-02;
}
}


int n_salta;
if(segunda_vez==1){
  FILE * NOM;
  char nomAr[410];
  sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",first_theory,last_theory,si_bias,si_densidad);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le\n"
,&n_salta,&co,&param_bd[0],&param_bd[1],&param_bd[2],&param_bd[3],&param_bd[4],&param_bd[5],&param_bd[6],&param_bd[7],&param_bd[8],&param_bd[9],&param_bd[10])!=EOF);
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

//      CargaSigma(id_theory);//	pone en sig_actual el sigma correpondiente a id_theory
//      Spline_interp sigm(radio_xi,sig_actual);
//      alpha1[id_theory]=param_bd[Npd+0]*sigm.interp(param_bd[Npd+1]);

    if(si_bias==1){
      for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
        bias_aux=bias(sigma_bias[id_theory][id_stack],id_theory);
        chi2+=pow((bias_aux - medida_bias[id_theory][id_stack])/error_bias[id_theory][id_stack],2.0);
      }
    }

    if(si_densidad==1){
//      CargaCorrelacion(id_theory);//	pone en xi_actual la correlacion correpondiente a id_theory
//      Spline_interp correl(radio_xi,xi_actual);

//      CargaSigma(id_theory);//	pone en sig_actual el sigma correpondiente a id_theory
//      Spline_interp sigm(radio_xi,sig_actual);
//      alpha1[id_theory]=param_bd[Npd+0]*sigm.interp(param_bd[Npd+1]);


    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
      ParametrosPerfil(id_theory,id_stack);
      for(int i=0;i<bines_radio_densidad[id_stack];i++)
        if((i!=19)&&(i!=20))
        chi2+=pow( (perfil(radio_densidad[i])-densidad[id_theory][id_stack][i]) /error_densidad[id_theory][id_stack][i],2.0);
    }
    }

  }
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
  LeeBias2();		//	read data
  carga();		//	carga: parametros iniciales mcmc bias+densidad
  carga_espectro_99();	//	carga: espectro lcdm en 99 y 100

  calcula_correlation_fuction=1;	//calcula la funcion de correlacion cuando entre a main_rodrigo
  for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
    if(first_theory<3)	{param_cosmo[0]=log10(parametro[id_theory]);	MG=1;}
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
//        printf("%le %le\n",param_bd[aux_param][aux_coeff],delta_param_bd[aux_param][aux_coeff]);
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
                if(delta_param_bd[aux_param][aux_coeff]!=0.)
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
for(int id_theory=first_theory;id_theory<=last_theory;id_theory++){
  sprintf(nomAr,"../%c/params_best-fit_bd_%d_%d.dat",NameTheories[id_theory],first_theory,last_theory);
  PLO=fopen(nomAr,"w+");
  fprintf(PLO,"b0=%le\n",b1);
  fprintf(PLO,"b1=%le\n",b2);
  fprintf(PLO,"b2=%le\n",b3);
  for(int id_stack=first_stack;id_stack<=last_stack;id_stack++){
    ParametrosPerfil(id_theory,id_stack);
    fprintf(PLO,"a%d=%le\n",id_stack,alpha);
    fprintf(PLO,"b%d=%le\n",id_stack,beta);
    fprintf(PLO,"A%d=%le\n",id_stack,A);
    fprintf(PLO,"B%d=%le\n",id_stack,B);
  }
  fclose(PLO);
}






cout<<"terminó"<<endl;
  return 0;
}







