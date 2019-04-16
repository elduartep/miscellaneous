////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////	bias
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


//	lee el resultado de cosmo8_desciende_bias_densidad.c
//	incompleto
void CargaAjustesBiasDensidad(void){
  printf("leyendo ajuste bias-densidad...");fflush(stdout);
  FILE * NOM;
  int n_salta_aux;
  char nomAr[200];
  sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",first_theory,last_theory,0,1);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le\n",&n_salta_aux,&co,
&param_bd[0],&param_bd[1],&param_bd[2],&param_bd[3],&param_bd[4],&param_bd[5],&param_bd[6],&param_bd[7],&param_bd[8],&param_bd[9],&param_bd[10])!=EOF);
  fclose(NOM);
  printf("hecho\n");fflush(stdout);
}



//	cuando voy a ajustar la densidad debo usar el resultado del bias
void CargaAjustesBias(void){
  printf("leyendo ajuste bias...\n");fflush(stdout);
  FILE * NOM;
  double aux;
  int n_salta_aux;
  char nomAr[200];
  sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",first_theory,last_theory,1,0);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le\n",&n_salta_aux,&co,
&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&b0,&b1,&b2)!=EOF);
  fclose(NOM);
  printf("b0=%le b1=%le b2=%le\n",b0,b1,b2);
  printf("hecho\n");fflush(stdout);
}




inline double bias(double sigma){
if(si_bias==1)
  return param_bd[Npd+0] - param_bd[Npd+1]*pow(sigma,-param_bd[Npd+2]);
else
  return b0 - b1 * pow(sigma,-b2);
}






void LeeBias(void){
  int i,j;
  double aux;
  printf("leyendo bias...");fflush(stdout);
  char NomArch[300];
  FILE * Ud;

for(int id_theory=0;id_theory<=6;id_theory++){
  sprintf(NomArch,"../%c/bv_03.txt",NameTheories[id_theory]);
  Ud=fopen(NomArch,"r");
  for(j=0;j<bin_stack;j++)
    fscanf(Ud,"%lf %lf %le %le %le %le\n",&radio_bias[id_theory][j], &aux, &medida_bias[id_theory][j], &error_bias[id_theory][j],&aux,&aux);

  fclose(Ud);
}
  printf("hecho\n");fflush(stdout);
}



void LeeBias2(void){
  int i,j;
  double aux,error;
  printf("leyendo bias...");fflush(stdout);
  char NomArch[300];
  FILE * Ud;

for(int id_theory=0;id_theory<=6;id_theory++){
  sprintf(NomArch,"ajuste_BIAS_%c.txt",NameTheories[id_theory]);
  Ud=fopen(NomArch,"r");
  for(j=0;j<bin_stack;j++)
    fscanf(Ud,"%le %le %le %le\n",&radio_bias2[id_theory][j], &aux, &medida_bias2[id_theory][j], &error_bias2[id_theory][j]);
  fclose(Ud);
}

  printf("hecho\n");fflush(stdout);
/*
//	imprimiendo la media de las dos medidas
for(int id_theory=0;id_theory<=6;id_theory++){
  sprintf(NomArch,"BIAS_%c.txt",NameTheories[id_theory]);
  Ud=fopen(NomArch,"w+");
  for(j=0;j<bin_stack;j++){
    aux=(medida_bias[id_theory][j]+medida_bias2[id_theory][j])*0.5;
    error=pow(error_bias[id_theory][j],2)*0.25 + pow(error_bias2[id_theory][j],2)*0.25;
    error+=pow(medida_bias[id_theory][j]-aux,2)*0.5 + pow(medida_bias2[id_theory][j]-aux,2)*0.5;
    medida_bias[id_theory][j]=aux;
    error_bias[id_theory][j]=sqrt(error);
    fprintf(Ud,"%le %le %le\n",radio_bias[id_theory][j], medida_bias[id_theory][j], error_bias[id_theory][j]);
  }
  fclose(Ud);
}
*/

}




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////	densidad
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////




void ParametrosSuppress(int t, int s){
//ley de potencia
//A=param_bd[4]*pow(radio_stack[t][s],-param_bd[5]);
//B=param_bd[6]*pow(radio_stack[t][s],+param_bd[7]);

//recta y parabola
A=param_bd[4]+pow(radio_stack[t][s]-6.9,2)*param_bd[5];
B=param_bd[6]+radio_stack[t][s]*param_bd[7];

//cuando todo es libre, completamente libre, incluyendo el bias
//A=param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+0];
//B=param_bd[4+3*(last_theory-first_theory+1)*(s-first_stack)+3*(t-first_theory)+1];
}



double Suppress(double y){
return exp(-A*pow(y,-B));
}

//	calculados usando cosmo8b_desiende_bias_densidad.c
double alpha=7.593649e+00;
double beta=9.381333e+00;
double ca=1.415600e+00;
double cb=1.106479e+00;

double PerfilUniversal(double r){
//return -dc/(1.+pow(r/param_bd[0],param_bd[1]));

//fijo
//return -dc*(1.-pow(r/ca,alpha))/(1.+pow(r/cb,beta));

//normal
return -dc*(1.-pow(r/param_bd[2],param_bd[0]))/(1.+pow(r/param_bd[3],param_bd[1]));

//return -dc*0.5*( (1.-pow(r/2.3,2))/(1.+pow(r/1.09,20))
  //       + (1.-pow(r/param_bd[2],param_bd[0]))/(1.+pow(r/param_bd[3],param_bd[1])) );

//return -dc*(1.-pow(r/3,2))/(1.+pow(r/1.09,20));
//         + (1.-pow(r/param_bd[2],param_bd[0]))/(1.+pow(r/param_bd[3],param_bd[1])) );

}



void CargaCorrelacion(int id_theory){
  for(int i=0;i<max_id_xi;i++)
    xi_actual[i]=xi[id_theory][i];
}




//	tomado de		mcmc_todo_libre5nuevofast_2.c
// lee el archivo que contiene los perfiles
void LeePerfiles(void){
  int i,j;
  double aux;
  printf("Leyendo perfiles...\n");
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
  fclose(Ud);
}


#if uniformiza_densidad>0
//	no me interesa el perfil interno
//	voy a ponerles el error igual a la varianza de la media de todos los casos
  double inner_profile[20],error_inner_profile[20];
  for(i=0;i<20;i++){
    inner_profile[i]=0.;
    error_inner_profile[i]=0.;}
  int cuantos_casos=0;
  for(int id_theory=0;id_theory<=6;id_theory++){
    for(j=0;j<bin_stack;j++){
      cuantos_casos++;
      for(i=0;i<20;i++){
        inner_profile[i]+=densidad[id_theory][j][i];
  }}}
  for(i=0;i<20;i++)
    inner_profile[i]/=cuantos_casos;
  for(int id_theory=0;id_theory<=6;id_theory++)
    for(j=0;j<bin_stack;j++)
      for(i=0;i<20;i++)
        error_inner_profile[i]+=pow(inner_profile[i]-densidad[id_theory][j][i],2);
  for(i=0;i<20;i++){
    error_inner_profile[i]/=cuantos_casos;
    printf("%d %le %le\n",i,sqrt(error_inner_profile[i]),inner_profile[i]);}
  for(int id_theory=0;id_theory<=6;id_theory++)
    for(j=0;j<bin_stack;j++)
      for(i=0;i<20;i++)
        error_densidad[id_theory][j][i]=3.*sqrt(error_inner_profile[i]);
#endif

  cout<<"echo"<<endl;
}









