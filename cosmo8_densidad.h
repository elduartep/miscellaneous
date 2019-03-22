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
  char nomAr[210];
  if(MG==1)
    sprintf(nomAr,"../l/lista_bias_densidad_%d_%d.txt",0,3);
  else
    sprintf(nomAr,"../l/lista_bias_densidad_%d_%d.txt",3,6);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le %le\n",&n_salta_aux,&co,
&param_bd[0],&param_bd[1],&param_bd[2],&param_bd[3],&param_bd[4],&param_bd[5],&param_bd[6],&param_bd[7],&param_bd[8],&param_bd[9],&param_bd[10],&param_bd[11])!=EOF);
  fclose(NOM);
  printf("hecho\n");fflush(stdout);
}






//	f(R)
//	c=-0.15895 + 1.75e-4 log10(3e-10 + fR0)
//	symmetron
//	c=-0.1705 + 6.333e-4 z + ..zz
//	bias= a + b sigma_2 + c sigam_4
void ParametrosBias(int id_theory){
double x,x2;
  if(MG==1){
    x=log10(pow(10.,param_cosmo[0])+param_bd[Npd+2]);
    b3=param_bd[Npd+0]+param_bd[Npd+1]*x;}
  else if(MG==2){
    x=param_cosmo[0];
    b3=param_bd[Npd+0]+param_bd[Npd+1]*x+param_bd[Npd+2]*x*x;}
  b2=param_bd[Npd+3]+param_bd[Npd+4]*b3;
  b1=param_bd[Npd+5];
}



//	lo mismo de arriba, solo que para desciende
void coeficientes(int t){
double x,x2;
  if(first_theory==0){
    x=log10(parametro[t]+param_bd[Npd+2]);
    b3=param_bd[Npd+0]+param_bd[Npd+1]*x;}
  else if(first_theory==3){
    x=parametro[t];
    b3=param_bd[Npd+0]+param_bd[Npd+1]*x+param_bd[Npd+2]*x*x;}
  b2=param_bd[Npd+3]+param_bd[Npd+4]*b3;
  b1=param_bd[Npd+5];
}






inline double bias(double sigma){
return b1 + b2*pow(sigma,-2.)+ b3*pow(sigma,-4.);
}






//	tomado de		mcmc_todo_libre5nuevofast_2.c
// lee el archivo que contiene los perfiles
void LeeBias(void){
  int i,j;
  double aux;
  printf("leyendo bias...");fflush(stdout);
  char NomArch[300];
  FILE * Ud;

for(int id_theory=0;id_theory<=6;id_theory++){
  sprintf(NomArch,"../%c/bv_03.txt",NameTheories[id_theory]);

  Ud=fopen(NomArch,"r");
  for(j=0;j<bin_bias;j++)
    fscanf(Ud,"%lf %lf %le %le\n"
      ,&radio_bias[id_theory][j], &aux, &medida_bias[id_theory][j], &error_bias[id_theory][j]);

  fclose(Ud);
}
  printf("hecho\n");fflush(stdout);
}




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////	densidad
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////




void ParametrosPerfil(int t, int s){
 A=param_bd[4]*pow(radio_stack[t][s]/6.9,param_bd[5]);
 B=param_bd[6]*pow(radio_stack[t][s]/6.9,param_bd[7]);
}




void PerfilTeorico(double r){
double universal= -dc*(1.-pow(r/param_bd[2],param_bd[0]))/(1.+pow(r/param_bd[3],param_bd[1]));
//	en el caso del side-descendent, dados los parametros del presente paso puedo calcular el perfil de todas las teorias sin intepolar
return universal;
}






//	tomado de		mcmc_todo_libre5nuevofast_2.c
// lee el archivo que contiene los perfiles
void LeePerfiles(void){
  int i,j;
  double aux;
  cout<<"Leyendo perfiles...";
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
    for(i=0;i<bin_radio_densidad;i++){//printf("%d %le %le %le\n",j,radio_densidad[i],densidad[j][i],error_densidad[j][i]);fflush(stdout);
      if(error_densidad[id_theory][j][i]<0.025){
         error_densidad[id_theory][j][i]=0.025;
      }
    }
//	fijando el valor de un punto
densidad[id_theory][j][0+j]=-dc;		///////////////////
error_densidad[id_theory][j][0+j]=0.005;		///////////////////
  }
#endif
  fclose(Ud);
}
  cout<<"echo"<<endl;
}









