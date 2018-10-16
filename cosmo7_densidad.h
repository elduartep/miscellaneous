//	lee los parametros ajustados por mcmc_todo_libre5nuevofast_2.c
//	tomado de		mcmc_todo_libre5nuevofast_2.c
void CargaAjustesDensidad(void){
  printf("leyendo ajuste densidad...");fflush(stdout);
  FILE * NOM;
  int n_salta_aux;
  char nomAr[210];
  if(MG==1)
    sprintf(nomAr,"../l/lista_todo_libre_%d_%d.txt",0,3);
  else
    sprintf(nomAr,"../l/lista_todo_libre_%d_%d.txt",3,6);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",&n_salta_aux,&co,&param_densidad[0][0],&param_densidad[0][1],&param_densidad[0][2],&param_densidad[0][3],&param_densidad[0][4],&param_densidad[1][0],&param_densidad[1][1],&param_densidad[1][2],&param_densidad[1][3],&param_densidad[1][4],&param_densidad[2][0],&param_densidad[2][1],&param_densidad[2][2],&param_densidad[2][3],&param_densidad[2][4])!=EOF);
  fclose(NOM);
  printf("hecho\n");fflush(stdout);
}






void ParametrosPerfil(int t,int s){
double x,x2;

  if(MG==1){
      x=log10(param_densidad[0][2]+pow(10.,param_cosmo[0]));
     x2=x*x;
  alpha=param_densidad[0][0] +param_densidad[0][1]*x
        +param_densidad[0][3]*sigma_stack[t][s];
   beta=param_densidad[0][4]*alpha;
      B=param_densidad[1][0]+param_densidad[1][1]*pow(sigma_stack[t][s],-param_densidad[1][2]);
        A=(param_densidad[2][0] +param_densidad[2][1]*x
         +param_densidad[2][3]*sigma_stack[t][s] 
         +param_densidad[2][4]*pow(sigma_stack[t][s],2)) *pow(B,1./param_densidad[0][4]);}
  else if(MG==2){
      x=param_cosmo[0];
     x2=x*x;
  alpha=param_densidad[0][0] +param_densidad[0][1]*x +param_densidad[0][2]*x2
        +param_densidad[0][3]*sigma_stack[t][s];
   beta=param_densidad[0][4]*alpha;
      B=param_densidad[1][0]+param_densidad[1][1]*pow(sigma_stack[t][s],-param_densidad[1][2]);
      A=(param_densidad[2][0] +param_densidad[2][1]*x +param_densidad[2][2]*x2
         +param_densidad[2][3]*sigma_stack[t][s] 
         +param_densidad[2][4]*pow(sigma_stack[t][s],2)) *pow(B,1./param_densidad[0][4]);}
}




inline double perfil(double r){
return -dc*(1.-A*pow(r,alpha))/(1.+B*pow(r,beta));
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









