//	lee el resultado de cosmo7_desciende_bias.c
//	incompleto
void CargaAjustesBias(void){
  printf("leyendo ajuste bias...");fflush(stdout);
  FILE * NOM;
  int n_salta_aux;
  char nomAr[210];
  if(MG==1)
    sprintf(nomAr,"../l/lista_bias_%d_%d.txt",0,3);
  else
    sprintf(nomAr,"../l/lista_bias_%d_%d.txt",3,6);
  NOM=fopen(nomAr,"r");
//sin parametros libres que dependen de la teoria
//  while(fscanf(NOM,"%d %le %le %le %le\n",
//&n_salta_aux,&co,&param_bias[0],&param_bias[1],&param_bias[2])!=EOF);
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le\n",&n_salta_aux,&co,
&param_bias[0],&param_bias[1],&param_bias[2],&param_bias[3],&param_bias[4],&param_bias[5])!=EOF);
  fclose(NOM);
  printf("hecho\n");fflush(stdout);
}





//	sin parametros libres que dependen de la teoria
//	bias= a + b sigma_2 + c sigam_4
/*void ParametrosBias(int id_theory){
  b1=param_bias[0];
  b2=param_bias[1];
  b3=param_bias[2];
}*/



//	f(R)
//	c=-0.15895 + 1.75e-4 log10(3e-10 + fR0)
//	symmetron
//	c=-0.1705 + 6.333e-4 z + ..zz
//	bias= a + b sigma_2 + c sigam_4
void ParametrosBias(int id_theory){
double x,x2;
  if(MG==1){
    x=log10(pow(10.,param_cosmo[0])+param_bias[2]);
    b3=param_bias[0]+param_bias[1]*x;}
  else if(MG==2){
    x=param_cosmo[0];
    b3=param_bias[0]+param_bias[1]*x+param_bias[2]*x*x;}
  b2=param_bias[3]+param_bias[4]*b3;
  b1=param_bias[5];
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









