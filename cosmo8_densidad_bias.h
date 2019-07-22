////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////	bias
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


//	lee el resultado de cosmo8_desciende_bias_densidad.c
//	incompleto
void CargaAjustesBiasDensidad(void){
  printf("leyendo ajuste bias-densidad...\n");fflush(stdout);
  FILE * NOM;
  int n_salta_aux;
  char nomAr[200];
  if(MG==1)
    sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",0,3,0,1);
  else
    sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",3,6,0,1);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",&n_salta_aux,&co,
&param_densidad[0],&param_densidad[1],&param_densidad[2],&param_densidad[3],&param_densidad[4],&param_densidad[5],&param_densidad[6],&param_densidad[7],&param_densidad[8],&param_densidad[9],&param_densidad[10],&param_densidad[11],&param_densidad[12])!=EOF);
  fclose(NOM);
printf("[0]=%le\n",param_densidad[0]);
printf("[12]=%le\n",param_densidad[12]);
  printf("hecho\n");fflush(stdout);
}



//	cuando voy a ajustar la densidad debo usar el resultado del bias
void CargaAjustesBias(void){
  printf("leyendo ajuste bias...\n");fflush(stdout);
  FILE * NOM;
  double aux;
  int n_salta_aux;
  char nomAr[200];
  if(MG==1)
    sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",0,3,1,0);
  else
    sprintf(nomAr,"lista_bias_densidad_t%d-%d_b%d_d%d.txt",3,6,1,0);
  NOM=fopen(nomAr,"r");
  while(fscanf(NOM,"%d %le %le %le %le %le\n",&n_salta_aux,&co,&b1,&b2,&b0,&b4)!=EOF);
  fclose(NOM);
  printf("b1=%le b2=%le b0=%le b4=%le\n",b1,b2,b0,b4);

  printf("hecho\n");fflush(stdout);
}








inline double bias(double sigma,int t){
  double x,x2;
  if(first_theory==0)
      x=log10(parametro[t]);//pow(10.,param_cosmo[t]));
  else if(first_theory==3)
      x=parametro[t];//param_cosmo[t];
  b3= param_bd[0] +param_bd[1]*x;
  b0=param_bd[2];
  b4=param_bd[3];
  return b0 + b3 *pow(sigma,-2.)+ b4 *pow(sigma,-4.);
}

inline double bias_cosmo(double sigma,int t){
  double x,x2;
  x=param_cosmo[0];
  b3 = b1 +b2*x;
  return b0 + b3 *pow(sigma,-2.) +b4 *pow(sigma,-4.);
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
//C=param_bd[6]*pow(radio_stack[t][s],+param_bd[7]);

//recta y parabola
//A=param_bd[4]+pow(radio_stack[t][s]-6.9,2)*param_bd[5];
//B=param_bd[6]+radio_stack[t][s]*param_bd[7];

//cuando todo es libre, completamente libre, incluyendo el bias
//A=param_bd[t] *radio_stack[t][s] /(1.0+ radio_stack[t][s]);
//C=param_bd[t+7]*tanh(radio_stack[t][s]-7.0)+param_bd[14];
}


double all(double r,int id_theory,int id_stack){
//exponential (plus bump)
A=param_bd[3]*pow(radio_stack[id_theory][id_stack]-param_bd[4],3)+param_bd[5];
C=alpha0[id_theory]*pow(sin(radio_stack[id_theory][id_stack]/param_bd[2]/pi),alpha0[id_theory]*0.1);
//C=2.0*alpha0[id_theory]*param_bd[2]-pow((radio_stack[id_theory][id_stack]-param_bd[10])*alpha0[id_theory],2)/param_bd[2];
return exp(-pow(A*r,-C));
}



double Suppress(double y,int id_theory,int id_stack){
A=param_bd[4]*pow(radio_stack[id_theory][id_stack],-param_bd[5]);
B=param_bd[6]*pow(radio_stack[id_theory][id_stack],+param_bd[7]);
return exp(-A*pow(y,-B));
}
double alpha=7.593649e+00;
double beta=9.381333e+00;
double ca=1.415600e+00;
double cb=1.106479e+00;
double PerfilUniversal(double r,int id_theory,int id_stack){
//	universal (plus supppres)
return -dc*(1.-pow(r/param_bd[0],param_bd[1]))/(1.+pow(r/param_bd[2],param_bd[3]));
//return -dc*(1.-pow(r/ca,alpha))/(1.+pow(r/cb,beta));



/*	//double power law
double a,alpha,b,beta;
//alpha = alpha0[id_theory]*(1.-exp(-pow(radio_stack[id_theory][id_stack]/param_bd[2],param_bd[3])));
alpha = alpha0[id_theory]*(param_bd[2]-param_bd[3]*pow(radio_stack[id_theory][id_stack]/param_bd[4]-1.,2));
    a = param_bd[id_theory+10]*param_bd[9] +pow(param_bd[id_theory+10] /radio_stack[id_theory][id_stack],param_bd[8]);
 beta = alpha+param_bd[5];
    b = param_bd[6]+param_bd[7]/radio_stack[id_theory][id_stack];
return (pow(r/a,alpha) +pow(r/b,beta)) /(1.0 +pow(r/b,beta));
*/

}


void ParametrosPerfil(int t,int s){
double x,x2;
  if(first_theory==0){
      x=log10(param_bd[2]+parametro[t]);
  alpha= param_bd[0] +param_bd[1]*x
        +param_bd[3]*sigma_bias[t][s];
   beta= param_bd[4]*alpha;
      B= param_bd[6]*pow(sigma_bias[t][s],-param_bd[7]);
      A=(param_bd[8]+param_bd[9]*x
        +param_bd[11]*log10(sigma_bias[t][s]) ) *pow(B,1./param_bd[4]);}
  else if(first_theory==3){
      x=parametro[t];
     x2=x*x;
  alpha= param_bd[0] +param_bd[1]*x// +param_bd[2]*x2
        +param_bd[3]*sigma_bias[t][s];
   beta= param_bd[4]*alpha;
      B= param_bd[6]*pow(sigma_bias[t][s],-param_bd[7]);
      A=(param_bd[8]+param_bd[9]*x //+param_bd[10]*x2
        +param_bd[11]*log10(sigma_bias[t][s]) ) *pow(B,1./param_bd[4]);}
}
inline double perfil(double r){
return -dc*(1.-A*pow(r,alpha))/(1.+B*pow(r,beta));
}
void ParametrosPerfilCosmo(int t,int s){
double x,x2;
  if(MG==1){
      x=log10(param_densidad[2]+pow(10.,param_cosmo[0]));
  alpha= param_densidad[0] +param_densidad[1]*x
        +param_densidad[3]*sigma_bias[t][s];
   beta= param_densidad[4]*alpha;
      B= param_densidad[6]*pow(sigma_bias[t][s],-param_densidad[7]);
      A=(param_densidad[8]+param_densidad[9]*x
        +param_densidad[11]*log10(sigma_bias[t][s])) *pow(B,1./param_densidad[4]);}
  else if(MG==2){
      x=param_cosmo[0];
     x2=x*x;
  alpha= param_densidad[0] +param_densidad[1]*x +param_densidad[2]*x2
        +param_densidad[3]*sigma_bias[t][s];
   beta= param_densidad[4]*alpha;
      B= param_densidad[6]*pow(sigma_bias[t][s],-param_densidad[7]);
      A=(param_densidad[8]+param_densidad[9]*x +param_densidad[10]*x2
        +param_densidad[11]*log10(sigma_bias[t][s])) *pow(B,1./param_densidad[4]);}
//printf("MG=%d t=%d s=%d a=%le b=%le A=%le B=%le\n",MG,t,s,alpha,beta,A,B);
}



/*
void ParametrosPerfil(int t,int s){
double x,x2;
  if(first_theory==0){
      x=log10(param_bd[2]+parametro[t]);
  alpha= param_bd[0] +param_bd[1]*x
        +param_bd[3]*sigma_bias[t][s];
   beta= param_bd[4]*alpha;
      B= param_bd[5]+param_bd[6]*pow(sigma_bias[t][s],-param_bd[7]);
      A=(param_bd[8]+param_bd[9]*x
        +param_bd[11]*sigma_bias[t][s] 
        +param_bd[12]*pow(sigma_bias[t][s],2)) *pow(B,1./param_bd[4]);}
  else if(first_theory==3){
      x=parametro[t];
     x2=x*x;
  alpha= param_bd[0] +param_bd[1]*x +param_bd[2]*x2
        +param_bd[3]*sigma_bias[t][s];
   beta= param_bd[4]*alpha;
      B= param_bd[5]+param_bd[6]*pow(sigma_bias[t][s],-param_bd[7]);
      A=(param_bd[8]+param_bd[9]*x +param_bd[10]*x2
        +param_bd[11]*sigma_bias[t][s]
        +param_bd[12]*pow(sigma_bias[t][s],2)) *pow(B,1./param_bd[4]);}
}
inline double perfil(double r){
return -dc*(1.-A*pow(r,alpha))/(1.+B*pow(r,beta));
}
void ParametrosPerfilCosmo(int t,int s){
double x,x2;
  if(MG==1){
      x=log10(param_densidad[2]+pow(10.,param_cosmo[0]));
  alpha= param_densidad[0] +param_densidad[1]*x
        +param_densidad[3]*sigma_bias[t][s];
   beta= param_densidad[4]*alpha;
      B= param_densidad[5]+param_densidad[6]*pow(sigma_bias[t][s],-param_densidad[7]);
      A=(param_densidad[8]+param_densidad[9]*x
        +param_densidad[11]*sigma_bias[t][s] 
        +param_densidad[12]*pow(sigma_bias[t][s],2)) *pow(B,1./param_densidad[4]);}
  else if(MG==2){
      x=param_cosmo[0];
     x2=x*x;
  alpha= param_densidad[0] +param_densidad[1]*x +param_densidad[2]*x2
        +param_densidad[3]*sigma_bias[t][s];
   beta= param_densidad[4]*alpha;
      B= param_densidad[5]+param_densidad[6]*pow(sigma_bias[t][s],-param_densidad[7]);
      A=(param_densidad[8]+param_densidad[9]*x +param_densidad[10]*x2
        +param_densidad[11]*sigma_bias[t][s]
        +param_densidad[12]*pow(sigma_bias[t][s],2)) *pow(B,1./param_densidad[4]);}
//printf("MG=%d t=%d s=%d a=%le b=%le A=%le B=%le\n",MG,t,s,alpha,beta,A,B);
}
*/


void CargaCorrelacion(int id_theory){
  for(int i=0;i<max_id_xi;i++)
    xi_actual[i]=xi[id_theory][i];
}


void CargaSigma(int id_theory){
  for(int i=0;i<max_id_xi;i++)
    sig_actual[i]=sig[id_theory][i];
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


#if uniformiza_densidad>0
  for(j=0;j<bin_stack;j++){
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

}

/*
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
#endif*/


  cout<<"echo"<<endl;
}









