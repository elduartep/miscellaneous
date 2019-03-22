//aqui lo unico que necesito es calcular la derivada de sigma e interpolarla
//interpolar el propio sigma
//usar estos dos valores y la formula para la abundancia que ajusté previamente
//dar como resultado la abundancia para los radios en los que medí la abundancia en las simulaciones


void ParametrosAbundancia(int id_theory){
double x;

  if(MG==1){
    x=log10(pow(10.,param_cosmo[0])+param_abundancia[2]);
    gama=param_abundancia[0]+param_abundancia[1]*x;}
  else{
    x=param_cosmo[0];
    gama=param_abundancia[0]+param_abundancia[1]*x/(1.+param_abundancia[2]*exp(-x*x));}
}



//	A nu^alpha ( 1 + b nu^beta) exp( -c nu^gama )
inline double FuncionRadio(double sigma){
double nu=1.686/sigma;
return param_abundancia[3]*pow(sigma,-gama)*(1.+pow(nu,param_abundancia[4]))
        *exp(-param_abundancia[5]*pow(sigma,-2.));
}



// lee el archivo que contiene la función de radio que fué calculada por FuncionRadio.c
void LeeAbundancia(void){
  printf("Leyendo funcion radio...");fflush(stdout);
  int i;
  char NomArch[300];
  char line[450];
  FILE * Ud;
  for(int id_theory=0;id_theory<=6;id_theory++){
    sprintf(NomArch,"../%c/dnv_esferas_%s%s",NameTheories[id_theory],prefijo,caso);
    Ud=fopen(NomArch,"r");

    for(i=0;i<bin_abundancia;i++){
      fscanf(Ud,"%le %le %le\n",
        &radio_abundancia[id_theory][i], &abundancia[id_theory][i], &error_abundancia[id_theory][i]);
    }
    fclose(Ud);
  }
  printf("echo\n");fflush(stdout);
}







void CargaAjustesAbundancia(void){
  printf("leyendo ajuste abundancia...");fflush(stdout);
  FILE * NOM;
  char nomAr[410];
  double aux;int au;


  if(MG==1)
    sprintf(nomAr,"../l/list_abundancia_0_3.txt");
  else if(MG==2)
    sprintf(nomAr,"../l/list_abundancia_3_6.txt");printf("%s\n",nomAr);
  NOM=fopen(nomAr,"r");
  while((fscanf(NOM,"%d %le %le %le %le %le %le %le\n",&au,&aux
    ,&param_abundancia[0],&param_abundancia[1],&param_abundancia[2]
    ,&param_abundancia[3],&param_abundancia[4],&param_abundancia[5]))!=EOF);

  fclose(NOM);
  printf("echo\n");fflush(stdout);
}






void calcula_sigma_derivadas(int id_theory){
Int i;
//	calcula sigma y derivada para todos los radios en que tenemos abundancia
double w,dw,wm,dwm,sigma2;int j;
  for(j=0;j<bin_abundancia;j++){
    sigma2 = 0.0;
    dsigma_abundancia[id_theory][j] = 0.0;
    for(i=0;i<nk-1;i++){///	integral
      W_dW(radio_abundancia[id_theory][j]*K[i],w,dw);
      W_dW(radio_abundancia[id_theory][j]*K[i+1],wm,dwm);
       sigma2     += (K[i+1]-K[i])*(P[i]*pow(K[i],2)*pow(w,2) + P[i+1]*pow(K[i+1],2)*pow(wm,2));
      dsigma_abundancia[id_theory][j] += (K[i+1]-K[i])*(P[i]*pow(K[i],3)*w*dw + P[i+1]*pow(K[i+1],3)*wm*dwm);
    }
    sigma_abundancia[id_theory][j]=sqrt(sigma2)/(2.0*pi);
    dsigma_abundancia[id_theory][j]/=(4.0*pi*pi);
    }

}






void calcula_abundancia_teorica(int id_theory){
int i;
for(i=1;i<bin_abundancia;i++)
  teoria_abundancia[id_theory][i]=
        -3.*dsigma_abundancia[id_theory][i]/
         (4.*pi*pow(sigma_abundancia[id_theory][i]*radio_abundancia[id_theory][i],2))
        *FuncionRadio(sigma_abundancia[id_theory][i]);
}








