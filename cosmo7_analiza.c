//coge los histogramas que calculó cosmo5.c para ver cuales la media, mejor ajuste, error
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;



int MG;//=1;		//	teoria de gravedad para calcular sigma: 1=f(R)		2=Symm
int si_densidad;//=1;	//	quiero usar densidad?
int si_abundancia;//=0;	//	quiero usar abundancia?
const	int      esferas=1;
int solo_analiza=0;
int  segunda_vez=0;

const int     theories=7;
int first_theory;//=3;
int  last_theory;//=3;
char NameTheories[]={'4','5','6','l','A','B','D'};
char NameTheory[]={'f','s'};
char NameObs[]={'a','p','b'};
int dash[]={2,5,1};
int color[]={1,4,7};
int grueso[]={1,2,5};
double porcentage[]={0.6827,0.9545,0.9973};


  const int	 NP=1000;
  double param_cosmo[NP];
  double chi2[NP];
  double integral[NP];










int main(int argc,char **argv){



  FILE * GNU;
  FILE * TOD;
  char nomAr[210];
  char nomOb[210];
  double gsup,ginf;
  
  sprintf(nomAr,"histo.help");
  TOD=fopen(nomAr,"w+");

  sprintf(nomAr,"cosmo7.gnu");
  GNU=fopen(nomAr,"w+");

  for(MG=1;MG<3;MG++){		//	clase de teoria
  for(first_theory=(MG-1)*3;first_theory<=(MG-1)*3+3;first_theory++){//	cual teoria
  last_theory=first_theory;


fprintf(GNU,"reset\n");
fprintf(GNU,"set terminal pngcairo size 1000,750 enhanced font 'TimesNewRoman,22' dashlength 1.7\n");
fprintf(GNU,"unset ytics\n");
fprintf(GNU,"set yr[0:1.05]\n");
if(MG==2) fprintf(GNU,"set xlabel 'z_{SSB}'\n");
if(MG==1) fprintf(GNU,"set xlabel 'log_{10}|f_{R0}|'\n");
if(first_theory==3)
fprintf(GNU,"set key t r\n");
else
fprintf(GNU,"unset key\n");
//fprintf(GNU,"set key t l\n");
if((MG==1)&&(first_theory==3))
fprintf(GNU,"set xtics 0.2\nset mxtics 4\n");
else
fprintf(GNU,"set xtics 0.1\nset mxtics 2\n");
if((MG==1)&&(first_theory==3)){
  fprintf(GNU,"set xlabel '|f_{R0}|'\n");
  fprintf(GNU,"set format x '%%2.0t{/Symbol \264}10^{%%L}'\n");
  fprintf(GNU,"set xtics 2e-8\n set mxtics 4\n");
}

  for(int obs=0;obs<3;obs++){
  if(obs==0){si_densidad=0;	si_abundancia=1;	sprintf(nomOb,"abundance");}
  if(obs==1){si_densidad=1;	si_abundancia=0;	sprintf(nomOb,"density profile");}
  if(obs==2){si_densidad=1;	si_abundancia=1;	sprintf(nomOb,"both");}

  FILE * NOM;


  sprintf(nomAr,"lista_%d_%d_%d_%d.txt",first_theory,MG,si_densidad,si_abundancia);
  NOM=fopen(nomAr,"r");


  int i,j;
  if(MG==1){
    if(first_theory==0){      ginf=-4.2;   gsup=-3.8;}
    if(first_theory==1){      ginf=-5.2;   gsup=-4.8;}
    if(first_theory==2){      ginf=-6.2;   gsup=-5.8;}
    if(first_theory==3){      ginf=0.;     gsup=7.e-8;}
  }
  else{
    if(first_theory==3){      ginf=0.;	gsup=0.4;}
    if(first_theory==4){      ginf=0.8;	gsup=1.2;}
    if(first_theory==5){      ginf=1.8;	gsup=2.2;}
    if(first_theory==6){      ginf=2.8;	gsup=3.2;}
  }

  double chi2min;
  int min;
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cosmo[i]);
    if(i==0){	integral[i]=exp(-0.5*chi2[i]);		chi2min=chi2[i];}
    else	integral[i]=integral[i-1]+exp(-0.5*chi2[i]);
    if(chi2min>=chi2[i]){	chi2min=chi2[i];	min=i;}
  }
  fclose(NOM);

  int med=0;
  while(integral[med]<0.5*integral[NP-1])
    med++;
  printf("%c%c %le %le ",NameTheory[MG-1],NameTheories[first_theory],param_cosmo[min],param_cosmo[med]);fflush(stdout);
for(j=0;j<3;j++){
  int izq=0;
  int der=NP-1;
  double fraccion=integral[der];
  while(fraccion>porcentage[j]*integral[NP-1]){
    if((first_theory==3)&&(si_abundancia==0)){
      fraccion=integral[der];
      der--;
    }
    else{
      while(chi2[der]>chi2[izq])
        der--;
      fraccion=integral[der]-integral[izq];
      izq++;
    }
  }
if((first_theory==3)&&(si_abundancia==0))
  printf("%le ",param_cosmo[der]);
else
  printf("%le %le ",param_cosmo[izq],param_cosmo[der]);
}
printf("\n");


if(obs!=0)
fprintf(GNU,"re");
if((MG==1)&&(first_theory==3))
fprintf(GNU,"plot '%s' u (10.**$2):(exp(-0.5*($1-%le))) w l lw 2 dt %d lc %d t '%s'\n",nomAr,chi2[min],dash[obs],color[obs],nomOb);
else
fprintf(GNU,"plot '%s' u 2:(exp(-0.5*($1-%le))) w l lw 2 dt %d lc %d t '%s'\n",nomAr,chi2[min],dash[obs],color[obs],nomOb);
fprintf(TOD,"%le\n",param_cosmo[min]);

}	//	observables

fprintf(GNU,"set xr[%le:%le]\n",ginf,gsup);
fprintf(GNU,"set output 'histo_%c%c.png'\n",NameTheory[MG-1],NameTheories[first_theory]);
fprintf(GNU,"replot\n");
fprintf(GNU,"reset\n\n\n");

}	//	teoria

}	//	clase de teoria
fclose(GNU);
fclose(TOD);
  cout<<"terminó"<<endl;
  return 0;
}

