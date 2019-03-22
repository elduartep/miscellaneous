//coge los histogramas que calculó cosmo5.c para ver cuales la media, mejor ajuste, error
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "../l/nr3.h"
#include "../l/cosmo7_parametros.h"

int cruzado=1;

//int MG;//=1;		//	teoria de gravedad para calcular sigma: 1=f(R)		2=Symm
//int si_densidad;//=1;	//	quiero usar densidad?
//int si_abundancia;//=0;	//	quiero usar abundancia?
//int si_bias;//=0;	//	quiero usar bias?
const	int      esferas=1;
int solo_analiza=0;
int  segunda_vez=0;

//const int     theories=7;
int first_theory;//=3;
int  last_theory;//=3;
//char NameTheories[]={'D','B','A','l','6','5','4'};
const int param[]={3,2,1,0,-6,-5,-4};
char NameTheories_correct[]={'4','5','6','l','A','B','D'};
char NameTheory[]={'f','s'};
char NameObs[]={'n','d','b','a'};
int dash[]={2,3,5,1};
int color[]={1,2,4,7};
int grueso[]={2,3,2,2};
double porcentage[]={0.6827,0.9545,0.9973};



  const int	 NP=1000;
  double param_cos[NP];
  double chi2[NP];
  double integral[NP];






double ar(double x){
return roundf(x * 100) / 100;}



int main(int argc,char **argv){


  int dof[4]={0};
    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++)
      for(int i=0;i<bines_radio_densidad[id_stack];i++)
        dof[1]++;
    for(int i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      dof[0]++;
    for(int id_bin=first_bin_bias;id_bin<=last_bin_bias;id_bin++)
      dof[2]++;
  dof[3]=dof[0]+dof[1]+dof[2];


  FILE * GNU;
  FILE * TOD;
  FILE * TEX;
  char nomAr[210];
  char nomOb[210];
  double gsup,ginf;
  
  sprintf(nomAr,"histo.help");
  TOD=fopen(nomAr,"w+");

  sprintf(nomAr,"cosmo72.gnu");
  GNU=fopen(nomAr,"w+");



  for(MG=1;MG<3;MG++){		//	clase de teoria
  for(first_theory=(MG-1)*3;first_theory<=(MG-1)*3+3;first_theory++){//	cual teoria
  last_theory=first_theory;


  if(MG==1){	//	para las cosmologicas cruzadas
    if(first_theory==0){      ginf=-4.8;   gsup=-2.8;}
    if(first_theory==1){      ginf=-6.0;   gsup=-4.;}
    if(first_theory==2){      ginf=-10.0;   gsup=-5.;}
    if(first_theory==3){      ginf=0.;     gsup=7.e-8;}
  }
  else{
    if(first_theory==3){      ginf=0.;	gsup=1.5;}
    if(first_theory==4){      ginf=0.;	gsup=2.;}
    if(first_theory==5){      ginf=1.;	gsup=3.;}
    if(first_theory==6){      ginf=1.6;	gsup=3.6;}
  }

/*  if(MG==1){//	para las cosmologias correctas
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
  }*/

  sprintf(nomAr,"table_%c_%c.tex",NameTheory[MG-1],NameTheories[6-first_theory]);
  TEX=fopen(nomAr,"w");

fprintf(GNU,"reset\n");
fprintf(GNU,"set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7\n");
fprintf(GNU,"unset ytics\n");
fprintf(GNU,"set yr[0:1.05]\n");

if(MG==2) fprintf(GNU,"set xlabel 'z_{SSB}'\n");
if(MG==1) fprintf(GNU,"set xlabel 'log_{10}|f_{R0}|'\n");
if(cruzado==1){
if(MG==2) fprintf(GNU,"set label '|f_{R0}|=10^{%d}' at %.1lf,%.1lf\n",param[last_theory],ginf+0.05*(gsup-ginf),0.9);
if(MG==1) fprintf(GNU,"set label 'z_{SSB}=%d' at %.1lf,%.1lf\n",param[last_theory],ginf+0.05*(gsup-ginf),0.9);}

if((first_theory==2)&&(MG==1))
fprintf(GNU,"set key t r\n");
else if((first_theory==4)&&(MG==2))
fprintf(GNU,"set key t r\n");
else
fprintf(GNU,"unset key\n");
//fprintf(GNU,"set key t l\n");
if((MG==1)&&(first_theory==3))
fprintf(GNU,"set xtics 0.5\nset mxtics 5\n");
else
fprintf(GNU,"set xtics 0.2\nset mxtics 2\n");
if((MG==1)&&(first_theory==3)){
  fprintf(GNU,"set xlabel '|f_{R0}| {/Symbol \264} 10^8'\n");
  if(cruzado==1)
  fprintf(GNU,"set label 'z_{SSB}=0' at %.1lf,%.1lf\n",ginf+0.05*(gsup-ginf),0.9);
  fprintf(GNU,"set format x '%%1.0t'\n");
  fprintf(GNU,"set xtics 1e-8\n set mxtics 2\n");
}

  for(int obs=0;obs<4;obs++){
  if(obs==0){si_densidad=0; si_abundancia=1; si_bias=0; sprintf(nomOb,"abundance");
    fprintf(TEX,"%s & ",nomOb);}
  if(obs==1){si_densidad=1; si_abundancia=0; si_bias=0; sprintf(nomOb,"density profile");
    fprintf(TEX,"%s & ",nomOb);}
  if(obs==2){si_densidad=0; si_abundancia=0; si_bias=1; sprintf(nomOb,"bias");
    fprintf(TEX,"%s & ",nomOb);}
  if(obs==3){si_densidad=1; si_abundancia=1; si_bias=1; sprintf(nomOb,"all");
    fprintf(TEX,"%s & ",nomOb);}

  FILE * NOM;

  sprintf(nomAr,"lista2_%d_%d_%d_%d_%d.txt",first_theory,MG,si_densidad,si_abundancia,si_bias);
  NOM=fopen(nomAr,"r");



  int i,j;


  double chi2min;
  int min;

//	buscando el minimo, es necesario antes de calcular la integra porque normaliza los dados
//	de otro nodo casi todo contribuye con 0, y el pico es 10^-un_monton
if((MG==1)&&(first_theory==3))
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cos[i]);
    if(i==0){ chi2min=chi2[0]; min=0;}
    if(chi2min>chi2[i]){	chi2min=chi2[i];	min=i;}
  }
else
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cos[i]);
    if(i==0){ chi2min=chi2[0]; min=0;}
    if(chi2min>chi2[i]){	chi2min=chi2[i];	min=i;}
  }
  rewind(NOM);

//	ahora si estimamos la integral de la distribucion
if((MG==1)&&(first_theory==3))
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cos[i]);
    if(i==0)	integral[i]=exp(-0.5*(chi2[i]-chi2min))*pow(10.,param_cos[i]);
    else 	integral[i]=integral[i-1]+exp(-0.5*(chi2[i]-chi2min))*(pow(10.,param_cos[i])-pow(10.,param_cos[i-1]));
  }
else
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cos[i]);
    if(i==0)	integral[i]=exp(-0.5*(chi2[i]-chi2min));
    else	integral[i]=integral[i-1]+exp(-0.5*(chi2[i]-chi2min));
  }


  fclose(NOM);
if(obs==0) printf("\n");
printf("chimin=%lf param=%le id=%d ",chi2min/dof[obs],param_cos[min],min);

  int med=0;
  while(integral[med]<0.5*integral[NP-1])
    med++;
  if((MG==1)&&(first_theory==3)){
    printf("%c%c %le %le ",NameTheory[MG-1],NameTheories[6-first_theory],
      pow(10.,param_cos[min]),pow(10,param_cos[med]));
    fprintf(TEX,"%.2lf & %.2lf ",ar(pow(10.,param_cos[min]+8.)),ar(pow(10,param_cos[med]+8.)));}
  else{
    printf("%c%c %le %le ",NameTheory[MG-1],NameTheories[6-first_theory],param_cos[min],param_cos[med]);
    fflush(stdout);
    fprintf(TEX,"%.2lf & %.2lf ",ar(param_cos[min]),ar(param_cos[med]));
    fflush(stdout);}


for(j=0;j<3;j++){
  int izq=0;
  int der=NP-1;
  double fraccion=integral[der];
  if((first_theory==3)){
    while(integral[der]>porcentage[j]*integral[NP-1]){
      der--;
    }
  if(MG==1){
    printf("%le ",pow(10.,param_cos[der]));
    fprintf(TEX,"& <%.2lf ",ar(pow(10.,param_cos[der]+8.)));
}
  else{
    printf("%le ",param_cos[der]);
    fprintf(TEX,"& <%.2lf ",ar(param_cos[der]));}

  }
  else{
    while(fraccion>porcentage[j]*integral[NP-1]){
      while(chi2[der]>chi2[izq])
        der--;
      fraccion=integral[der]-integral[izq];
      izq++;
    }
    printf("(%le %le %le) ",0.5*(-param_cos[izq]+param_cos[der]),param_cos[izq],param_cos[der]);
    fprintf(TEX,"& $\\pm$ %.2lf ",ar(0.5*(-param_cos[izq]+param_cos[der])));
  }
}
printf("\n");
fprintf(TEX,"\\\\ \n");


if(obs!=0)
fprintf(GNU,"re");
if((MG==1)&&(first_theory==3))
fprintf(GNU,"plot '%s' u (10.**$2):(exp(-0.5*($1-%le))) w l dt %d lc %d lw %d t '%s'\n",nomAr,chi2[min],dash[obs],color[obs],grueso[obs],nomOb);
else
fprintf(GNU,"plot '%s' u 2:(exp(-0.5*($1-%le))) w l dt %d lc %d lw %d t '%s'\n",nomAr,chi2[min],dash[obs],color[obs],grueso[obs],nomOb);
fprintf(TOD,"%le\n",param_cos[min]);

}	//	observables

fprintf(GNU,"set xr[%le:%le]\n",ginf,gsup);
fprintf(GNU,"set output 'histo_%c%c.png'\n",NameTheory[MG-1],NameTheories[6-first_theory]);
fprintf(GNU,"replot\n");
fprintf(GNU,"reset\n\n\n");
fprintf(TEX,"\n");
fclose(TEX);
}	//	teoria
printf("\n");
}	//	clase de teoria
fclose(GNU);
fclose(TOD);


  cout<<"terminó"<<endl;
  return 0;
}

