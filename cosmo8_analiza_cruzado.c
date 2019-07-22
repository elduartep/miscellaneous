//coge los histogramas que calculó cosmo8.c para ver cuales la media, mejor ajuste, error
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "../l/nr3.h"
#include "../l/cosmo8_parametros.h"

int cruzado;


const	int      esferas=1;
int solo_analiza=0;
int  segunda_vez=0;

const int param[]={-4,-5,-6,0,1,2,3};
char NameTheories_correct[]={'4','5','6','l','A','B','D'};
char NameTheory[]={'f','s'};
char NameObs[]={'n','d','b','a'};
int dash[]={2,3,5,1};
int color[]={1,2,4,7};
int grueso[]={2,3,2,2};
double porcentage[]={0.6827,0.9545,0.9973};



  const int	 NP=500;
  double param_cos[NP];
  double chi2[NP];
  double integral[NP];

  double chi2_dof[2][7];
  double delta_cr[2][7];
  double sigma_cr[2][7];




double ar(double x){
return roundf(x * 100) / 100;}



int main(int argc,char **argv){


  int dof[4]={0};
    for(int id_stack=first_stack;id_stack<=last_stack;id_stack++)
      for(int i=0;i<bines_radio_densidad[id_stack];i++)
        dof[1]++;
    dof[1]-=14;
    for(int i=first_bin_abundancia;i<=last_bin_abundancia;i++)
      dof[0]++;
    for(int id_bin=first_stack;id_bin<=last_stack;id_bin++)
      dof[2]++;
  dof[3]=dof[0]+dof[1]+dof[2];


  FILE * GNU;
  FILE * TOD;
  FILE * TEX;
  FILE * CHX;
  char nomAr[210];
  char nomOb[210];
  double gsup,ginf;

  ////////////////////////////////////////////////////////////////////////////////
  for(cruzado=2;cruzado>0;cruzado--){
  
  sprintf(nomAr,"histo.help");
  TOD=fopen(nomAr,"w+");

  if(cruzado==1)
  sprintf(nomAr,"cosmo7c.gnu");
  else
  sprintf(nomAr,"cosmo7.gnu");
  GNU=fopen(nomAr,"w+");


  

  for(MG=1;MG<3;MG++){		//	clase de teoria
  for(first_theory=(MG-1)*3;first_theory<=(MG-1)*3+3;first_theory++){//	cual teoria
  last_theory=first_theory;



  if(MG==1){//	para las cosmologias correctas
    if(first_theory==0){      ginf=-5.;   gsup=-3.;}
    if(first_theory==1){      ginf=-6.;   gsup=-4.;}
    if(first_theory==2){      ginf=-7.;   gsup=-5.;}
    if(first_theory==3){      ginf=-10;   gsup=-5.5;}
  }
  else{
    if(first_theory==3){      ginf=0.;	gsup=1.;}
    if(first_theory==4){      ginf=0.;	gsup=2.;}
    if(first_theory==5){      ginf=1.;	gsup=3.;}
    if(first_theory==6){      ginf=2.;	gsup=4.;}
  }

  if(cruzado==1)
  sprintf(nomAr,"table_%c_%c.tex",NameTheory[MG-1],NameTheories[6-first_theory]);
  else
  sprintf(nomAr,"table_%c_%c.tex",NameTheory[MG-1],NameTheories[first_theory]);
  TEX=fopen(nomAr,"w");

fprintf(GNU,"reset\n");
fprintf(GNU,"set terminal pngcairo size 1000,600 enhanced font 'TimesNewRoman,22' dashlength 1.7\n");
fprintf(GNU,"unset ytics\n");
fprintf(GNU,"set yr[0:1.05]\n");

if(MG==2) fprintf(GNU,"set xlabel 'z_{SSB}'\n");
if(MG==1) fprintf(GNU,"set xlabel 'log_{10}|f_{R0}|'\n");

////	label
if(cruzado==1){
if(MG==2) fprintf(GNU,"set label '|f_{R0}|=10^{%d}' at %.1lf,%.1lf\n",param[6-last_theory],ginf+0.8*(gsup-ginf),0.5);
if(MG==1) fprintf(GNU,"set label 'z_{SSB}=%d' at %.1lf,%.1lf\n",param[6-last_theory],ginf+0.8*(gsup-ginf),0.5);}
else{
if(first_theory==3) fprintf(GNU,"set label '{/Symbol L}CDM' at %.1lf,%.1lf\n",ginf+0.8*(gsup-ginf),0.4);
else if(MG==1) fprintf(GNU,"set label '|f_{R0}|=10^{%d}' at %.1lf,%.1lf\n",param[last_theory],ginf+0.8*(gsup-ginf),0.4);
else if(MG==2) fprintf(GNU,"set label 'z_{SSB}=%d' at %.1lf,%.1lf\n",param[last_theory],ginf+0.8*(gsup-ginf),0.4);
}

////	key
if(cruzado==1){
  if((first_theory==2)&&(MG==1))
  fprintf(GNU,"set key t r\n");
  else if((first_theory==4)&&(MG==2))
  fprintf(GNU,"set key t r\n");
  else
  fprintf(GNU,"unset key\n");}
else{
  if((first_theory==3)&&(MG==1))
  fprintf(GNU,"set key t r\n");
  else if((first_theory==3)&&(MG==2))
  fprintf(GNU,"set key t r\n");
  else
  fprintf(GNU,"unset key\n");}


fprintf(GNU,"set xtics 0.5\nset mxtics 5\n");


//////////////////////////////////////////////////////////////
	//	begin observables
//////////////////////////////////////////////////////////////

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

  if(cruzado==1)
  sprintf(nomAr,"lista8c_%d_%d_%d_%d_%d.txt",first_theory,MG,si_densidad,si_abundancia,si_bias);
  else
  sprintf(nomAr,"lista8_%d_%d_%d_%d_%d.txt",first_theory,MG,si_densidad,si_abundancia,si_bias);
  NOM=fopen(nomAr,"r");



  int i,j;


  double chi2min;
  int min;

//	buscando el minimo, es necesario antes de calcular la integra porque normaliza los dados
//	de otro nodo casi todo contribuye con 0, y el pico es 10^-un_monton
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cos[i]);
    if(i==0){ chi2min=chi2[0]; min=0;}
    if(chi2min>chi2[i]){	chi2min=chi2[i];	min=i;}
  }
  rewind(NOM);

//	ahora si estimamos la integral de la distribucion
  for(i=0;i<NP;i++){
    fscanf(NOM,"%le %le\n",&chi2[i],&param_cos[i]);
    if(i==0)	integral[i]=exp(-0.5*(chi2[i]-chi2min));
    else	integral[i]=integral[i-1]+exp(-0.5*(chi2[i]-chi2min));
  }


  fclose(NOM);
if(obs==0) printf("\n");
printf("chimin=%lf param=%le id=%d ",chi2min/dof[obs],param_cos[min],min);
chi2_dof[(MG-1+cruzado)%2][first_theory]=chi2min/dof[obs]*2;	//	saving chi2
//	computing delta_sigma
if(obs==0)      delta_cr[(MG-1+cruzado)%2][first_theory]=param_cos[min];
if(obs==1)	delta_cr[(MG-1+cruzado)%2][first_theory]-=param_cos[min];

  int med=0;
  while(integral[med]<0.5*integral[NP-1])
    med++;
    printf("%c%c %le %le ",NameTheory[MG-1],NameTheories[6-first_theory],param_cos[min],param_cos[med]);
    fflush(stdout);
    fprintf(TEX,"%.2lf & %.2lf ",ar(param_cos[min]),ar(param_cos[med]));
    fflush(stdout);


for(j=0;j<3;j++){
  int izq=0;
  int der=NP-1;
  double fraccion=integral[der];
  if((first_theory==3)){
    while(integral[der]>porcentage[j]*integral[NP-1]){
      der--;
    }
    printf("%le ",param_cos[der]);
    fprintf(TEX,"& <%.2lf ",ar(param_cos[der]));
    //	max of the two sigma
    if((j==0)&&(obs==0))
      sigma_cr[(MG-1+cruzado)%2][first_theory]=ar(param_cos[der]);
    if((j==0)&&(obs==1))
      if(sigma_cr[(MG-1+cruzado)%2][first_theory]<ar(param_cos[der]))
        sigma_cr[(MG-1+cruzado)%2][first_theory]=ar(param_cos[der]);
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
    //	max of the two sigma
    if((j==0)&&(obs==0))
      sigma_cr[(MG-1+cruzado)%2][first_theory]=ar(0.5*(-param_cos[izq]+param_cos[der]));
    if((j==0)&&(obs==1))
      if(sigma_cr[(MG-1+cruzado)%2][first_theory]<ar(0.5*(-param_cos[izq]+param_cos[der])))
        sigma_cr[(MG-1+cruzado)%2][first_theory]=ar(0.5*(-param_cos[izq]+param_cos[der]));
  }
  if((j==0)&&(obs==1)){double auxx=delta_cr[(MG-1+cruzado)%2][first_theory]/sigma_cr[(MG-1+cruzado)%2][first_theory];
           delta_cr[(MG-1+cruzado)%2][first_theory]=abs(auxx);}

}
printf("\n");
fprintf(TEX,"\\\\ \n");


if(obs!=0)
fprintf(GNU,"re");
fprintf(GNU,"plot '%s' u 2:(exp(-0.5*($1-%le))) w l dt %d lc %d lw %d t '%s'\n",nomAr,chi2[min],dash[obs],color[obs],grueso[obs],nomOb);
fprintf(TOD,"%le\n",param_cos[min]);

}
//////////////////////////////////////////////////////////////
	//	end observables
//////////////////////////////////////////////////////////////






fprintf(GNU,"set xr[%le:%le]\n",ginf,gsup);
if(cruzado==1)
fprintf(GNU,"set output 'histo_%c%c.png'\n",NameTheory[MG-1],NameTheories[6-first_theory]);
else
fprintf(GNU,"set output 'histo_%c%c.png'\n",NameTheory[MG-1],NameTheories[first_theory]);
fprintf(GNU,"replot\n");
fprintf(GNU,"reset\n\n\n");
fprintf(TEX,"\n");
fclose(TEX);
}	//	teoria
printf("\n");
}	//	clase de teoria
fclose(GNU);
fclose(TOD);


}	//	cruzado


  sprintf(nomAr,"chi2.tex");
  CHX=fopen(nomAr,"w+");
  fprintf(CHX,"$|f_{R0}|=10^{-4}$ & \\textbf{%.2lf} & %.2lf \\\\ \n",chi2_dof[0][0],chi2_dof[0][6]);	printf("%le\n",chi2_dof[0][6]/chi2_dof[0][0]);
  fprintf(CHX,"$|f_{R0}|=10^{-5}$ & \\textbf{%.2lf} & %.2lf \\\\ \n",chi2_dof[0][1],chi2_dof[0][5]);	printf("%le\n",chi2_dof[0][5]/chi2_dof[0][1]);
  fprintf(CHX,"$|f_{R0}|=10^{-6}$ & \\textbf{%.2lf} & %.2lf \\\\ \n",chi2_dof[0][2],chi2_dof[0][4]);	printf("%le\n",chi2_dof[0][4]/chi2_dof[0][2]);
  fprintf(CHX,"$\\Lambda$CDM & %.2lf & %.2lf \\\\ \n",chi2_dof[1][3],chi2_dof[0][3]);
  fprintf(CHX,"$z_{SSB}=1$ & %.2lf & \\textbf{%.2lf} \\\\ \n",chi2_dof[1][2],chi2_dof[1][4]);	printf("%le\n",chi2_dof[1][2]/chi2_dof[1][4]);
  fprintf(CHX,"$z_{SSB}=2$ & %.2lf & \\textbf{%.2lf} \\\\ \n",chi2_dof[1][1],chi2_dof[1][5]);	printf("%le\n",chi2_dof[1][1]/chi2_dof[1][5]);
  fprintf(CHX,"$z_{SSB}=3$ & %.2lf & \\textbf{%.2lf} \\\\ \n",chi2_dof[1][0],chi2_dof[1][6]);	printf("%le\n",chi2_dof[1][0]/chi2_dof[1][6]);
  fclose(CHX);

  sprintf(nomAr,"delta.tex");
  CHX=fopen(nomAr,"w+");
  fprintf(CHX,"$|f_{R0}|=10^{-4}$ & \\textbf{%.2lf} & %.2lf \\\\ \n",delta_cr[0][0],delta_cr[0][6]);
  fprintf(CHX,"$|f_{R0}|=10^{-5}$ & \\textbf{%.2lf} & %.2lf \\\\ \n",delta_cr[0][1],delta_cr[0][5]);
  fprintf(CHX,"$|f_{R0}|=10^{-6}$ & \\textbf{%.2lf} & %.2lf \\\\ \n",delta_cr[0][2],delta_cr[0][4]);
  fprintf(CHX,"$\\Lambda$CDM & %.2lf & %.2lf \\\\ \n",delta_cr[1][3],delta_cr[0][3]);
  fprintf(CHX,"$z_{SSB}=1$ & %.2lf & \\textbf{%.2lf} \\\\ \n",delta_cr[1][2],delta_cr[1][4]);
  fprintf(CHX,"$z_{SSB}=2$ & %.2lf & \\textbf{%.2lf} \\\\ \n",delta_cr[1][1],delta_cr[1][5]);
  fprintf(CHX,"$z_{SSB}=3$ & %.2lf & \\textbf{%.2lf} \\\\ \n",delta_cr[1][0],delta_cr[1][6]);
  fclose(CHX);



  cout<<"terminó"<<endl;
  return 0;
}

