#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../l/parametros.h"

//	split the halo catalog into high and low concentration by using the r200 to rx00 ratio as the concentration parameter



int main(){
double min_mass_low[7]={30000.};
double min_mass_high[7]={30000.};
double max_mass_low[7]={0.};
double max_mass_high[7]={0.};


for(int t=0;t<7;t++){
  char catalogo_halos[500];
  char catalogo_halos_low[500];
  char catalogo_halos_high[500];
  sprintf(catalogo_halos,"../%c/halos_m_%s%s",NameTheories[t],prefijo,caso);	//	halos with mass estimation for different thresholds
  sprintf(catalogo_halos_low,"../%c/halos_lowc_%s%s",NameTheories[t],prefijo,caso);
  sprintf(catalogo_halos_high,"../%c/halos_highc_%s%s",NameTheories[t],prefijo,caso);
  FILE *IN;
  FILE *LO;
  FILE *HI;
  IN=fopen(catalogo_halos,"r");
  LO=fopen(catalogo_halos_low,"w+");
  HI=fopen(catalogo_halos_high,"w+");

  double x,y,z,m[9],dm[9],r[9],dr[9];
  int NumHalos=0;
  double c,mean=0.;
  while(fscanf(IN,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
&x,&y,&z,&m[0],&dm[0],&r[0],&dr[0],&m[1],&dm[1],&r[1],&dr[1],&m[2],&dm[2],&r[2],&dr[2],&m[3],&dm[3],&r[3],&dr[3],&m[4],&dm[4],&r[4],&dr[4],&m[5],&dm[5],&r[5],&dr[5],&m[6],&dm[6],&r[6],&dr[6],&m[7],&dm[7],&r[7],&dr[7],&m[8],&dm[8],&r[8],&dr[8])!=EOF){
    c=r[2]/r[0];
    if(c<1.){
      mean+=c;
      NumHalos++;
    }  
  }
  rewind(IN);
  mean/=NumHalos;
  printf("mean=%le\n",mean);
  while(fscanf(IN,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
&x,&y,&z,&m[0],&dm[0],&r[0],&dr[0],&m[1],&dm[1],&r[1],&dr[1],&m[2],&dm[2],&r[2],&dr[2],&m[3],&dm[3],&r[3],&dr[3],&m[4],&dm[4],&r[4],&dr[4],&m[5],&dm[5],&r[5],&dr[5],&m[6],&dm[6],&r[6],&dr[6],&m[7],&dm[7],&r[7],&dr[7],&m[8],&dm[8],&r[8],&dr[8])!=EOF){
    c=r[2]/r[0];
    if(c<1){
      if(c<mean){
        fprintf(HI,"%lf %lf %lf %lf %lf\n",x,y,z,m[0],r[0]);
        if(min_mass_high[t]>m[0])	min_mass_high[t]=m[0];
        if(max_mass_high[t]<m[0])	max_mass_high[t]=m[0];
      }
      else{
        fprintf(LO,"%lf %lf %lf %lf %lf\n",x,y,z,m[0],r[0]);
        if(min_mass_low[t]>m[0])	min_mass_low[t]=m[0];
        if(max_mass_low[t]<m[0])	max_mass_low[t]=m[0];
      }
    }
  }

  fclose(IN);
  fclose(LO);
  fclose(HI);
}

double min_mass=0.;
double max_mass=30000.;
for(int t=0;t<7;t++){
  if(min_mass<min_mass_low[t])	min_mass=min_mass_low[t];
  if(min_mass<min_mass_high[t])	min_mass=min_mass_high[t];
  if(max_mass>max_mass_low[t])	max_mass=max_mass_low[t];
  if(max_mass>max_mass_high[t])	max_mass=max_mass_high[t];
}
printf("sup mass = %le\n",min_mass);
printf("inf mass = %le\n",max_mass);



// solo falta imprimir la lista de halos, etc, como lo habia hecho en FINDER
// para poder usar esta salida dentro de los demas codigos

  return 0;
}


