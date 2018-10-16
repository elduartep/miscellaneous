//Calcula el numero de voids a partir de las posiciones
//usando spherical-overdensity>200
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../l/parametros.h"



const int   acumulada=0;
const int   bin=15;                    //	numero de bines radiales log

const float rho_cr=2.77526627;          //	[10¹¹ h² Ms / Mpc³]
const float pi=4.*atan(1.);             //	numero pi


const int esferas=1;
const int JN=8;

int main(int argc, char **argv){
for(int t=0;t<=6;t++){


  char inFile[560];
  char outFile[560];
  char outFile2[560];
  if(esferas==1){
    sprintf(inFile,"../%c/esferas_%s%s",NameTheories[t],prefijo,caso);
    sprintf(outFile,"../%c/dnv_esferas_%s%s",NameTheories[t],prefijo,caso);
    sprintf(outFile2,"../%c/nv_esferas_%s%s",NameTheories[t],prefijo,caso);}
  else{
    sprintf(inFile,"union_%s%s",prefijo,caso);
    sprintf(outFile,"dnv_union_%s%s",prefijo,caso);
    sprintf(outFile2,"nv_union_%s%s",prefijo,caso);}



  int i, j, k, l;
  int NumVoids=0;
  int id;

  float x,y,z,radio;
  float min,max;

  float nv[JN+1][bin];
  float snv[JN+1][bin];
  float Rad[JN+1][bin];
  float sRad[JN+1][bin];
  float sigma_nv[bin];

  for(i=0;i<bin;i++){
    sigma_nv[i]=0.;
    for(j=0;j<=JN;j++){
      nv[j][i]=0.;
      snv[j][i]=0.;
      Rad[j][i]=0.;
      sRad[j][i]=0.;
    }
  }


///////////////////////////////////		masa 	maxima y minima

  printf("Leyendo el catalogo para establecer radio max y min\n");
  FILE * L;
  FILE * E;
  FILE * F;
  L=fopen(inFile,"r");
  max=0.;
  min=Lbox;
  while(fscanf(L,"%f %f %f %f\n",&x, &y, &z, &radio)!=EOF){
    if(radio<min)
      min=radio;
    if(radio>max)
      max=radio;
    NumVoids++;}
  fclose(L);
  printf("Se encontraron %d voids en el catálogo\n",NumVoids);

  printf("radio_max=%le radio_min=%le\n",max,min);


///////////////////////		definicion de los bines radiales
  float Rmax =max*1.01;
  float r    =min*0.99;

//////////////////////  la misma escala para todas las teorias
  r=radio_min_todos*0.99;
  Rmax=radio_max_todos*1.01;

  float c    =pow(Rmax/r,1./bin);	//	constante de proporcionalidad log
  float aux1 =1./log(c);
  float aux2,aux3;
  int id_bin,id_JN,id_i,id_j,id_k;

//////////////////////		lleno los bines, contador y radio
  L=fopen(inFile,"r");
  for(i=0;i<NumVoids;i++){
    fscanf(L,"%f %f %f %f\n",&x, &y, &z, &radio);
  if((radio<Rmax)&&(radio>r)){
    aux2=aux1*log(radio/r);
    id=floor(aux2);		//	es el bin radial de este void

    if(x<Lbox/2)	id_i=0;
    else		id_i=1;
    if(y<Lbox/2)	id_j=0;
    else		id_j=1;
    if(z<Lbox/2)	id_k=0;
    else		id_k=1;
    id_JN=id_k+2*id_j+4*id_i;

    if(l==JN) l=0;
    for(j=0;j<JN;j++){		//	submuestras
//      if(j!=id_JN){						//	todos las submuestras menos l
      if(j==id_JN){						//	solo la submuestra l
        nv[j][id]+=1.;		//	contador entero del número de voids
        Rad[j][id]+=radio;
      }
    }
    nv[JN][id]+=1.;			//	total
    Rad[JN][id]+=radio; 
  }
  }
  fclose(L);

  for(i=0;i<bin;i++)
    printf("%le ",Rad[JN][i]/nv[JN][i]);
  printf("\n");

	//	como nv es una densidad, y para j\in[0,JN) solo tube en cuenta una fraccion de
	//	objetos igual a float(JN-1)/JN, entonces debo multiplicarlos por JN/float(JN-1)
//  float cor=8.;//float(JN)/(JN-1);





	//	NO acumulada

    for(i=0;i<bin;i++){
      for(j=0;j<=JN;j++){
        if(nv[j][i]>0)
          Rad[j][i]/=nv[j][i];		//	media de radios en el bin i
//        nv[j][i]*=Rad[j][i]/(r*pow(c,i*1.)*(c-1.));	//funcion de radio
      }
    }
  	//	errores abundancia: jacknife
    for(i=0;i<bin;i++){		//	bines
      for(j=0;j<JN;j++){
        sigma_nv[i]+=pow(nv[j][i]*JN-nv[JN][i],2);
      }
    }
    for(i=0;i<bin;i++)
      sigma_nv[i]*=1./float(JN*JN);

    //	imprine abundancia
    E=fopen(outFile,"w+");  
    for(i=bin-1;i>=0;i--){

    float factor=pow(Lbox,-3.)/log(c);
    double error=nv[JN][i];
    if(error<sigma_nv[i])	error=sigma_nv[i];
//      fprintf(E,"%f %e %e %e\n",
//      r*pow(c,1.*i+0.5),nv[JN][i]*factor,sqrt(nv[JN][i])*factor,sqrt(sigma_nv[i])*factor);
      fprintf(E,"%f %e %e\n",
      r*pow(c,1.*i+0.5),nv[JN][i]*factor,sqrt(nv[JN][i])*factor);
    }

    fclose(E);
  


}
return 0;
}

