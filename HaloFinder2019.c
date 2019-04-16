#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../l/parametros.h"

const double rho_cr=2.77526627;		//	[10¹¹ h² Ms / Mpc³]
const int   nc=np;
const int   NumPart=np*np*np;		//	numero total de particulas
const int   NumCel=nc*nc*nc;		//	numero total de particulas
const double Mp=pow(Lbox/np,3.)*Om*rho_cr;	//	masa de 1 particula [10¹¹ Ms/h]
const double Mmin=40*Mp;			//	masa minima de halo [10¹¹ Ms/h] 40 part/halo_min
const double pi=4.*atan(1.);		//	numero pi
const int   NumDelta=9;
const double delta[NumDelta]={200.*Om*rho_cr,300.*Om*rho_cr,400.*Om*rho_cr,600.*Om*rho_cr,800.*Om*rho_cr,1200.*Om*rho_cr,1800.*Om*rho_cr,2400.*Om*rho_cr,3200.*Om*rho_cr};

//const double limite=pow(1.*np/Lbox,3.)*Mmin/(Om*rho_cr);//	densidad limite para ser centro de halo
const int  bin=10000;			//	numero de bines radiales log

int *Ocu;			//	numero de ocupacion de cada celda
int *Rho;			//	densidad de las particulas
int *Masa;			//	masa de las particulas (=0 si \in halo)
int *lista;			//	lista de Id de las particulas dentro de Rmax
int **Id;			//	indice de las particulas en cada celda
double *X,*Y,*Z;			//	coordenadas de todas las particulas





int calc_den=1;			//	1: si calcula	0: lee del archivo previo




void allocar(void){
  if(!(Ocu = (int*) calloc (NumCel , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Rho = (int*) calloc (NumPart , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Masa = (int*) malloc (NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(lista = (int*) malloc (NumPart * sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(X = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Y = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Z = (double*) malloc (NumPart * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      fflush(stderr);
      exit(0);
    }
  printf("alloca...hecho\n");
  fflush(stdout);
}






void allocaId(void){
  printf("alloca Id...\n");
  fflush(stdout);
  Id = (int **) malloc(NumCel * sizeof(int *));
  for(int p = 0; p < NumCel; p++) { 
   if(Ocu[p]>0)
    Id[p] = (int *) malloc(Ocu[p] * sizeof(int));
   else
    Id[p] = (int *) malloc( sizeof(int));
   }
  printf("alloca Id...hecho\n");
  fflush(stdout);
}










int main(int argc, char **argv)
{
  printf("comienza\n");
  fflush(stdout);
  omp_set_num_threads(16);

  double x,y,z;
  double fx,fy,fz;
  int dx,dy,dz,dxx,dyy,dzz;
  int i,j,k,ijk,l,p;

  allocar();

  char catalogo_materia[500];
  char catalogo_halos[500];
  char particulas_halos[500];
  char archivo_densidad[500];
  sprintf(catalogo_materia,"%s%s",prefijo,caso);
  sprintf(catalogo_halos,"halos_m_%s%s",prefijo,caso);
  sprintf(archivo_densidad,"densidad_%s%s",prefijo,caso);

  printf("lee xyz...\n");
  fflush(stdout);
  FILE * LEE;
  LEE = fopen(catalogo_materia,"r");	//halo.txt
  // crea el vector de posiciones y cuenta las particulas que hay en cada celda
  for (p=0;p<NumPart;p++)
    {
    fscanf(LEE,"%lf %lf %lf %lf %lf %lf\n",&x,&y,&z,&fx,&fy,&fz);
    X[p]=x*nc/Lbox;
    Y[p]=y*nc/Lbox;
    Z[p]=z*nc/Lbox;
    i=floor(X[p]+0.5);	//	tengo que sumar 0.5 porque las células de densidad
    j=floor(Y[p]+0.5);	//	están centradas en los vértices enteros
    k=floor(Z[p]+0.5);	//	no en puntos tales como 0.5, 15.5, 364.5
    ijk=i%nc+(j%nc)*nc+(k%nc)*nc*nc;
    Ocu[ijk]++;
    }
  fclose (LEE);
  printf("lee xyz...hecho\n");
  fflush(stdout);

  allocaId();

  for(p=0;p<NumCel;p++)
    Ocu[p]=0;

  for(p = 0; p < NumPart; p++) { 
    i=floor(X[p]+0.5);	//	tengo que sumar 0.5 porque las células de densidad
    j=floor(Y[p]+0.5);	//	están centradas en los vértices enteros
    k=floor(Z[p]+0.5);	//	no en puntos tales como 0.5, 15.5, 364.5
    ijk=i%nc+(j%nc)*nc+(k%nc)*nc*nc;
    Id[ijk][Ocu[ijk]]=p;
    Ocu[ijk]++;
  }
  printf("llena Id...hecho\n");
  fflush(stdout);




  // Calcula densidades
  double r,r2;					//	radio para buscar vecinos
  r=pow(3*Mmin/(4*pi*delta[0]),1./3)*nc/Lbox;	//	en unidades de nc
  r2=r*r;
  printf("radio mínimo = %lf\n",r*Lbox/nc);
  fflush(stdout);
  double dist;


  double aux_i,aux_j,aux_k,aux_aux;

if(calc_den==1){
  for(p=0;p<NumPart;p++){
    x=X[p];
    y=Y[p];
    z=Z[p];
    double aux_r=r+sqrt(2.);
    dx=floor(x-aux_r+0.5);
    dy=floor(y-aux_r+0.5);
    dz=floor(z-aux_r+0.5);
    dxx=floor(x+aux_r+0.5);
    dyy=floor(y+aux_r+0.5);
    dzz=floor(z+aux_r+0.5);
    aux_aux=r2+3./4+r*sqrt(3.);
    for(k=dz;k<=dzz;k++){
      if(k<0) fz=-1.*nc;
      else if(k>=nc) fz=1.*nc;
      else fz=0.;
      for(j=dy;j<=dyy;j++){
        if(j<0) fy=-1.*nc;
        else if(j>=nc) fy=1.*nc;
        else fy=0.;
        for(i=dx;i<=dxx;i++){
          if(i<0) fx=-1.*nc;
          else if(i>=nc) fx=1.*nc;
          else fx=0.;

	  ijk=(i+nc)%nc+((j+nc)%nc)*nc+((k+nc)%nc)*nc*nc;
	  if(Ocu[ijk]>0){

            dist=pow(x-1.*i,2.)+pow(y-1.*j,2.)+pow(z-1.*k,2.);

            if(dist<aux_aux){	//	célula que tiene partículas dentro del rádio requerido

	      aux_i=x-fx;
	      aux_j=y-fy;
	      aux_k=z-fz;
	      int id_aux;
	      #pragma omp parallel
	      {
	      double dist_priv;
	      int count_priv=0;
	      #pragma omp for
	      for(l=0;l<Ocu[ijk];l++){
	        id_aux=Id[ijk][l];
	        dist_priv=pow(aux_i-X[id_aux],2.)+pow(aux_j-Y[id_aux],2.)+pow(aux_k-Z[id_aux],2.);
	        if (dist_priv<r2)
	           count_priv++;}
	      #pragma omp critical
	      {
	      Rho[p]+=count_priv;
	      }
	      }
	    }	//	cierra el if de la distancia
	  }	//	cierra el idel número de ocupación
	}
      }
    }
  }		//	cierra el ciclo sobre las partículas
  FILE * DEN;
  DEN = fopen(archivo_densidad,"w+");
  for(p=0;p<NumPart;p++)
    fprintf(DEN,"%d\n",Rho[p]);
  fflush(DEN);
  fclose(DEN);
  printf("calculó la densidad\n");
  fflush(stdout);
}		//	cierra el if




else{
  FILE * DEN;
  DEN = fopen(archivo_densidad,"r");
  for(p=0;p<NumPart;p++)
    fscanf(DEN,"%d\n",&Rho[p]);
  fclose(DEN);
  printf("leyó la densidad\n");
}
















  int cont[bin];
  double vol[bin];			//	Mp/vol de cada bin
  double rad[bin];			//	radio de cada bin
  double Rmax=4.5*nc/Lbox;		//	radio maximo para formar halos
  double Rmax2=Rmax*Rmax;
  double c=pow(Rmax/r,1./bin);		//	constante de proporcionalidad log
  double aux1=0.5/log(c),aux2;
  int PartHalo;				//	numero de particulas dentro de Rmax
  int aux_int;

  // volumen y radio de cada bin
//  #pragma omp parallel for
  for(i=0;i<bin;i++){
    vol[i]=Mp*3/(4*pi*pow(r*Lbox/nc,3.)*pow(c,3.*i));	//	en [h/Mpc]^3
    rad[i]=r*Lbox/nc*pow(c,1.*i);}			//	en Mpc/h
  printf("creó radio y volumen inverso de cada bin\n");
  fflush(stdout);

  FILE * ESC;
  ESC = fopen(catalogo_halos,"w+");
  FILE * WRI;
  WRI = fopen(particulas_halos,"w+");

  int centro=1;
  double rho[bin];
  rho[0]=1.01*delta[0];			//	para que entre al while
  printf("entra al while\n");
  fflush(stdout);

  while(centro>=0){

    // Busco la particula mas densa
    int rho_max=ceil(Mmin/Mp);
    centro=-1;
    #pragma omp parallel
    {
    double rho_max_priv=rho_max;
    int centro_priv=-1;
    #pragma omp for
    for(i=0;i<NumPart;i++){
      if((rho_max_priv<Rho[i])&&(Rho[i]>=0)){
        rho_max_priv=Rho[i];
        centro_priv=i;}}
    #pragma omp critical
    {
      if(rho_max<rho_max_priv){
        rho_max=rho_max_priv;
        centro=centro_priv;}
    }
    }

    // Densidad física de la partícula con más vecinos
    rho[0]=vol[0]*Rho[centro];

    if(rho[0]<delta[0]){
      centro=-1;
      printf("en este punto deberia salir del while\n");}
    else{
      #pragma omp parallel for
      for(j=0;j<bin;j++)
        cont[j]=0;
      double dista;
      double aux_Rmax=Rmax+sqrt(2.);

      PartHalo=0;
      x=X[centro];
      y=Y[centro];
      z=Z[centro];
      dx=floor(x-aux_Rmax+0.5);
      dy=floor(y-aux_Rmax+0.5);
      dz=floor(z-aux_Rmax+0.5);
      dxx=floor(x+aux_Rmax+0.5);
      dyy=floor(y+aux_Rmax+0.5);
      dzz=floor(z+aux_Rmax+0.5);
 
      aux_aux=Rmax2+3./4+Rmax*sqrt(3.);

      for(k=dz;k<=dzz;k++){
        if(k<0) fz=-1.*nc;
        else if(k>=nc) fz=1.*nc;
        else fz=0.;
        for(j=dy;j<=dyy;j++){
          if(j<0) fy=-1*nc;
          else if(j>=nc) fy=1.*nc;
          else fy=0.;
          for(i=dx;i<=dxx;i++){
            if(i<0) fx=-1.*nc;
            else if(i>=nc) fx=1.*nc;
            else fx=0.;

            ijk=(i+nc)%nc+((j+nc)%nc)*nc+((k+nc)%nc)*nc*nc;
	    if(Ocu[ijk]>0){

              dist=pow(x-1.*i,2.)+pow(y-1.*j,2.)+pow(z-1.*k,2.);
    	      if(dist<aux_aux){

	        aux_i=x-fx;
	        aux_j=y-fy;
	        aux_k=z-fz;

                for(l=0;l<Ocu[ijk];l++){
                  int aux_id=Id[ijk][l];
                  if(Rho[aux_id]!=-1){
                    dista=pow(aux_i-X[aux_id],2.)+pow(aux_j-Y[aux_id],2.)+pow(aux_k-Z[aux_id],2.);
                    if (dista<Rmax2){
                      if(dista<=r2)
                        aux2=0.1;
                      else
                        aux2=aux1*log(dista/r2);
                      aux_int=floor(aux2);	//	es el bin radial de esta particula
                      cont[aux_int]+=1;
                      Masa[aux_id]=aux_int;
                      lista[PartHalo]=aux_id;
                      PartHalo++;		//	masa = bin radial
	  	    } 
		  }
	        }
	      }		//	cierra el if de la distancia
	    }		//	cierra el if del número de ocupación
	  }
	}
      }			//	cierra los tres for

      // número de partículas en el rádio de búsqueda
      PartHalo--;

      // cuenta las particulas en cada bin
      // y decide cual es el radio (el indice del bin) que tienen densidad adecuada

      //  densidad del bin radial más pequeño
      rho[0]=vol[0]*cont[0];

      if(rho[0]>delta[0]){

        //	a partir de los conteos en las capas calculo la densidad media
        for(l=1;l<bin;l++){
          cont[l]+=cont[l-1];
          rho[l]=vol[l]*cont[l];}
        if(rho[bin-1]>delta[0]){
          printf("error, faltaron bines radiales para cubrir un halo\n");
          printf("%f\n",rho[0]);
          fflush(stdout);
          exit(0);}

        //imprimo centro
        fprintf(ESC,"%le %le %le ",X[centro]*Lbox/nc,Y[centro]*Lbox/nc,Z[centro]*Lbox/nc);
        int lmin,Lmax,l0;
        //imprimo radios y masas segun la lista de Delta
        for(l=0;l<NumDelta;l++){
          lmin=0,Lmax=bin-1;
          if((rho[lmin]>delta[l])&&(rho[Lmax]<delta[l])){
            while(rho[lmin]>delta[l]) lmin++;
            while(rho[Lmax]<delta[l]) Lmax--;}
          if(l==0)l0=lmin;
          fprintf(ESC,"%le %le %le %le ",(cont[lmin]+cont[Lmax])*Mp*0.5,abs(cont[lmin]-cont[Lmax])*Mp*0.5,(rad[lmin]+rad[Lmax])*0.5,abs(rad[lmin]-rad[Lmax])*0.5);}
        fprintf(ESC,"\n");fflush(ESC);


        // ahora que tengo el radio del halo debo anular la masa de sus particulas
        // para que no contribuyan a los demas halos
        for(i=0;i<PartHalo;i++){
          if(Masa[lista[i]]<=l0){	//Masa[lista[l]] is the radial bin where the particle [i] is, whilw lmin[0] is the radial bin defining the halo200
            Rho[lista[i]]=-1;
	   }}
      }
      else{
        Rho[centro]=-2;}

      } // if rho>=delta
    }   // while centro>=0
  fclose(ESC);
  fclose(WRI);



// solo falta imprimir la lista de halos, etc, como lo habia hecho en FINDER
// para poder usar esta salida dentro de los demas codigos

  return 0;
}


