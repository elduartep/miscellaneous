#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parametros.h"
#include "nr3.h"
#include "ran.h"


//	11 de Abril de 2018
//	versión anterior: VoidFinderPrueba2.c
//	voy a eliminar el sobrelapamiento
//	creo que así, la cantidad de voids pequeños va a disminuir
//	la idea es poder modelar la abundancia usando la teoria esferica
//	el problema en *2.c es que parece que la cantidad de voids aumenta potencialmente cuando disminuye el tamaño

//	actualmente estoy trabajando solo con los voids esfericos
//	así que no voy a entrar en la parte de enriquecimiento


//	la primera modificación va a estar en la función: criterio(...)
//	que es la que me diría si hay sobrelapamiento con otro void o nó








//	no se olvide que las células están centradas en los vértices enteros del sistema cartesiano
//	Nó están en puntos tales como 0.5, 1.5, 203.5

//	Todos los arreglos son globales, así que puedo modificarlos desde cualquier lugar

//	radios en unidades de nc
//	volumen y densidade en unidades físicas


const int   calcula=0;			//	=1 calcula densidad celular, =0 lee el archivo existente


  char catalogo_voids[500];		//	voids==conjunto de huecos que se sobrelapan a la anterior
  char catalogo_particulas[500];	//	particulas de la simulacion
  char densidad_radios[500];		//	radio de las esferas centradas en la grilla cartesiana



const int   nc=np;			//	munero 1D de celulas, mi ~KDtree
const int   nh=np;			//	numero 1D de pixeles para huecos, grilla cartesiana
const int   NumPart=np*np*np;		//	numero total de particulas
const int   NumCel=nc*nc*nc;		//	numero total de celulas para el ~KDtree
const int   NumVoid=10000000;		//	numero maximo de voids que esperamos encontrar
const int   nc2=nc*nc;
const int   nc24=nc2/4;

const float rho_cr=2.77526627;		//	[10¹¹ h² Ms / Mpc³]
const float Mp=pow(Lbox,3.)*Om*rho_cr/NumPart;	//	masa de 1 particula [10¹¹ Ms/h]
const float pi=4.*atan(1.);		//	numero pi
const float limite=0.2*Om*rho_cr;	//	sobredensidad de los voids en [10¹¹ Ms h² / Mpc³]
const float halo=200.*Om*rho_cr;	//	sobredensidad de los halo en [10¹¹ Ms h² / Mpc³]

//	en http://arxiv.org/pdf/1403.5499v2.pdf solo toman en cuenta voids con radio > 2 rp
//	donde rp es la distancia media entre particulas
//	y a la hora de calcular el perfil de densidad solamente muestran el resultado para r > rp
const float rad_min_dos=0.2;		//	radio minimo HD en unidades de nc, centrados en randoms
const float rad_min_uno=2.;		//	radio minimo pixelizado en unidades de nc, red cartesiana

//	para el crecimiento de las esferas centradas en la red
//  radio maximo de busqueda en unidades de nc, equivale a 32Mpc     /4=8Mpc
const int   incremento=(int)(8.*nc/Lbox);//	cuando no es suficiente el radio siguiente
const int   otro=4*incremento;		//	radio de busqueda en unidades de nc
const int   semi=2;			//	=2, entonces radio semi-entero
const int   bin=semi*otro;		//	numero de bines radiales log, es el maximo
const int   semi_otro_max=768;		//	es el maximo semi*otro para los nc(128) y Lbox(256) que estamos usando


//	para refinar los radios estimados en centros de red y centrados en randoms
const float despl_max=0.05;		//	paso maximo al buscar el centro, fracción del radio
const float despl_min=0.01;		//	paso minimo = error del centro, fracción del radio

int   *Ocu;			//	numero de ocupacion de cada celda
int  **Id;			//	indice de las particulas en cada celda

float *Rho;			//	densidad de las celdas
float *Rad;			//	radio de la máxima esfera subdensa centrada en cada celula

float *X,*Y,*Z;			//	coordenadas de todas las particulas

float *Radv;			//	radio de los voids
float *Xv,*Yv,*Zv;		//	coordenadas del centro de los voids

int   *ID_random;		//	<0:random >0:pixel
int   *IDV;			//	se incrementa cada vez que hallo un nuevo void
int    NumHuecos;		//	candidad de esferas subdensas

float *cont;	//	masa acumulada de cada capa de células, ref: partículas en cada bin
float *vol;	//	Mp/vol de cada bin

int NumVoids;					//	número de familias: volumen y elip



float rad_rad[9];	//	variables para comparar con los vecinos
float rho_rho[9];
float rho_rho_fail[9];
float x_x[9];
float y_y[9];
float z_z[9];



int   id_void=-1;		//	número de voids encontrados
float rad_min_actual;		//	es el radio minimo de cada fase: 2 o 0.5

float radio_maximo_para_refinar;




void archivos(void){
  sprintf(catalogo_voids,  	    "voids_esfericos_%s%s",prefijo,caso); // voids esfericos con densidad 0.2
  sprintf(catalogo_particulas,	                    "%s%s",prefijo,caso); // particulas de la simulacion
  sprintf(densidad_radios, 	             "radios_%s%s",prefijo,caso); // radio de las esferas centradas en la grilla cartesiana

}







//	el número máximo de voids que espero encontrar es fijo
//	en este momento ya tengo los huecos de la primera fase
//	y voy a seguir encontrando más
//	por lo tanto no necesito contar cuantos llevo hasta ahora

//	lee el archivo de respaldo (el que no tiene HD)



void allocar(void){
  printf("\nallocar...\n");fflush(stdout);
  if(!(Rho = (float*) calloc (NumCel , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Rho.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Ocu = (int*) calloc (NumCel , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory Ocu.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Rad = (float*) calloc (NumCel , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Rad.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(X = (float*) malloc (NumPart * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory X.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Y = (float*) malloc (NumPart * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Y.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Z = (float*) malloc (NumPart * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Z.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Radv = (float*) calloc (NumVoid , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Radv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Xv = (float*) malloc (NumVoid * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Xv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yv = (float*) malloc (NumVoid * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Yv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zv = (float*) malloc (NumVoid * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Zv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(cont = (float*) malloc (bin * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory cont.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(vol = (float*) malloc (bin * sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory vol.\n");
      fflush(stderr);
      exit(0);
    }
  printf("allocar...hecho\n");
  fflush(stdout);
}








void allocaId(void){
  printf("\nalloca Id...\n");
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






void free_pixelizado(void){
  printf("\nfree_pixelizado...\n");fflush(stdout);
  free(Rho);			//	densidad CIC de las células
  free(Rad);			//	radio de la esfera subdensa más grande centrada en cada célula
  printf("free_pixelizado...hecho\n");fflush(stdout);
}











void anula(float x,float y,float z,float radio){
  // Dados x, y, z y rad_aux, anula Rad[ijk] para todas las células ijk in void (+ rad_min_actual)

  int i,j,k,ii,jj,kk,iii,jjj,kkk;
  ii=floor(x-radio+0.5-rad_min_actual);
  jj=floor(y-radio+0.5-rad_min_actual);
  kk=floor(z-radio+0.5-rad_min_actual);
  iii=floor(x+radio+0.5+rad_min_actual);
  jjj=floor(y+radio+0.5+rad_min_actual);
  kkk=floor(z+radio+0.5+rad_min_actual);
  float dy,dz;

  for(k=kk;k<=kkk;k++){
    dz=z-1.*k;
    for(j=jj;j<=jjj;j++){
      dy=y-1.*j;
//      #pragma omp parallel
      {
      float dist_priv,dx_priv;
      int ijk_priv;
//      #pragma omp for
      for(i=ii;i<=iii;i++){
        dx_priv=x-1.*i;
	dist_priv=dx_priv*dx_priv+dy*dy+dz*dz;
	if(dist_priv<(radio*radio)){		//	deja una franja sin anular cerca al borde?
	  ijk_priv=((i+nc)%nc)+nc*((j+nc)%nc)+nc2*((k+nc)%nc);
	  Rad[ijk_priv]=-1.;
	  Rho[ijk_priv]=halo;
	}
      }
      }
    }
  }

}














//	dada una posición x,y,z
//	me dice cual es la distancia al vorde del void más cercano
void criterio(float x, float y, float z){

  float radio_maximo=nc*10.;	//	es la variable que voy a retornar
  int aux_void;

  #pragma omp parallel
  {
  float tx,ty,tz;
  float radio_maximo_priv=nc*10.;			//	radio_maximo_para_refina_aux
  float aux_priv;
  #pragma omp for
  for(aux_void=0;aux_void<=id_void;aux_void++){
    tx = Xv[aux_void]-x;
         if (tx<-nc*0.5)	tx+=1.*nc;
    else if (tx>nc*0.5)		tx-=1.*nc;

    ty = Yv[aux_void]-y;
         if (ty<-nc*0.5)	ty+=1.*nc;
    else if (ty>nc*0.5)		ty-=1.*nc;

    tz = Zv[aux_void]-z;
         if (tz<-nc*0.5)	tz+=1.*nc;
    else if (tz>nc*0.5)		tz-=1.*nc;

    aux_priv=sqrt(tx*tx+ty*ty+tz*tz)-Radv[aux_void];	// esta es la distancia al vorde del void aux_void
    if(radio_maximo_priv>aux_priv)
       radio_maximo_priv=aux_priv;
  }
  #pragma omp critical
  {
  if(radio_maximo>radio_maximo_priv)
     radio_maximo=radio_maximo_priv;
  }
  }// pragma

  radio_maximo_para_refinar=radio_maximo;

}













// busca la esfera subdensa más grande alrededor de x,y,z.
// tiene como límite de búsqueda radio_maximo_para_refinar
// Devuelve su rádio y densidad
// esta funcion solo se llama cuando radio_maximo_para_refinar > rad_min_actual
void crece_refina(float x,float y,float z,float &radio,float &rho,float &rho_fail){
  int i,j,k,ii,jj,kk,iii,jjj,kkk,ijk,l;
  float dx,dy,dz,r,Rmax;
  float dista,c;//,radio_minimo;
  int permanencia=1;
  float const_vol=Mp*3./(4.*pi*pow(Lbox/nc,3));

//	14 de marzo de 2018
//	parece injusto que se busque en un rango fijo de 5 Mpc pata voids de radio 20 y para voids de radio 2
//	para voids grandes parece bien buscar entre 20% más y 20% menos
//	pero pava voids pequeños, r=2, dicha franja está entre 1.6 y 2.4, lo cual parece muy pequeño
//	y dado que la cantidad de voids pequeños es grande, no es bueno que estos tamaños no estén optimizados
//	voy a probar con 20% para ver que pasa, si el tiemp de calculo no se incrementa mucho lo dejo así

  Rmax=radio*1.2;
  if(Rmax>radio_maximo_para_refinar)	//	11 de Abrild e 2018
     Rmax=radio_maximo_para_refinar;
  r=Rmax*0.694444;
  if(r<rad_min_actual)
     r=rad_min_actual;
  if(r>=Rmax){
    printf("\nerror\nlos limites para crece refina xon (%f,%f)\n",r,Rmax);fflush(stdout);}

while(permanencia>0){

  float Rmax2=Rmax*Rmax;
  float r2   =r*r;
  int caso=-1;

  // volumen y radio de cada bin
  c=(Rmax-r)/(bin-1);
  #pragma omp parallel for
  for(i=0;i<bin;i++){
    vol[i]=const_vol*pow(r+i*c,-3);	//	en [Ms h² / Mpc³]
    cont[i]=0.;						//	# de partículas en cada bin
  }

  float aux=Rmax+sqrt(3.);

  ii=floor(x-aux+0.5);	//	el (3/2)^0.5 se debe a que en la diagonal de 45º
  jj=floor(y-aux+0.5);	//	las células a esta distancia corregida
  kk=floor(z-aux+0.5);	//	tienen partículas a la distancia rad+2
  iii=floor(x+aux+0.5);
  jjj=floor(y+aux+0.5);
  kkk=floor(z+aux+0.5);

  float aux_aux=Rmax2+3./4+Rmax*sqrt(3.);

  for(k=kk;k<=kkk;k++){
    dz=z-1.*k;
    for(j=jj;j<=jjj;j++){
      dy=y-1.*j;
      for(i=ii;i<=iii;i++){
        ijk=((i+nc)%nc)+((j+nc)%nc)*nc+((k+nc)%nc)*nc2;

	if(Ocu[ijk]>0){
          dx=x-1.*i;
           
          dista=dx*dx+dy*dy+dz*dz;

          if(dista<aux_aux){  //	celula que tiene partículas dentro del rádio requerido

//	    #pragma omp parallel
//	    {
	    int id_priv;
	    int int_priv;
	    float dist_priv;
	    float aux2_priv;
            float dx_priv;
            float dy_priv;
            float dz_priv;
//	    int cont_priv[bin]={0};
//	    #pragma omp for
            for(l=0;l<Ocu[ijk];l++){
              id_priv=Id[ijk][l];
              dx_priv=x-X[id_priv];
                   if(dx_priv<-nc*0.5) dx_priv+=1.*nc;
              else if(dx_priv>nc*0.5)  dx_priv-=1.*nc;

              dy_priv=y-Y[id_priv];
                   if(dy_priv<-nc*0.5) dy_priv+=1.*nc;
              else if(dy_priv>nc*0.5)  dy_priv-=1.*nc;

              dz_priv=z-Z[id_priv];
                   if(dz_priv<-nc*0.5) dz_priv+=1.*nc;
              else if(dz_priv>nc*0.5)  dz_priv-=1.*nc;

              dist_priv=dx_priv*dx_priv + dy_priv*dy_priv + dz_priv*dz_priv;

              if (dist_priv<Rmax2){
                if(dist_priv<=r2)
                  aux2_priv=-0.1;
                else
                  aux2_priv=(sqrt(dist_priv)-r)/c;
                int_priv=ceil(aux2_priv);	//	es el bin radial de esta particula
//                cont_priv[int_priv]+=1;	//	contador entero del número de partículas
                cont[int_priv]+=1.;	//	contador entero del número de partículas
	      }
	    }
//	    #pragma omp critical
//	    {
//	    for(l=0;l<bin;l++)
//	    cont[l]+=1.*cont_priv[l];
//	    }
//	    }
	  }	//	if de la distancia al centro
	}	//	if del número de ocupación
      }
    }
  }

  for(l=1;l<bin;l++){
    cont[l]+=cont[l-1];}

  l=bin-1;
  rho=cont[l]*vol[l];	//	densidad media hasta el ultimo nivel en unidades [Ms h² / Mpc³]
  rho_fail=cont[0]*vol[0];

  if(rho>=limite){
    while((rho>limite)&&(l>0)){
      l--;
      rho=cont[l]*vol[l];}
    if(rho<=limite){
      rho_fail=cont[l]*vol[l];	// para verificar, debe ser menor que limite
      l++;
      rho=cont[l]*vol[l];
      radio=r+c*l;
      permanencia=-1;
    }
    else{
      if(r>rad_min_actual){	// si puedo disminuir -> disminuyo, no salgo
        r*=0.833334;
        if(r<rad_min_actual)
           r=rad_min_actual;
        Rmax=r*1.44;
        permanencia=1;}
      else{// si no puedo disminuir -> salgo
        radio=r;
        rho=limite*0.5;
        rho_fail=1.;
        permanencia=-1;}
    }
  }
  else{
    if(Rmax<radio_maximo_para_refinar){ // si puedo aumentar -> aumento, no salgo
      Rmax*=1.2;
      if(Rmax>radio_maximo_para_refinar)
         Rmax=radio_maximo_para_refinar;
      r=Rmax*0.694444;
      permanencia=1;}
    else{// si no puedo aumentar -> salgo
      radio=Rmax;
      rho=limite*0.5;
      rho_fail=0.;
      permanencia=-1;}
  }

}	//	cierra el while


// los muy vacios van a salir como rho=limite*0.99 y rho_fail=0.
// los muy densos van a salir como rho=limite*0.99 y rho_fail=1.

//printf("final crece_refina %f %f %f %f\n",x,y,z,radio);fflush(stdout);
}


















void refino(int l,float &x,float &y,float &z, float &radio, float &rho, float &rho_fail){
  int i,j,k,p,q;

  float escala=radio*despl_max;		//	mayor desplazamiento del centro durante la busqueda
  float presicion=radio*despl_min;	//	error al determinar la posicion del centro

  if (escala>0.5) escala=0.5;
  if (presicion>0.125) presicion=0.125;


  i=l%nc;
  j=(l-i)%nc2;
  k=(l-i-j)/nc2;
  j/=nc;
  x_x[0]=i*1.;
  y_y[0]=j*1.;
  z_z[0]=k*1.;



  float sx,sy,sz;
  float aux_rad;

  rad_rad[0]=radio;
  rho_rho[0]=rho;
  rho_rho_fail[0]=rho_fail;

criterio(x_x[0],y_y[0],z_z[0]);	//	calcula el radio maximo para refinar
if(radio_maximo_para_refinar<rad_min_actual){//  me dice si está o nó dentro de un hueco
  rho_rho[0]=limite*0.5;rad_rad[0]=radio_maximo_para_refinar;}
else
  crece_refina(x_x[0],y_y[0],z_z[0],rad_rad[0],rho_rho[0],rho_rho_fail[0]);//	modifica rho_rho[0] y rad_rad[0]

  while(escala>=presicion){
    q=1;
    for(sz=z_z[0]-escala;sz<=z_z[0]+escala*1.5;sz+=2.*escala){
      for(sy=y_y[0]-escala;sy<=y_y[0]+escala*1.5;sy+=2.*escala){
        for(sx=x_x[0]-escala;sx<=x_x[0]+escala*1.5;sx+=2.*escala){
          x_x[q]=(float) fmod(1.*sx+1.*nc,1.*nc);
          y_y[q]=(float) fmod(1.*sy+1.*nc,1.*nc);
          z_z[q]=(float) fmod(1.*sz+1.*nc,1.*nc);
          rad_rad[q]=radio;
          criterio(x_x[q],y_y[q],z_z[q]);
          if(radio_maximo_para_refinar<rad_min_actual){
            rho_rho[q]=limite*0.5;	rad_rad[q]=radio_maximo_para_refinar;}
          else
            crece_refina(x_x[q],y_y[q],z_z[q],rad_rad[q],rho_rho[q],rho_rho_fail[q]);
          q++;
	}
      }
    }
// los muy vacios van a tener rho=limite*0.5 y rho_fail=0.
// los muy densos van a tener rho=limite*0.5 y rho_fail=1.
    aux_rad=0.;
    p=0;
    for(q=8;q>=0;q--){
      if((aux_rad<=rad_rad[q])&&(rho_rho[q]>limite)){
        aux_rad=rad_rad[q];
        p=q;
      }
    }
/*    aux_rho=1.;
    for(q=8;q>=0;q--){
      if(aux_rad==rad_rad[q]){
        if(aux_rho>=rho_rho[q]){
          aux_rho=rho_rho[q];
          p=q;
	}
      }
    }
*/
    if(p==0){   //	si el centro es mejor candidato que sus vecinos reduce el radio de busqueda
      escala*=0.5;
    }
    else{	//	si un vecino es mejor candidato, lo convierto en el centro de busqueda
      x_x[0]=x_x[p];
      y_y[0]=y_y[p];
      z_z[0]=z_z[p];
      rad_rad[0]=rad_rad[p];
      rho_rho[0]=rho_rho[p];
      rho_rho_fail[0]=rho_rho_fail[p];
//      escala*=0.5;              //	cuando quiera desplazarme debo comentar esta línea
    }


///
///
///
//	si radio<radio_minimo entonces salgo, y anulo todo lo que pueda

//printf("p=%d ",p);fflush(stdout);
  }             //	cierra el while de escala



  radio=rad_rad[0];
  rho=rho_rho[0];
  rho_fail=rho_rho_fail[0];
  x=x_x[0];
  y=y_y[0];
  z=z_z[0];
// los muy vacios van a tener rho=limite*0.5 y rho_fail=0.
// los muy densos van a tener rho=limite*0.5 y rho_fail=1.

//printf("\nrefino: id=%d l=%d i=%d j=%d k=%d x=%lf y=%lf z=%lf radio=%lf rho=%lf rho_fail=%lf\n"
//       ,id_void,l,i,j,k,x,y,z,radio,rho/Om/rho_cr,rho_fail/Om/rho_cr);fflush(stdout);

}

















inline int fi(int indice,int radio){
  return -radio+indice%(2*radio+1);
  }
inline int fj(int indice,int radio){
  int aux=floor(1.*indice/(2*radio+1));
  return -radio+aux%(2*radio+1);
  }
inline int fk(int indice,int radio){
  return -radio+floor(1.*indice*pow(2.*radio+1.,-2.));
  }













//	lee xyz, y llena los IDs
void lee(){
  printf("lee xyz...\n");
  fflush(stdout);

  int i,j,k,ijk,p;
  float x,y,z,u,v,w;

  FILE * LEE;
  LEE = fopen(catalogo_particulas,"r");

  // llena el vector de posiciones y cuenta las particulas que hay en cada celda
  for (p=0;p<NumPart;p++)
    {
    fscanf(LEE,"%f %f %f %f %f %f\n",&x,&y,&z,&u,&v,&w);
    X[p]=fmod(x*nc/Lbox,1.*nc);
    Y[p]=fmod(y*nc/Lbox,1.*nc);
    Z[p]=fmod(z*nc/Lbox,1.*nc);
    i=floor(X[p]+0.5);	//	debo sumar 0.5 porque las células están centradas
    j=floor(Y[p]+0.5);	//	en los vértices enteros del sistema de coordenadas
    k=floor(Z[p]+0.5);	//	Nó en puntos tales como 1.5, 45.5, 178.5
    ijk=(i%nc)+(j%nc)*nc+(k%nc)*nc*nc;
    Ocu[ijk]++;
    }
  fclose (LEE);
  printf("lee xyz...hecho\n");
  fflush(stdout);

  allocaId();

  for(p=0;p<NumCel;p++)
    Ocu[p]=0;

  printf("Llena ID...\n");
  fflush(stdout);
  for(p=0;p<NumPart;p++)//loop sobre las artículas
    {
    i=floor(X[p]+0.5);	//	debo sumar 0.5 porque las células están centradas
    j=floor(Y[p]+0.5);	//	en los vértices enteros del sistema de coordenadas
    k=floor(Z[p]+0.5);	//	Nó en puntos tales como 1.5, 45.5, 178.5
    ijk=(i%nc)+(j%nc)*nc+(k%nc)*nc2;
    Id[ijk][Ocu[ijk]]=p;
    Ocu[ijk]++;
    }

  printf("Llena Id...hecho\n");
  fflush(stdout);

}












//	calcula da densidad celular a partir de las partículas
void densidad(){
  printf("\ncalcula la densidad de las celdas...\n");
  fflush(stdout);

  int i,j,k,ii,jj,kk,p;
  float sx,sy,sz;
  float tx,ty,tz;

  for(p=0;p<NumPart;p++)//loop sobre las particulas para calcular la densidad
    {
    i=floor(X[p]);        j=floor(Y[p]);          k=floor(Z[p]);
    sx=X[p]-1.*i;         sy=Y[p]-1.*j;           sz=Z[p]-1.*k;
    tx=1.-sx;             ty=1.-sy;               tz=1.-sz;

    // condiciones de frontera periodicas
    ii=(i+1)%nc;          jj=(j+1)%nc;            kk=(k+1)%nc;

    // incrementando la densidad de las celdas adyacentes
    Rho[i+nc*j +nc2*k] +=Mp*tx*ty*tz;         Rho[ii+nc*j +nc2*k] +=Mp*sx*ty*tz;
    Rho[i+nc*jj+nc2*k] +=Mp*tx*sy*tz;         Rho[ii+nc*jj+nc2*k] +=Mp*sx*sy*tz;
    Rho[i+nc*j +nc2*kk]+=Mp*tx*ty*sz;         Rho[ii+nc*j +nc2*kk]+=Mp*sx*ty*sz;
    Rho[i+nc*jj+nc2*kk]+=Mp*tx*sy*sz;         Rho[ii+nc*jj+nc2*kk]+=Mp*sx*sy*sz;
    // Rho sería la masa física contenida en cada célula
    }

  printf("calcula densidad de las celdas...hecho\n");
  fflush(stdout);
}












// busca las esferas subdensas a partir de las densidades celulares
// busco la esfera subdensa mas grande centrada cada célula (en cada vértice entero del plano)
// le asigno el rádio (entero*1.)
// cambio: octubre 30: solo crece alrededor de las esferas subdensas
void crece(){
  printf("\nradio máximo de búsqueda usando células = %f [Mpc/h]\n",otro*Lbox/nc);
  printf("Encuentra el radio de la mayor esfera subdensa centrada en cada célula...\n");
  fflush(stdout);

  int i,j,k,l,p,ijk;
  float rho;

  FILE * RAD;
  RAD=fopen(densidad_radios,"w+");

  // incremento la densidad de todas las esferas
  int max=incremento;		//	radio máximo en unidades de nc
  int max2=incremento*incremento;

  for(k=0;k<nc;k++){				//	ciclo sobre las células
    for(j=0;j<nc;j++){
      for(i=0;i<nc;i++){

      ijk=i+nc*j+nc2*k;
      //	nuevo: optimiza: solo crezco esferas alrededor de pixeles subdensos
      if(Rho[ijk]>limite)
        Rad[ijk]=0.;

      else{

	int sale=1;
	while(sale>0){	//	buca para un radio, si no funciona lo aumenta,etc

	  // borro la información antigua
	  #pragma omp parallel for
 	  for(p=0;p<semi*max;p++){
	    cont[p]=0.;
	    vol[p]=0.;
	    }

	  #pragma omp parallel
	  {
	  int i_priv;
	  int j_priv;
	  int k_priv;
	  int ijk_priv;
	  int bin_priv;
	  int dist_priv;
	  int vol_priv[semi_otro_max]={0};
	  float mas_priv[semi_otro_max]={0.};
	  #pragma omp for
	  for(p=0;p<(2*max+1)*(2*max+1)*(2*max+1);p++){
	    k_priv=k+fk(p,max);
	    j_priv=j+fj(p,max);
	    i_priv=i+fi(p,max);
            dist_priv=(k_priv-k)*(k_priv-k)+(j_priv-j)*(j_priv-j)+(i_priv-i)*(i_priv-i);
            if(dist_priv<max2){		//	si es vecina de verdad...
	      bin_priv=floor(semi*sqrt(dist_priv));
              ijk_priv=(i_priv+nc)%nc+nc*((j_priv+nc)%nc)+nc2*((k_priv+nc)%nc);
	      mas_priv[bin_priv]+=Rho[ijk_priv];	//	masa física contenida en el bin l
	      vol_priv[bin_priv]++;		//	volumen en unidades de nc
	    }
	  }
	  #pragma omp critical
	  {
	  for(p=0;p<semi*max;p++){
	    cont[p]+=mas_priv[p];
	    vol[p]+=vol_priv[p]*1.;
	  }
	  }
	  }	//	sale del pragma

	  for(l=1;l<semi*max;l++){
	    cont[l]+=cont[l-1];	//	masa física de la esfera l
	    vol[l]+=vol[l-1];	//	volumen de la esfera l en unidades de nc
	  }

	  // hacer las cosas de fuera para dentro nos permite detectar voids con halos dentro
	  l=semi*max-1;
	  rho=cont[l]/vol[l]*pow(nc/Lbox,3);	//	densidad en [Ms h² / Mpc³]

	  if(rho<limite){	//	debo aumentar el radio de busqueda
	    max+=incremento;
	    max2=max*max;
	  }

	  else{		//	si todo esta bien
	    max=incremento;
	    max2=max*max;
	    sale=-1;
	    while((rho>=limite)&&(l>-1)){
	      rho=cont[l]/vol[l]*pow(nc/Lbox,3);	//	densidad en [Ms h² / Mpc³]
	      l--;
	    }
	  }

	  l++;	//	porque l=0 es la esfera de radio 1

	  Rad[ijk]=l*1./semi;		//	radio en unidades de nc
	}	//	cierra el while de sale
      }		//	else: si es subdensa
      fprintf(RAD,"%f\n",Rad[ijk]);
      }		//	i
    }		//	j
  }		//	k

  fclose(RAD);

  printf("ya tenemos los radios de las mayores esferas subdensas centradas en cada célula\n\n");
  fflush(stdout);
}

















void lee_radios(){
  printf("\nlee el archivo con radios y desidades de celulas...\n");
  fflush(stdout);

  float dist;

  FILE * RAD;
  RAD=fopen(densidad_radios,"r");
  for(int p=0;p<NumCel;p++){
    fscanf(RAD,"%f\n",&dist);
    Rad[p]=dist;
  }
  fclose(RAD);
  printf("lee el archivo con radios y desidades de celulas...hecho\n");
  fflush(stdout);
}












//	busca el contro de celula con el radio más grande de subdensidad 0.2
int busco(void){
  float radio=rad_min_actual;
//  float rho=halo;
  int l=-1;
  int i;
  #pragma omp parallel
  {
  float radio_priv=radio;	//	=radio_min_actual
//  float rho_priv=rho;
  int   l_priv=l;		//	=-1
  #pragma omp for
  for(i=0;i<NumCel;i++){	//	solo busco el radio más grande
    if(radio_priv<Rad[i]){
      radio_priv=Rad[i];	//	radio en unidades de nc
      l_priv=i;			//	es el índice de la célula con la esfera más grande
    }
//    else if(radio_priv==Rad[i]){
//      if(rho_priv>Den[i]){
//        rho_priv=Den[i];		//	densidad física
//        l_priv=i;			//	es el índice de la célula con la esfera más grande
//      }
//    }
  }
  #pragma omp critical
  {
  if(radio<radio_priv){
    radio=radio_priv;	//	radio en unidades de nc
//    rho=rho_priv;		//	densidad física
    l=l_priv;			//	es el índice de la célula con la esfera más grande
  }
//  else if(radio==radio_priv){
//    if(rho>rho_priv){
//      rho=rho_priv;		//	densidad física
//      l=l_priv;			//	es el índice de la célula con la esfera más grande
//    }
//  }
  }				//	ya tengo el radio máximo
  }
  return l;
}












//	libera X Y Z ID
void free_xyz(void){
  free(X);
  free(Y);
  free(Z);
  free(Xv);
  free(Yv);
  free(Zv);
  free(Radv);
  free(Ocu);
  for(int p=0;p<NumCel;p++)
    free(Id[p]);
  free(Id);
}










//	realloca Xv Yv Zv, y llena dos ID
void realloca(void){
  printf("\nrealloca V...\n");
  fflush(stdout);
  if(!(Xv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Xv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Yv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Yv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Zv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Zv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Radv = (float*) calloc (NumHuecos , sizeof(float))))
    {
      fprintf(stderr, "failed to allocate memory Radv.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(IDV = (int*) calloc (NumHuecos , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory IDV.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(ID_random = (int*) calloc (NumHuecos , sizeof(int))))
    {
      fprintf(stderr, "failed to allocate memory ID_random.\n");
      fflush(stderr);
      exit(0);
    }
  printf("realloca V...hecho\n");
  fflush(stdout);
}


















int main(int argc, char **argv){
  printf("comienza\n");
  fflush(stdout);
  omp_set_num_threads(24);


  FILE * ESC;	//	catalogo de trabajo
  int l;
  float x,y,z,radio,rho,rho_fail;

  archivos();		//	carga el nombre de los archivos

  allocar();
  lee();
  densidad();


  rad_min_actual=rad_min_uno;

  if(calcula==1)
    crece();
  else
    lee_radios();

  // Busco la esfera más grande
  l=busco();	//	es el índice de la célula con la mayor esfera subdensa
   printf("ya buscó la célula con el rádio más grande\n");
  fflush(stdout);

  // de aquí en adelante x, y, z son las coordenadas del centro del cada void
  radio=Rad[l];printf("su radio es %f\n",Rad[l]);fflush(stdout);
  rho=Rho[l];
  rho_fail=Rho[l];
  Rad[l]=-1.;

  refino(l,x,y,z,radio,rho,rho_fail);		//	refina la esfera de la célula l
// los muy vacios van a tener rho=limite*0.99 y rho_fail=0.
// los muy densos van a tener rho=limite*0.99 y rho_fail=1.

  printf("radio refinado del void más grande = %f\n",radio);
  printf("densidad refinada del void más grande = %f\n",rho);
  printf("centro refinado del void más grande = (%f, %f, %f)\n\n",x,y,z);  fflush(stdout);

  // escribo en los dos catalogos (trabajo y respaldo)
  ESC = fopen(catalogo_voids,"w+");
  fprintf(ESC,"%f %f %f %f %f %f\n",x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc,rho/Om/rho_cr,rho_fail/Om/rho_cr);
  fclose(ESC);

  id_void++;
  Xv[id_void]=x;	//	en unidades de nc
  Yv[id_void]=y;
  Zv[id_void]=z;
  Radv[id_void]=radio;

  anula(x,y,z,radio);		//	anula el rádio de las celdas \in void
				//	y aumenta la densidad de las células \in void
  //  crece_celulas(x,y,z,radio);	//	recalcula los radios celulares vecinos


  printf("ya encontramos, refinamos, anulamos e imprimímos el Void 0\n\n");
  fflush(stdout);








  int condicion=1;		//	permanencia dentro del ciclo principal
  printf("entra al ciclo principal\n");
  fflush(stdout);

  while(condicion>=0){		//	ciclo principal

    // Busco la siguiente esfera mas grande
    condicion=busco();
    //	ya identifiqué la siguiente esfera más grande, condicion es el índice de su celda central
    //  si no hay más esferas con estas características, la condición debe ser -1

    if(condicion==-1){
      printf("en este punto deberia salir del ciclo principal\n\n");
      fflush(stdout);
    }

    else{
      radio=Rad[condicion];
      rho=Rho[condicion];
      rho_fail=Rho[condicion];
      refino(condicion,x,y,z,radio,rho,rho_fail);	//	refina la esfera de la célula l=condicion
// los muy vacios van a tener rho=limite*0.99 y rho_fail=0.
// los muy densos van a tener rho=limite*0.99 y rho_fail=1.
      if(rho>limite){
        ESC = fopen(catalogo_voids,"a");
        fprintf(ESC,"%f %f %f %f %f %f\n",
                     x*Lbox/nc,y*Lbox/nc,z*Lbox/nc,radio*Lbox/nc,rho/Om/rho_cr,rho_fail/Om/rho_cr);
        fclose(ESC);}
      id_void++;
      Xv[id_void]=x;	//	en unidades de nc
      Yv[id_void]=y;
      Zv[id_void]=z;
      Radv[id_void]=radio;
      Rad[condicion]=-1.;
      anula(x,y,z,radio);		//	anula el rádio de las celdas \in void
    }	//	cierra el else
  }	//	cierra el ciclo principal




  printf("de pura casualidad llega a este punto? 2?\n");fflush(stdout);
  free_xyz();	//	libera X Y Z ID
  free_pixelizado();	//	libera radio y densidad de celulas

  printf("este es el fin del camino\n");fflush(stdout);



  return 0;

}



