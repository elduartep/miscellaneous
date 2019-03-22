//Calcula el espectro de potencias para z=0
//tiene en cuenta distorciones de redshift pixelizadas
//Calcula el bias lineal de voids y halos

//paquetes gnu
#include<iostream>
#include<fstream>
#include<cmath>
#include <alloca.h>
//OPENMP
#include <omp.h>
//FFTW
#include "fftw3.h"
//numerical recipes
#include "../l/nr3.h"
#include "../l/gamma.h"
#include "../l/incgammabeta.h"
#include "../l/fitab.h"
#include "../l/ran.h"
#include "../l/deviates.h"
#include "../l/interp_1d.h"

#include "../l/RANDOM.h"

#include "../l/parametros.h"


using namespace std;

const int puntual=1;

const int distorcion=0;
const int espectro=1;
const int bias_halos=0;
const int bias_voids=1;

const int solo_analiza_bias_voids=0;
const int bin_lineal=1;

const int espectro_voids=0;
const int espectro_halos=0;

const double k_lineal=0.25;

const double rho_cr=2.77526627;		//	[10¹¹ h² Ms / Mpc³]
const int    nc=256;
const int    NumPart=np*np*np;		//	numero total de particulas
const int    NumCel=nc*nc*nc;		//	numero total de particulas
const int    nc2=nc*nc;
const double Mp=pow(Lbox/np,3.)*Om*rho_cr;	//	masa de 1 particula [10¹¹ Ms/h]

const double pi=4.*atan(1.);		//	numero pi

//	k_Nyquist   = (2 pi / Lbox )(nc / 2) = pi nc / Lbox
//	k_fisico = k_entero 2 pi / Lbox == k_entero dpin
//	x_fisico = x_entero Lbox / nc
//	f_fisico = f_fftw Lbox / nc
const double dpin=2*pi/Lbox;
const double kmin=0.0*dpin*0.99;	//	= 0.1*k_fisico(k_entero=1)
const double kmax=0.5*nc*dpin*1.01;	//	= k_Nyquist


  int nc21=(nc/2+1);
  int nc21nc=nc21*nc;

//	P(k)= <dd*> / (2 pi)³
//	P(k)~ <dd*> / (nc)³


//	delta = (delta rho)/(bar rho) = (rho -bar rho)/(bar rho) = (rho/(bar rho)) - 1
//	bar rho = MT/VT = np³ M1 / Lbox³
//	rho = sum fx fy fz M1
//	rho / (bar rho) = sum fx fy fz (Lbox/np)³ == sum fx fy fz f


// en el articulo del bi-espectro usaron una caja de (2048 Mpc/h)^3
// las simulaciones de D. Mota son de (256 Mpc/h)^3
// osea, 8^3 veces más pequeña
// entonces vamos a usar intervalos para k, 8 veces más grandes

// yo usé cajas de (400 Mpc/h), entonces voy a usar intervalos 5 veces más grandes

const int    bin=64;              	//	bines logaritmicos para el espectro
const int    binb=10;              	//	bines logaritmicos para el bias
//const int    bin_k=16;              	//	numero de bines k (0.09, 0.41) 0.02
//const int    bin_kp=10;              	//	numero de bines k'(0.005,0.105) 0.01
//const double k_min=0.09*5;
//const double k_max=0.41*5;
//const double kp_min=0.09*5;
//const double kp_max=0.41*5;
const double limite=0.2*Om*rho_cr;	//	sobredensidad de los voids en [10¹¹ Ms h² / Mpc³]


const int nt=24;                //	numero de nucleos


////////////////////////////            principales variables
double *rho;			//	contraste de densidad
double *x;			//	auxiliar para rsd
double *vx;			//	auxiliar para rsd
double *MP, *Mk, *EP, *Ek;	//      variables para calcular la media del espectro
double *MP2, *EP2;		//      variables para calcular la media del espectro mixto
double *MPb, *EPb;

//      variables para calcular la media del espectro mixto
double *X, *Y, *Z, *R;		//	posición y radio de los halos
double *MPb2, *EPb2;		//      variables para calcular la media del espectro mixto
double *MPb3, *EPb3;		//      variables para calcular la media del espectro mixto
double *MPb4, *EPb4;		//      variables para calcular la media del espectro mixto
double *MPb5, *EPb5;		//      variables para calcular la media del espectro mixto
double *MPb6, *EPb6;		//      variables para calcular la media del espectro mixto
long int *cont;			//      variable para calcular la media del espectro
long int *cont2;			//      variable para calcular la media del espectro
double *carga;			//	masa/radio medio de cada binb
double *Ecarga;			//	dispersión masa/radio medio de cada binb
int    *cantidad;		//	cantidad de halos/voids en cada bin de masa/radio
int    *fuera;			//	fuera o dentro de un void
double cmax,cmin,c;		//	para los bines binb de masa/radio
int NumVoids;			//	número de voids
int NumHalos;			//	número de halos

const char ejes[3]={'X','Y','Z'};	//	nombre de los ejes para aplicar redshift-distortions

double corte,pendiente,segunda,scorte;





fftw_complex *delta;		//      potencial gravitacional
fftw_complex *deltab;		//      potencial gravitacional
fftw_complex *deltac;		//      potencial gravitacional
fftw_plan dft,dftb,dftc;	//      plan principal











inline double kernel(double k){
  double aux=k*Lbox*0.5/nc;
  return pow(Lbox,3.)*pow(aux/sin(aux),4.);
}





inline int fi(int indice,int radio){
  return -radio+indice%(2*radio+1);}

inline int fj(int indice,int radio){
  int aux=floor(1.*indice/(2*radio+1));
  return -radio+aux%(2*radio+1);}

inline int fk(int indice,int radio){
  return -radio+floor(1.*indice*pow(2.*radio+1.,-2.));}





void allocar(){
  if(!(rho = (double*) malloc (NumCel * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory rho.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(delta = (fftw_complex*) fftw_malloc (nc2 * (nc/2+1) * sizeof(fftw_complex))))
    {
      fprintf(stderr, "failed to allocate memory delta.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(deltac = (fftw_complex*) fftw_malloc (nc2 * (nc/2+1) * sizeof(fftw_complex))))
    {
      fprintf(stderr, "failed to allocate memory deltac.\n");
      fflush(stderr);
      exit(0);
    }

  if(!(dft = fftw_plan_dft_r2c_3d(nc, nc, nc, rho, delta, FFTW_MEASURE)))
    {
      fprintf(stderr, "failed to allocate memory dft.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(dftc = fftw_plan_dft_r2c_3d(nc, nc, nc, rho, deltac, FFTW_MEASURE)))
    {
      fprintf(stderr, "failed to allocate memory dftb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(cont  = (long int*) malloc(sizeof(long int) * bin)))
    {
      fprintf(stderr, "failed to allocate memory cont.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(cont2  = (long int*) malloc(sizeof(long int) * bin)))
    {
      fprintf(stderr, "failed to allocate memory cont2.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MP = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MP.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EP = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EP.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Mk  = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory Mk.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Ek  = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory Ek.\n");
      fflush(stderr);
      exit(0);
    }
}











void allocarb(){
  if(!(deltab = (fftw_complex*) fftw_malloc (nc2 * (nc/2+1) * sizeof(fftw_complex))))
    {
      fprintf(stderr, "failed to allocate memory deltab.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(dftb = fftw_plan_dft_r2c_3d(nc, nc, nc, rho, deltab, FFTW_MEASURE)))
    {
      fprintf(stderr, "failed to allocate memory dftb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MP2 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(cantidad = (int*) malloc(sizeof(int) * bin )))
    {
      fprintf(stderr, "failed to allocate memory cantidad.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EP2 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MPb = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EPb = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MPb2 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EPb2 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MPb3 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EPb3 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MPb4 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EPb4 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MPb5 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EPb5 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(MPb6 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory MPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(EPb6 = (double*) malloc(sizeof(double) * bin )))
    {
      fprintf(stderr, "failed to allocate memory EPb.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(carga = (double*) malloc(sizeof(double) * binb )))
    {
      fprintf(stderr, "failed to allocate memory carga.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Ecarga = (double*) malloc(sizeof(double) * binb )))
    {
      fprintf(stderr, "failed to allocate memory Ecarga.\n");
      fflush(stderr);
      exit(0);
    }
}






void alloca_voids(){
  if(!(fuera = (int*) malloc(sizeof(int) * NumCel )))
    {
      fprintf(stderr, "failed to allocate memory fuera.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(X = (double*) malloc(sizeof(double) * NumVoids )))
    {
      fprintf(stderr, "failed to allocate memory X.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Y = (double*) malloc(sizeof(double) * NumVoids )))
    {
      fprintf(stderr, "failed to allocate memory Y.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Z = (double*) malloc(sizeof(double) * NumVoids )))
    {
      fprintf(stderr, "failed to allocate memory Z.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(R = (double*) malloc(sizeof(double) * NumVoids )))
    {
      fprintf(stderr, "failed to allocate memory R.\n");
      fflush(stderr);
      exit(0);
    }

}







void allocard(){
  if(!(x = (double*) malloc (NumCel * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory x.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(vx = (double*) malloc (NumCel * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory vx.\n");
      fflush(stderr);
      exit(0);
    }
}























void max_min_c(const char *NomArch){
  printf("busca maximo y minimo...\n");
  fflush(stdout);
  float a1, a2, a3, a4, a5;
  int i;

  cmin=100000000.,cmax=0.;

  FILE * IN;
  IN=fopen(NomArch,"r");
  NumHalos=0;
  while(fscanf(IN,"%f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5)!=EOF){
    if(cmin>double(a4)) cmin=double(a4);
    if(cmax<double(a4)) cmax=double(a4);
  NumHalos++;
  }
  fclose(IN);

  printf("número de halos = %d\n",NumHalos);
  fflush(stdout);

  printf("cmax=%lf cmin=%lf para este caso\n",cmax,cmin);fflush(stdout);

  cmax=masa_max_todos;
  cmin=masa_min_todos;

  printf("cmax=%lf cmin=%lf globales escogidos a mano, y usados de aqui en adelante\n",cmax,cmin);fflush(stdout);


//  cmin*=0.99;
//  cmax*=1.01;
  c   = pow(cmax/cmin,1./binb);       //	constante de proporcionalidad log
  double aux1 = 1./log(c);
  double aux2;
  int bin_priv;

  for(i=0;i<binb;i++){
    cantidad[i]=0;
    carga[i]=0.;
    Ecarga[i]=0.;
  }

  IN=fopen(NomArch,"r");	//	carga media de cada bin
  for(i=0;i<NumHalos;i++){
    fscanf(IN,"%f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5);
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2=aux1*log(double(a4)/cmin);
      bin_priv=floor(aux2);	//	es el bin de esta particula
      carga[bin_priv]+=double(a4);
      cantidad[bin_priv]++;
    }
  }
  fclose(IN);

  for(i=0;i<binb;i++){
    carga[i]/=cantidad[i];
  }

  IN=fopen(NomArch,"r");	//	dispersion de la carga media en cada bin
  for(i=0;i<NumHalos;i++){
    fscanf(IN,"%f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5);
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2=aux1*log(double(a4)/cmin);
      bin_priv=floor(aux2);	//	es el bin de esta particula
      Ecarga[bin_priv]+=pow(double(a4)-carga[bin_priv],2.);
    }
  }
  fclose(IN);

  for(i=0;i<binb;i++){
    aux2=Ecarga[i]/(cantidad[i]*(cantidad[i]-1));
    Ecarga[i]=sqrt(aux2);
  }

  printf("carga_max=%f, carga_min=%f\n",cmax,cmin);
  printf("busca maximo y minimo...echo\n");
  fflush(stdout);
}


















void max_min_cv(const char *NomArch){
  printf("busca maximo y minimo...\n");
  fflush(stdout);
  float a1, a2, a3, a4;
  int i;

  cmin=100000000.,cmax=0.;

  FILE * IN;
  IN=fopen(NomArch,"r");
  NumVoids=0;
  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){
    if(cmin>double(a4)) cmin=double(a4);
    if(cmax<double(a4)) cmax=double(a4);
    NumVoids++;
  }
  fclose(IN);

  printf("número de voids = %d\n",NumVoids);
  fflush(stdout);

  printf("cmax=%lf cmin=%lf para este caso\n",cmax,cmin);fflush(stdout);

  cmax=radio_max_todos;
  cmin=radio_min_todos;

  printf("cmax=%lf cmin=%lf globales escogidos a mano, y usados de aqui en adelante\n",cmax,cmin);fflush(stdout);

//  cmin*=0.99;
//  cmax*=1.01;
  c   = pow(cmax/cmin,1./binb);       //	constante de proporcionalidad log
  double aux1 = 1./log(c);
  double aux2;
  int bin_priv;


  alloca_voids();

  IN=fopen(NomArch,"r");
  for(i=0;i<NumVoids;i++){	//	guarda el catalogo de voids
    fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4);
    X[i]=double(a1*nc)/Lbox;
    Y[i]=double(a2*nc)/Lbox;
    Z[i]=double(a3*nc)/Lbox;
    R[i]=double(a4*nc)/Lbox;
  }
  fclose(IN);

  for(i=0;i<binb;i++){
    cantidad[i]=0;
    carga[i]=0.;
    Ecarga[i]=0.;
  }

  for(i=0;i<NumVoids;i++){	//	radio medio de cada bin
    a4=float(R[i]*Lbox)/nc;
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2=aux1*log(double(a4)/cmin);
      bin_priv=floor(aux2);	//	es el bin de esta particula
      carga[bin_priv]+=double(a4);
      cantidad[bin_priv]++;
    }
  }

  for(i=0;i<binb;i++){
    carga[i]/=cantidad[i];
  }

  for(i=0;i<NumVoids;i++){	//	dispersion en el radio medio de cada bin
    a4=float(R[i]*Lbox)/nc;
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2=aux1*log(double(a4)/cmin);
      bin_priv=floor(aux2);	//	es el bin de esta particula
      Ecarga[bin_priv]+=pow(double(a4)-carga[bin_priv],2.);
    }
  }

  for(i=0;i<binb;i++){
    aux2=Ecarga[i]/(cantidad[i]*(cantidad[i]-1));
    Ecarga[i]=sqrt(aux2);
  }

  printf("carga_max=%f, carga_min=%f\n",cmax,cmin);
  printf("busca maximo y minimo...echo\n");
  fflush(stdout);
}
















void densidad(const char *NomArch, double dr, int cuantos){
  printf("calcula la densidad...\n");
  fflush(stdout);
  float a0, a1, a2, a3, a4, a5, a6;
  float dx, dy, dz, tx, ty, tz;
  int i,j,k,ii,jj,kk,n=0;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]=0.;

  FILE * IN;
  IN=fopen(NomArch,"r");
  //loop sobre las particulas para calcular la densidad de las celdas

  if(cuantos==6)
  while(fscanf(IN,"%f %f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5,&a6)!=EOF){	// coordenadas físicas
    n++;
    a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
    a1=fmod(1.*a4+dr,1.*nc);	a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
    i=floor(1.*a1);		j=floor(1.*a2);			k=floor(1.*a3);
    dx=a1-1.*i;			dy=a2-1.*j;			dz=a3-1.*k;
    tx=1.-dx;			ty=1.-dy;			tz=1.-dz;
    //condiciones de frontera periodicas
    ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
    //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
    rho[i+nc*j +nc2*k ]+=1.*tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=1.*dx*ty*tz;
    rho[i+nc*jj+nc2*k ]+=1.*tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=1.*dx*dy*tz;
    rho[i+nc*j +nc2*kk]+=1.*tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=1.*dx*ty*dz;
    rho[i+nc*jj+nc2*kk]+=1.*tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=1.*dx*dy*dz;}

  else if(cuantos==5)
  while(fscanf(IN,"%f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5)!=EOF){	// coordenadas físicas
    n++;
    a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
    a1=fmod(1.*a4+dr,1.*nc);	a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
    i=floor(1.*a1);		j=floor(1.*a2);			k=floor(1.*a3);
    dx=a1-1.*i;			dy=a2-1.*j;			dz=a3-1.*k;
    tx=1.-dx;			ty=1.-dy;			tz=1.-dz;
    //condiciones de frontera periodicas
    ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
    //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
    rho[i+nc*j +nc2*k ]+=1.*tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=1.*dx*ty*tz;
    rho[i+nc*jj+nc2*k ]+=1.*tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=1.*dx*dy*tz;
    rho[i+nc*j +nc2*kk]+=1.*tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=1.*dx*ty*dz;
    rho[i+nc*jj+nc2*kk]+=1.*tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=1.*dx*dy*dz;}

  else if(cuantos==4)
  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){	// coordenadas físicas
    n++;
    a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
    a1=fmod(1.*a4+dr,1.*nc);	a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
    i=floor(1.*a1);		j=floor(1.*a2);			k=floor(1.*a3);
    dx=a1-1.*i;			dy=a2-1.*j;			dz=a3-1.*k;
    tx=1.-dx;			ty=1.-dy;			tz=1.-dz;
    //condiciones de frontera periodicas
    ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
    //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
    rho[i+nc*j +nc2*k ]+=1.*tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=1.*dx*ty*tz;
    rho[i+nc*jj+nc2*k ]+=1.*tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=1.*dx*dy*tz;
    rho[i+nc*j +nc2*kk]+=1.*tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=1.*dx*ty*dz;
    rho[i+nc*jj+nc2*kk]+=1.*tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=1.*dx*dy*dz;}
  fclose(IN);

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)
    rho[i]/=n;

  printf("calcula la densidad...echo\n");
  fflush(stdout);
}












void densidad_halos(const char *NomArch, double dr){
  printf("calcula la densidad...\n");
  fflush(stdout);
  float a0, a1, a2, a3, a4, a5, a6;
  float dx, dy, dz, tx, ty, tz;
  int i,j,k,ii,jj,kk,n=0;
  double aux_m;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]=0.;

  FILE * IN;
  IN=fopen(NomArch,"r");
  //loop sobre las particulas para calcular la densidad de las celdas

  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){	// coordenadas físicas
    a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
    a1=fmod(1.*a4+dr,1.*nc);	a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
    i=floor(1.*a1);		j=floor(1.*a2);			k=floor(1.*a3);
    dx=a1-1.*i;			dy=a2-1.*j;			dz=a3-1.*k;
    tx=1.-dx;			ty=1.-dy;			tz=1.-dz;
    aux_m=Mp/a4;
    //condiciones de frontera periodicas
    ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
    //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
    rho[i+nc*j +nc2*k ]+=aux_m*tx*ty*tz;	rho[ii+nc*j +nc2*k ]+=aux_m*dx*ty*tz;
    rho[i+nc*jj+nc2*k ]+=aux_m*tx*dy*tz;	rho[ii+nc*jj+nc2*k ]+=aux_m*dx*dy*tz;
    rho[i+nc*j +nc2*kk]+=aux_m*tx*ty*dz;	rho[ii+nc*j +nc2*kk]+=aux_m*dx*ty*dz;
    rho[i+nc*jj+nc2*kk]+=aux_m*tx*dy*dz;	rho[ii+nc*jj+nc2*kk]+=aux_m*dx*dy*dz;}
  fclose(IN);

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]/=NumHalos;

  printf("calcula la densidad...echo\n");
  fflush(stdout);
}



















void densidad_voids(const char *NomArch, double dr){
  printf("calcula la densidad...\n");
  fflush(stdout);
  float a0, a1, a2, a3, a4, a5, a6;
  double dx, dy, dz, tx, ty, tz;
  int i,j,k,ii,jj,kk,iii,jjj,kkk,p,q;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++){//loop sobre las celdas para borrar la densidad antigua
    rho[i]=0.;
    fuera[i]=1;}

  FILE * IN;
  IN=fopen(NomArch,"r");
  //	contribución de la matéria en los voids
  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){	// coordenadas físicas
    a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
    a1=fmod(1.*a4+dr,1.*nc);	a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
    i=floor(1.*a1);		j=floor(1.*a2);			k=floor(1.*a3);
    dx=1.*a1-1.*i;		dy=1.*a2-1.*j;			dz=1.*a3-1.*k;
    tx=1.-dx;			ty=1.-dy;			tz=1.-dz;
    //condiciones de frontera periodicas
    ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
    //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
    rho[i+nc*j +nc2*k ]+=tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=dx*ty*tz;
    rho[i+nc*jj+nc2*k ]+=tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=dx*dy*tz;
    rho[i+nc*j +nc2*kk]+=tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=dx*ty*dz;
    rho[i+nc*jj+nc2*kk]+=tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=dx*dy*dz;}
  fclose(IN);

  //	borrando los bordes de la densidad en los voids
  //	e 'inviertiendo' la densidad en los voids (dens=0 -> delta máximo, dens=maximo -> delta=0)

  for(q=0;q<NumVoids;q++){
    double x,y,z,r;
    double maxima=0.,suma=0.;
    int contador=0;
    x=X[q]+dr;
    y=Y[q]+dr;
    z=Z[q]+dr;
    r=R[q];
    i=floor(x);
    j=floor(y);
    k=floor(z);
    int max=ceil(R[q]);    

    #pragma omp parallel
    {
    int    i_priv;
    int    j_priv;
    int    k_priv;
    int    ijk_priv;
    int    cont_priv=0;
    double dist_priv;
    double max_priv=0.;
    double sum_priv=0.;
    double aux_priv;
    #pragma omp for
    for(p=0;p<(2*max+1)*(2*max+1)*(2*max+1);p++){
      k_priv=k+fk(p,max);
      j_priv=j+fj(p,max);
      i_priv=i+fi(p,max);
      dist_priv=pow(1.*k_priv-z,2)+pow(1.*j_priv-y,2)+pow(1.*i_priv-x,2);
      if(dist_priv<=r*r){
        ijk_priv=(i_priv+nc)%nc+nc*((j_priv+nc)%nc)+nc2*((k_priv+nc)%nc);
	sum_priv+=rho[ijk_priv];
	cont_priv++;
	if(max_priv<rho[ijk_priv])
	  max_priv=rho[ijk_priv];
      }
    }

    #pragma critical
    {
    if(maxima<max_priv)
      maxima=max_priv;
    suma+=sum_priv;
    contador+=cont_priv;
    }

    #pragma omp for
    for(p=0;p<(2*max+1)*(2*max+1)*(2*max+1);p++){
      k_priv=k+fk(p,max);
      j_priv=j+fj(p,max);
      i_priv=i+fi(p,max);
      dist_priv=pow(1.*k_priv-z,2)+pow(1.*j_priv-y,2)+pow(1.*i_priv-x,2);
      if(dist_priv<=r*r){
        ijk_priv=(i_priv+nc)%nc+nc*((j_priv+nc)%nc)+nc2*((k_priv+nc)%nc);
	aux_priv=(maxima-rho[ijk_priv])/(maxima*contador-suma);
	rho[ijk_priv]=aux_priv;
        fuera[ijk_priv]=0;
      }
    }
    }
  }

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)	//	loop sobre las celdas para nomalizar la densidad
    if(fuera[i]==1)
      rho[i]=0;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)	//	loop sobre las celdas para nomalizar la densidad
    rho[i]/=NumVoids;

  printf("calcula la densidad...echo\n");
  fflush(stdout);
}










void densidad_bias(const char *NomArch, int nbin, double dr, int cuantos){
printf("calcula la densidad para el bin %d del bias...\n",nbin);
fflush(stdout);
  float a0, a1, a2, a3, a4, a5, a6;
  float dx, dy, dz, tx, ty, tz;
  int i,j,k,ii,jj,kk,n;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]=0.;		//	no importa si es 0 o -1, la FFT desecha la parte constante

  double aux1 = 1./log(c);
  double aux2_priv;
  int bin_priv;

  FILE * IN;
  IN=fopen(NomArch,"r");
  n=0;	//	número de partículas en el catálogo (para calcular la sobredensidad)

  if(cuantos==5)
  while(fscanf(IN,"%f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5)!=EOF){	// coordenadas físicas
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2_priv=aux1*log(1.*a4/cmin);
      bin_priv=floor(aux2_priv);	//	es el bin de esta particula
      if(bin_priv==nbin){
        n++;			//	es el número de partículas en el bin
        a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
        a1=fmod(1.*a4+dr,1.*nc);a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
        i=floor(a1);		j=floor(a2);			k=floor(a3);
        dx=a1-1.*i;		dy=a2-1.*j;			dz=a3-1.*k;
        tx=1.-dx;		ty=1.-dy;			tz=1.-dz;
        //condiciones de frontera periodicas
        ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
        //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
        rho[i+nc*j +nc2*k ]+=1.*tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=1.*dx*ty*tz;
        rho[i+nc*jj+nc2*k ]+=1.*tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=1.*dx*dy*tz;
        rho[i+nc*j +nc2*kk]+=1.*tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=1.*dx*ty*dz;
        rho[i+nc*jj+nc2*kk]+=1.*tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=1.*dx*dy*dz;
      }
    }
  }

  else if(cuantos==4)
  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){	// coordenadas físicas
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2_priv=aux1*log(1.*a4/cmin);
      bin_priv=floor(aux2_priv);	//	es el bin de esta particula
      if(bin_priv==nbin){
        n++;			//	es el número de partículas en el bin
        a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
        a1=fmod(1.*a4+dr,1.*nc);a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
        i=floor(a1);		j=floor(a2);			k=floor(a3);
        dx=a1-1.*i;		dy=a2-1.*j;			dz=a3-1.*k;
        tx=1.-dx;		ty=1.-dy;			tz=1.-dz;
        //condiciones de frontera periodicas
        ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
        //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
        rho[i+nc*j +nc2*k ]+=1.*tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=1.*dx*ty*tz;
        rho[i+nc*jj+nc2*k ]+=1.*tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=1.*dx*dy*tz;
        rho[i+nc*j +nc2*kk]+=1.*tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=1.*dx*ty*dz;
        rho[i+nc*jj+nc2*kk]+=1.*tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=1.*dx*dy*dz;
      }
    }
  }
  fclose(IN);

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]/=n;		//	ahora si es rho/delta rho

printf("calcula la densidad para el bin %d del bias...echo\n",nbin);
fflush(stdout);
}






















void densidad_bias_halos(const char *NomArch, int nbin, double dr){
printf("calcula la densidad para el bin %d del bias...\n",nbin);
fflush(stdout);
  float a0, a1, a2, a3, a4, a5, a6;
  float dx, dy, dz, tx, ty, tz;
  double aux_m;
  int i,j,k,ii,jj,kk,n;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]=0.;		//	no importa si es 0 o -1, la FFT desecha la parte constante

  double aux1 = 1./log(c);
  double aux2_priv;
  int bin_priv;

  FILE * IN;
  IN=fopen(NomArch,"r");
  n=0;	//	número de partículas en el catálogo (para calcular la sobredensidad)

  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){	// coordenadas físicas
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2_priv=aux1*log(1.*a4/cmin);
      bin_priv=floor(aux2_priv);	//	es el bin de esta particula
      if(bin_priv==nbin){
        a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
        a1=fmod(1.*a4+dr,1.*nc);a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
        i=floor(a1);		j=floor(a2);			k=floor(a3);
        dx=a1-1.*i;		dy=a2-1.*j;			dz=a3-1.*k;
        tx=1.-dx;		ty=1.-dy;			tz=1.-dz;
        aux_m=Mp/a4;
        //condiciones de frontera periodicas
        ii=(i+1)%nc;		jj=(j+1)%nc;		kk=(k+1)%nc;
        //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
        rho[i+nc*j +nc2*k ]+=aux_m*tx*ty*tz;	rho[ii+nc*j +nc2*k ]+=aux_m*dx*ty*tz;
        rho[i+nc*jj+nc2*k ]+=aux_m*tx*dy*tz;	rho[ii+nc*jj+nc2*k ]+=aux_m*dx*dy*tz;
        rho[i+nc*j +nc2*kk]+=aux_m*tx*ty*dz;	rho[ii+nc*j +nc2*kk]+=aux_m*dx*ty*dz;
        rho[i+nc*jj+nc2*kk]+=aux_m*tx*dy*dz;	rho[ii+nc*jj+nc2*kk]+=aux_m*dx*dy*dz;
      }
    }
  }
  fclose(IN);

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)//loop sobre las celdas para borrar la densidad antigua
    rho[i]/=cantidad[nbin];		//	ahora si es rho/delta rho

printf("calcula la densidad para el bin %d del bias...echo\n",nbin);
fflush(stdout);
}

























void densidad_bias_voids(const char *NomArch, int nbin, double dr){
  printf("calcula la densidad...\n");
  fflush(stdout);
  float a0, a1, a2, a3, a4, a5, a6;
  double dx, dy, dz, tx, ty, tz;
  int i,j,k,ii,jj,kk,iii,jjj,kkk,p,q;

  double aux1 = 1./log(c);
  double aux2_priv;
  int bin_priv;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++){//loop sobre las celdas para borrar la densidad antigua
    rho[i]=0.;
    fuera[i]=1;}

  FILE * IN;
  IN=fopen(NomArch,"r");
  //	contribución de la matéria en los voids
  while(fscanf(IN,"%f %f %f %f\n",&a1,&a2,&a3,&a4)!=EOF){	// coordenadas físicas
    if((a4>cmin)&&(a4<cmax)){	//	filtro: solo halos con masa en el rango requerido
      aux2_priv=aux1*log(R[q]*Lbox/cmin/nc);
      bin_priv=floor(aux2_priv);					// es el bin de esta particula  
      if(bin_priv==nbin){
        a4=a1*nc/Lbox;		a5=a2*nc/Lbox;			a6=a3*nc/Lbox;
        a1=fmod(1.*a4+dr,1.*nc);	a2=fmod(1.*a5+dr,1.*nc);	a3=fmod(1.*a6+dr,1.*nc);
        i=floor(1.*a1);		j=floor(1.*a2);			k=floor(1.*a3);
        dx=1.*a1-1.*i;		dy=1.*a2-1.*j;			dz=1.*a3-1.*k;
        tx=1.-dx;			ty=1.-dy;			tz=1.-dz;
        //condiciones de frontera periodicas
        ii=(i+1)%nc;		jj=(j+1)%nc;			kk=(k+1)%nc;
        //incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
        rho[i+nc*j +nc2*k ]+=tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=dx*ty*tz;
        rho[i+nc*jj+nc2*k ]+=tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=dx*dy*tz;
        rho[i+nc*j +nc2*kk]+=tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=dx*ty*dz;
        rho[i+nc*jj+nc2*kk]+=tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=dx*dy*dz;}}}
  fclose(IN);

  //	borrando los bordes de la densidad en los voids
  //	e 'inviertiendo' la densidad en los voids (dens=0 -> delta máximo, dens=maximo -> delta=0)

  for(q=0;q<NumVoids;q++){
    if(((R[q]*Lbox/nc)>cmin)&&((R[q]*Lbox/nc)<cmax)){//	filtro: solo halos con masa en el rango requerido
    aux2_priv=aux1*log(R[q]*Lbox/cmin/nc);
    bin_priv=floor(aux2_priv);	//	es el bin de esta particula  
    if(bin_priv==nbin){
      double x,y,z,r;
      double maxima=0.,suma=0.;
      int contador=0;
      x=X[q]+dr;
      y=Y[q]+dr;
      z=Z[q]+dr;
      r=R[q];
      i=floor(x);
      j=floor(y);
      k=floor(z);
      int max=ceil(R[q]);    

      #pragma omp parallel
      {
      int    i_priv;
      int    j_priv;
      int    k_priv;
      int    ijk_priv;
      int    cont_priv=0;
      double dist_priv;
      double max_priv=0.;
      double sum_priv=0.;
      double aux_priv;
      #pragma omp for
      for(p=0;p<(2*max+1)*(2*max+1)*(2*max+1);p++){
        k_priv=k+fk(p,max);
        j_priv=j+fj(p,max);
        i_priv=i+fi(p,max);
        dist_priv=pow(1.*k_priv-z,2)+pow(1.*j_priv-y,2)+pow(1.*i_priv-x,2);
        if(dist_priv<=r*r){
          ijk_priv=(i_priv+nc)%nc+nc*((j_priv+nc)%nc)+nc2*((k_priv+nc)%nc);
	  sum_priv+=rho[ijk_priv];
	  cont_priv++;
	  if(max_priv<rho[ijk_priv])
	    max_priv=rho[ijk_priv];
        }
      }

      #pragma critical
      {
      if(maxima<max_priv)
        maxima=max_priv;
      suma+=sum_priv;
      contador+=cont_priv;
      }

      #pragma omp for
      for(p=0;p<(2*max+1)*(2*max+1)*(2*max+1);p++){
        k_priv=k+fk(p,max);
        j_priv=j+fj(p,max);
        i_priv=i+fi(p,max);
        dist_priv=pow(1.*k_priv-z,2)+pow(1.*j_priv-y,2)+pow(1.*i_priv-x,2);
        if(dist_priv<=r*r){
          ijk_priv=(i_priv+nc)%nc+nc*((j_priv+nc)%nc)+nc2*((k_priv+nc)%nc);
	  aux_priv=(maxima-rho[ijk_priv])/(maxima*contador-suma);
	  rho[ijk_priv]=aux_priv;
          fuera[ijk_priv]=0;
        }
      }
      }
    }
  }
  }

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)	//	loop sobre las celdas para nomalizar la densidad
    if(fuera[i]==1)
      rho[i]=0;

  #pragma omp parallel for
  for(i=0;i<NumCel;i++)	//	loop sobre las celdas para nomalizar la densidad
    rho[i]/=cantidad[nbin];

  printf("calcula la densidad...echo\n");
  fflush(stdout);
}























void densidad_distorcion(const char *NomArch,int eje){
  float a0, a1, a2, a3, a4, a5, a6;
  float dx, dy, dz, tx, ty, tz;
  int i,j,k,ii,jj,kk,iii,jjj,kkk,ijk,n,p;

  //loop sobre las celdas para borrar la densidad y velocidad antigua
  #pragma omp parallel for
  for(i=0;i<NumCel;i++){
    x[i]=0.;
    vx[i]=0.;}


  FILE * IN;
  IN=fopen(NomArch,"r");
  //loop sobre las particulas para calcular la densidad de las celdas
  for(p=0;p<NumPart;p++){
    fscanf(IN,"%f %f %f %f %f %f\n",&a1,&a2,&a3,&a4,&a5,&a6);
    a1*=nc/Lbox;		a2*=nc/Lbox;		a3*=nc/Lbox;
    i=floor(a1);		j=floor(a2);		k=floor(a3);
    dx=a1-1.*i;			dy=a2-1.*j;		dz=a3-1.*k;
    tx=1.-dx;			ty=1.-dy;		tz=1.-dz;
    //condiciones de frontera periodicas
    ii=(i+1)%nc;		jj=(j+1)%nc;		kk=(k+1)%nc;
    //incrementando la densidade de las celdas adyacentes
    x[i+nc*j +nc2*k ]+=Mp*tx*ty*tz;		x[ii+nc*j +nc2*k ]+=Mp*dx*ty*tz;
    x[i+nc*jj+nc2*k ]+=Mp*tx*dy*tz;		x[ii+nc*jj+nc2*k ]+=Mp*dx*dy*tz;
    x[i+nc*j +nc2*kk]+=Mp*tx*ty*dz;		x[ii+nc*j +nc2*kk]+=Mp*dx*ty*dz;
    x[i+nc*jj+nc2*kk]+=Mp*tx*dy*dz;		x[ii+nc*jj+nc2*kk]+=Mp*dx*dy*dz;

    //incrementando la velocidad de las celdas adyacentes
    if(eje==0) a0=a4;
    if(eje==1) a0=a5;
    if(eje==2) a0=a6;
    vx[i+nc*j +nc2*k ]+=a0*tx*ty*tz;		vx[ii+nc*j +nc2*k ]+=a0*dx*ty*tz;
    vx[i+nc*jj+nc2*k ]+=a0*tx*dy*tz;		vx[ii+nc*jj+nc2*k ]+=a0*dx*dy*tz;
    vx[i+nc*j +nc2*kk]+=a0*tx*ty*dz;		vx[ii+nc*j +nc2*kk]+=a0*dx*ty*dz;
    vx[i+nc*jj+nc2*kk]+=a0*tx*dy*dz;		vx[ii+nc*jj+nc2*kk]+=a0*dx*dy*dz;}
  fclose(IN);

  /////////////////////////		calculo de la densidad de los pixeles pequeños
  ijk=0;
  for(kkk=0;kkk<nc;kkk++)
    for(jjj=0;jjj<nc;jjj++)
      for(iii=0;iii<nc;iii++){//loop sobre las celdas pequeñas
	a1=1.*iii;		a2=1.*jjj;	a3=1.*kkk;
	if((x[ijk]!=0)){
	  if(eje==0){
	    a4=a1+vx[ijk]/x[ijk];
	    a1=fmod(a4+1.*nc,nc);
	    if((a1>nc) || (a1<0)) cout<<"error x="<<a4<<endl;}
	  if(eje==1){
	    a4=a2+vx[ijk]/x[ijk];
	    a2=fmod(a4+1.*nc,nc);
	    if((a2>nc) || (a2<0)) cout<<"error y="<<a4<<endl;}
	  if(eje==2){
	    a4=a3+vx[ijk]/x[ijk];
	    a3=fmod(a4+1.*nc,nc);
	    if((a3>nc) || (a3<0)) cout<<"error z="<<a4<<endl;}}
	i=floor(a1);		j=floor(a2);		k=floor(a3);
	dx=a1-1.*i;		dy=a2-1.*j;		dz=a3-1.*k;
	tx=1.-dx;		ty=1.-dy;		tz=1.-dz;
	//condiciones de frontera periodicas
	ii=(i+1)%nc;		jj=(j+1)%nc;		kk=(k+1)%nc;
	//incrementando la densidade de las celdas adyacentes	(solo las componentes reales)
	a1=x[ijk];
        rho[i+nc*j +nc2*k ]+=a1*tx*ty*tz;		rho[ii+nc*j +nc2*k ]+=a1*dx*ty*tz;
        rho[i+nc*jj+nc2*k ]+=a1*tx*dy*tz;		rho[ii+nc*jj+nc2*k ]+=a1*dx*dy*tz;
        rho[i+nc*j +nc2*kk]+=a1*tx*ty*dz;		rho[ii+nc*j +nc2*kk]+=a1*dx*ty*dz;
        rho[i+nc*jj+nc2*kk]+=a1*tx*dy*dz;		rho[ii+nc*jj+nc2*kk]+=a1*dx*dy*dz;
	ijk++;}
}

















void media(const char *NomArch){
  printf("calcula la media del espectro...\n");
  fflush(stdout);
  int i,j,k,ijk,j2,k2;
  double elk;

  ///////////////////////         definicion de los (bin-1)  bines en k
  double ck   =pow(kmax/kmin,1./bin);	//		  MPb2[id]+=aux_h;constante de proporcionalidad log
  double aux1 =1./log(ck);
  double aux2;

  //////////////////////		borro la información antigua
  for(i=0;i<bin;i++){
    MP[i]=0.;
    Mk[i]=0.;
    EP[i]=0.;
    Ek[i]=0.;
    cont[i]=0;
  }

  #pragma omp parallel
  {
  double auxa;
  #pragma omp for
  for(ijk=0;ijk<nc21nc*nc;ijk++){
    auxa=0.5*(delta[ijk][0]+deltac[ijk][0]);
    delta[ijk][0]=auxa;
    auxa=0.5*(delta[ijk][1]+deltac[ijk][1]);
    delta[ijk][1]=auxa;
  }
  }

  //////////////////////		lleno los bines, contador y radio
  #pragma omp parallel
  {
  int id,ia,ja,ka,i2a,j2a,k2a;
  double aux2a,elka;
  double MPa[bin]={0.},Mka[bin]={0.},conta[bin]={0.};
  #pragma omp for
  for(ijk=0;ijk<nc21nc*nc;ijk++){
    ia=ijk%nc21;
    ja=(ijk-ia)%nc21nc;
    ka=(ijk-ia-ja)/nc21nc;
    ja/=nc21;
    if((ja==0) || (ja==nc/2))	{j2a=ja;}	else	{j2a=nc-ja;}
    if(ja<j2a)			{j2a=ja;}	else	;
    if((ka==0) || (ka==nc/2))	{k2a=ka;}	else	{k2a=nc-ka;}
    if(ka<k2a)			{k2a=ka;}	else	;
    elka=dpin*sqrt(k2a*k2a+j2a*j2a+ia*ia);		//	en h/Mpc, elk=2pi/Lbox
    if((elka>kmin)&&(elka<kmax)){
      if(bin_lineal==0)
        aux2a=aux1*log(elka/kmin);
      else
        aux2a=bin*(elka-kmin)/(kmax-kmin);
      id=floor(aux2a);         	//	es el bin al que pertenece elk
      conta[id]++;             	//	contador de k's	  
      MPa[id]+=kernel(elka)*(delta[ijk][0]*delta[ijk][0]+delta[ijk][1]*delta[ijk][1]);
      Mka[id]+=elka;
    }
  }
  #pragma omp critical
  {
  for(i=0;i<bin;i++){
    MP[i]+=MPa[i];
    Mk[i]+=Mka[i];
    cont[i]+=conta[i];
  }
  }
  }	//	pragma


  //////////////////////		divido entre cont para obtener la média
  for(i=0;i<bin;i++){
    MP[i]/=cont[i];
    Mk[i]/=cont[i];
  }

  //////////////////////          calculo la desviación
  #pragma omp parallel
  {
  int id,ia,ja,ka,i2a,j2a,k2a;
  double aux2a,elka;
  double EPa[bin]={0.},Eka[bin]={0.};
  #pragma omp for
  for(ijk=0;ijk<nc21nc*nc;ijk++){
    ia=ijk%nc21;
    ja=(ijk-ia)%nc21nc;
    ka=(ijk-ia-ja)/nc21nc;
    ja/=nc21;
    if((ja==0) || (ja==nc/2))	{j2a=ja;}	else	{j2a=nc-ja;}
    if(ja<j2a)			{j2a=ja;}	else	;
    if((ka==0) || (ka==nc/2))	{k2a=ka;}	else	{k2a=nc-ka;}
    if(ka<k2a)			{k2a=ka;}	else	;
    elka=dpin*sqrt(k2a*k2a+j2a*j2a+ia*ia);		//	en h/Mpc, elk=2pi/Lbox
    if((elka>kmin)&&(elka<kmax)){
      if(bin_lineal==0)
        aux2a=aux1*log(elka/kmin);
      else
        aux2a=bin*(elka-kmin)/(kmax-kmin);
      id=floor(aux2a);         	//	es el bin al que pertenece elk
      EPa[id]+=pow(kernel(elka)*(delta[ijk][0]*delta[ijk][0]+delta[ijk][1]*delta[ijk][1])-MP[id],2.);
      Eka[id]+=pow(elka-Mk[id],2.);
    }
  }
  #pragma omp critical
  {
  for(i=0;i<bin;i++){
    EP[i]+=EPa[i];
    Ek[i]+=Eka[i];
  }
  }
  }	//	pragma


  //////////////////////		normalizo la desviación
  for(i=0;i<bin;i++){
    double aux_p=EP[i]/(cont[i]*(cont[i]-1));
    EP[i]=sqrt(aux_p);
    double aux_k=Ek[i]/(cont[i]*(cont[i]-1));
    Ek[i]=sqrt(aux_k);
  }

  //////////////////////		imprimo
  FILE *ES;
  ES=fopen(NomArch,"w+");
  for(i=0;i<bin;i++){
    fprintf(ES,"%le %le %le %le\n",Mk[i],MP[i],Ek[i],EP[i]);
  }
  fclose(ES);
  printf("calcula la media del espectro...echo\n");
  fflush(stdout);
}



















void media_bias(const char *NomArch,int nbin){
  printf("calcula el bias para el bin %d...\n",nbin);
  fflush(stdout);
  int i,j,k,ijk,j2,k2;
  double elk;

  ///////////////////////         definicion de los (bin-1)  bines en k
  double ck   =pow(kmax/kmin,1./bin);	//	constante de proporcionalidad log
  double aux1 =1./log(ck);
  double aux2;

  //////////////////////		media de los dos espectros
  #pragma omp parallel
  {
  double auxa;
  #pragma omp for
  for(ijk=0;ijk<nc21nc*nc;ijk++){
    auxa=0.5*(deltab[ijk][0]+deltac[ijk][0]);
    deltab[ijk][0]=auxa;
    auxa=0.5*(deltab[ijk][1]+deltac[ijk][1]);
    deltab[ijk][1]=auxa;
  }
  }	//	pragma

  printf("ya junto las dos muestras, ahora va a repartir entre los bines\n");
  fflush(stdout);

  //////////////////////		borro la información antigua
  for(i=0;i<bin;i++){
    MP[i]=0.;      EP[i]=0.;
    MP2[i]=0.;     EP2[i]=0.;
    MPb[i]=0.;     EPb[i]=0.;
    MPb2[i]=0.;    EPb2[i]=0.;
    MPb3[i]=0.;    EPb3[i]=0.;
    MPb4[i]=0.;    EPb4[i]=0.;
    MPb5[i]=0.;    EPb5[i]=0.;
    MPb6[i]=0.;    EPb6[i]=0.;
    cont[i]=0;     cont2[i]=0;
  }

  //////////////////////          lleno los bines, contador y radio
  #pragma omp parallel
  {
  int id,ia,ja,ka,i2a,j2a,k2a;
  double aux_h,aux_m,aux_hm,elka;
  double aux2a;
  double MPa[bin]={0.},MP2a[bin]={0.},MPba[bin]={0.},MPb2a[bin]={0.},MPb3a[bin]={0.},MPb4a[bin]={0.};
  double MPb5a[bin]={0.},MPb6a[bin]={0.},conta[bin]={0.},cont2a[bin]={0.};
  #pragma omp for
  for(ijk=0;ijk<nc21nc*nc;ijk++){
    ia=ijk%nc21;
    ja=(ijk-ia)%nc21nc;
    ka=(ijk-ia-ja)/nc21nc;
    ja/=nc21;
    if((ja==0) || (ja==nc/2))	{j2a=ja;}	else	{j2a=nc-ja;}
    if(ja<j2a)			{j2a=ja;}	else	;
    if((ka==0) || (ka==nc/2))	{k2a=ka;}	else	{k2a=nc-ka;}
    if(ka<k2a)			{k2a=ka;}	else	;
    elka=dpin*sqrt(k2a*k2a+j2a*j2a+ia*ia);		//	en h/Mpc, elk=2pi/Lbox
    if((elka>kmin)&&(elka<kmax)){
      if(bin_lineal==0)
        aux2a=aux1*log(elka/kmin);
      else
        aux2a=bin*(elka-kmin)/(kmax-kmin);
      id=floor(aux2a);         	//	es el bin al que pertenece elk
      conta[id]++;             	//	contador de k's	  
      aux_m=            pow(delta[ijk][0],2.) +pow(delta[ijk][1],2.);
      aux_h=           pow(deltab[ijk][0],2.) +pow(deltab[ijk][1],2.) -1./cantidad[nbin];
      aux_hm=    deltab[ijk][0]*delta[ijk][0] +deltab[ijk][1]*delta[ijk][1];
      MPa[id]+=aux_m;	//	<dm*dm>
      MPba[id]+=aux_hm;	//	<dm*dv+dv*dm>/2
      MPb2a[id]+=aux_h;	//	<dv*dv>
      if(aux_m!=0.){
        cont2a[id]++;
	MP2a[id]+=sqrt(aux_m);		//	<sqrt(dm*dm)>
	MPb3a[id]+=aux_hm/aux_m;	//	<(dm*dv+dv*dm)/(2 dm*dm)>
	MPb4a[id]+=aux_h/aux_m;		//	<dv*dv/dm*dm>
	MPb5a[id]+=sqrt(aux_h/aux_m);	//	<sqrt(dv*dv/dm*dm)>
	MPb6a[id]+=sqrt(aux_h);		//	<sqrt(dv*dv)>
      }
    }
  }
  #pragma omp critical
  {
  for(i=0;i<bin;i++){
    MP[i]+=MPa[i];
    MP2[i]+=MP2a[i];
    MPb[i]+=MPba[i];
    MPb2[i]+=MPb2a[i];
    MPb3[i]+=MPb3a[i];
    MPb4[i]+=MPb4a[i];
    MPb5[i]+=MPb5a[i];
    MPb6[i]+=MPb6a[i];
    cont[i]+=conta[i];
    cont2[i]+=cont2a[i];
  }
  }
  }	//	pragma

  printf("ahora va a calcular la dispersion\n");
  fflush(stdout);

  //////////////////////		divido entre cont para obtener la média
  for(i=0;i<bin;i++){
    MP[i]/=cont[i];
    MPb[i]/=cont[i];
    MPb2[i]/=cont[i];
    MP2[i]/=cont2[i];
    MPb3[i]/=cont2[i];
    MPb4[i]/=cont2[i];
    MPb5[i]/=cont2[i];
    MPb6[i]/=cont2[i];
  }
  fflush(stdout);

  //////////////////////          calculo la desviación
 #pragma omp parallel
  {
  int id,ia,ja,ka,i2a,j2a,k2a;
  double aux_h,aux_m,aux_hm,elka;
  double aux2a;
  double EPa[bin]={0.},EP2a[bin]={0.},EPba[bin]={0.},EPb2a[bin]={0.},EPb3a[bin]={0.},EPb4a[bin]={0.};
  double EPb5a[bin]={0.},EPb6a[bin]={0.};
  #pragma omp for
  for(ijk=0;ijk<nc21nc*nc;ijk++){
    ia=ijk%nc21;
    ja=(ijk-ia)%nc21nc;
    ka=(ijk-ia-ja)/nc21nc;
    ja/=nc21;
    if((ja==0) || (ja==nc/2))	{j2a=ja;}	else	{j2a=nc-ja;}
    if(ja<j2a)			{j2a=ja;}	else	;
    if((ka==0) || (ka==nc/2))	{k2a=ka;}	else	{k2a=nc-ka;}
    if(ka<k2a)			{k2a=ka;}	else	;
    elka=dpin*sqrt(k2a*k2a+j2a*j2a+ia*ia);		//	en h/Mpc, elk=2pi/Lbox
    if((elka>kmin)&&(elka<kmax)){
      if(bin_lineal==0)
        aux2a=aux1*log(elka/kmin);
      else
        aux2a=bin*(elka-kmin)/(kmax-kmin);
      id=floor(aux2a);         	//	es el bin al que pertenece elk
      aux_m=            pow(delta[ijk][0],2.) +pow(delta[ijk][1],2.);
      aux_h=           pow(deltab[ijk][0],2.) +pow(deltab[ijk][1],2.) -1./cantidad[nbin];
      aux_hm=    deltab[ijk][0]*delta[ijk][0] +deltab[ijk][1]*delta[ijk][1];
      EPa[id]+=pow(aux_m-MP[id],2.);
      EPba[id]+=pow(aux_hm-MPb[id],2.);
      EPb2a[id]+=pow(aux_h-MPb2[id],2.);
      if(aux_m!=0.){
	EP2a[id]+=pow(sqrt(aux_m)-MP2[id],2.);
	EPb3a[id]+=pow(aux_hm/aux_m-MPb3[id],2.);
	EPb4a[id]+=pow(aux_h/aux_m-MPb4[id],2.);
	EPb5a[id]+=pow(sqrt(aux_h/aux_m)-MPb5[id],2.);
	EPb6a[id]+=pow(sqrt(aux_h)-MPb6[id],2.);
      }
    }
  }
  #pragma omp critical
  {
  for(i=0;i<bin;i++){
    EP[i]+=EPa[i];
    EP2[i]+=EP2a[i];
    EPb[i]+=EPba[i];
    EPb2[i]+=EPb2a[i];
    EPb3[i]+=EPb3a[i];
    EPb4[i]+=EPb4a[i];
    EPb5[i]+=EPb5a[i];
    EPb6[i]+=EPb6a[i];
  }
  }
  }

  //////////////////////		normalizo la desviación
  for(i=0;i<bin;i++){
    double aux_p;
    aux_p=EP[i]/(cont[i]*(cont[i]-1));		EP[i]=sqrt(aux_p);
    aux_p=EPb[i]/(cont[i]*(cont[i]-1));		EPb[i]=sqrt(aux_p);
    aux_p=EPb2[i]/(cont[i]*(cont[i]-1));	EPb2[i]=sqrt(aux_p);

    aux_p=EP2[i]/(cont2[i]*(cont2[i]-1));	EP2[i]=sqrt(aux_p);
    aux_p=EPb3[i]/(cont2[i]*(cont2[i]-1));	EPb3[i]=sqrt(aux_p);
    aux_p=EPb4[i]/(cont2[i]*(cont2[i]-1));	EPb4[i]=sqrt(aux_p);
    aux_p=EPb5[i]/(cont2[i]*(cont2[i]-1));	EPb5[i]=sqrt(aux_p);
    aux_p=EPb6[i]/(cont2[i]*(cont2[i]-1));	EPb6[i]=sqrt(aux_p);
/*
    if(EP[i]<MP[i]/sqrt(cantidad[nbin]))	EP[i]=MP[i]/sqrt(cantidad[nbin]);
    if(EPb[i]<MPb[i]/sqrt(cantidad[nbin]))	EPb[i]=MPb[i]/sqrt(cantidad[nbin]);
    if(EPb2[i]<MPb2[i]/sqrt(cantidad[nbin]))	EPb2[i]=MPb2[i]/sqrt(cantidad[nbin]);

    if(EP2[i]<MP2[i]/sqrt(cantidad[nbin]))	EP2[i]=MP2[i]/sqrt(cantidad[nbin]);
    if(EPb3[i]<MPb3[i]/sqrt(cantidad[nbin]))	EPb3[i]=MPb3[i]/sqrt(cantidad[nbin]);
    if(EPb4[i]<MPb4[i]/sqrt(cantidad[nbin]))	EPb4[i]=MPb4[i]/sqrt(cantidad[nbin]);
    if(EPb5[i]<MPb5[i]/sqrt(cantidad[nbin]))	EPb5[i]=MPb5[i]/sqrt(cantidad[nbin]);
    if(EPb6[i]<MPb6[i]/sqrt(cantidad[nbin]))	EPb6[i]=MPb6[i]/sqrt(cantidad[nbin]);
*/
  }
  fflush(stdout);




  printf("ahora imprime las medidas\n");
  fflush(stdout);

  //////////////////////		imprimo medidas
  char nombre[400];
  FILE *ES;
  if(bin_lineal==0)
    sprintf(nombre,"medidas_%02d_%s",nbin,NomArch);
  else
    sprintf(nombre,"medidas_lin_%02d_%s",nbin,NomArch);
  ES=fopen(nombre,"w+");

  FILE *EZ;
  sprintf(nombre,"espectros_%02d_%s",nbin,NomArch);
  EZ=fopen(nombre,"w+");

  double bias,ebias;

  for(i=0;i<bin;i++){
    //	K Pmm Phm Phh
    fprintf(EZ,"%f %f %e %e %e %e %e %e\n",Mk[i],Ek[i],MP[i],EP[i],MPb[i],EPb[i],MPb2[i],EPb2[i]);
    fprintf(ES,"%f %f ",Mk[i],Ek[i]);

    bias=MPb[i]/MP[i];	//	<dm*dv+dv*dm>/2   /  <dm*dm>
    ebias=abs(bias)*sqrt(pow(EPb[i]/MPb[i],2.)+pow(EP[i]/MP[i],2.));
    fprintf(ES,"%e %e ",bias,ebias);

    bias=sqrt(MPb2[i]/MP[i]);//	sqrt(  <dv*dv>   /  <dm*dm>  )
    ebias=abs(bias)*sqrt(pow(EPb2[i]/MPb2[i],2.)+pow(EP[i]/MP[i],2.))*0.5;
    fprintf(ES,"%e %e ",bias,ebias);

//    bias=MPb6[i]/MP2[i];	//	<sqrt(dv*dv)>   /  <sqrt(dm*dm)>
//    ebias=bias*sqrt(pow(EPb6[i]/MPb6[i],2.)+pow(EP2[i]/MP2[i],2.));
//    fprintf(ES,"%e %e ",bias,ebias);

    bias=MPb3[i];		//	<(dm*dv+dv*dm)/(2 dm*dm)>
    ebias=EPb3[i];
    fprintf(ES,"%e %e ",bias,ebias);

    bias=sqrt(MPb4[i]);		//	sqrt( < dv*dv / dm*dm > )
    ebias=EPb4[i]/(0.5*bias);
    fprintf(ES,"%e %e ",bias,ebias);
	
//    bias=MPb5[i];		//	<sqrt(dv*dv/dm*dm)>
//    ebias=EPb5[i];
//    fprintf(ES,"%e %e ",bias,ebias);
    fprintf(ES,"\n");
  }
  fclose(ES);
  fclose(EZ);




  printf("ahora va a calcular e imprimir el bias\n");
  fflush(stdout);

  ////////////////////			calculando el bias lineal
  j=0;
  int b;
  for(b=0;b<bin;b++)
    if(Mk[b]<k_lineal)
      j++;
  i=0;
  for(b=0;b<bin;b++)
    if((Mk[b]>0.0)&&(Mk[b]<k_lineal))
      i++;
  j-=i;

  printf("j=%d i=%d Mk[j]=%e Mk[i]=%e\n",j,i,Mk[j],Mk[i]);  fflush(stdout);

  VecDoub K(i);
  VecDoub K2(i);
  VecDoub B(i);
  VecDoub S(i);

  printf("antes de interpolar\n");  fflush(stdout);

  //	imprime ajuste
  FILE *BI;
  BI=fopen(NomArch,"a");
  fprintf(BI,"%f %f ",carga[nbin],Ecarga[nbin]);

  for(b=0;b<i;b++){
    K[b]=Mk[b+j];
    K2[b]=K[b]*K[b];
    B[b]=MPb[b+j]/MP[b+j];	//	<dm*dv+dv*dm>/2   /  <dm*dm>
    S[b]=B[b]*sqrt(pow(EPb[b+j]/MPb[b+j],2.)+pow(EP[b+j]/MP[b+j],2.));}
  printf("antes de interpolar\n");  fflush(stdout);
  Fitab linea1(K2,B,S);
  fprintf(BI,"%e %e ",linea1.a,linea1.siga);
  printf("primer bias\n");  fflush(stdout);

  for(b=0;b<i;b++){
    if(MPb2[b+j]>0){
      B[b]=sqrt(MPb2[b+j]/MP[b+j]);//	sqrt(  <dv*dv>   /  <dm*dm>  )
      S[b]=B[b]*sqrt(pow(EPb2[b+j]/MPb2[b+j],2.)+pow(EP[b+j]/MP[b+j],2.))*0.5;}
    else{
      B[b]=0.;
      S[b]=sqrt(1./(cantidad[nbin]*MP[b+j]));}}
  printf("antes de interpolar\n");  fflush(stdout);
  Fitab linea2(K2,B,S);
  fprintf(BI,"%e %e ",linea2.a,linea2.siga);
  printf("segundo bias\n");  fflush(stdout);

/*  for(b=0;b<i;b++){
    B[b]=MPb6[b+j]/MP2[b+j];	//	<sqrt(dv*dv)>   /  <sqrt(dm*dm)>
    S[b]=B[b]*sqrt(pow(EPb6[b+j]/MPb6[b+j],2.)+pow(EP2[b+j]/MP2[b+j],2.));}
  printf("antes de interpolar\n");  fflush(stdout);
  Fitab linea3(K,B,S);
  fprintf(BI,"%e %e ",linea3.a,linea3.siga);
  printf("tercer bias\n");  fflush(stdout);

  for(b=0;b<i;b++){
    B[b]=MPb3[b+j];		//	<(dm*dv+dv*dm)/(2 dm*dm)>
    S[b]=EPb3[b+j];}
  printf("antes de interpolar\n");  fflush(stdout);
  Fitab linea4(K,B,S);
  fprintf(BI,"%e %e ",linea4.a,linea4.siga);
  printf("cuarto bias\n");  fflush(stdout);

  for(b=0;b<i;b++){
    B[b]=sqrt(MPb4[b+j]);	//	sqrt( < dv*dv / dm*dm > )
    S[b]=EPb4[b+j]/(0.5*B[b]);}
  printf("antes de interpolar\n");  fflush(stdout);
  Fitab linea5(K,B,S);
  fprintf(BI,"%e %e ",linea5.a,linea5.siga);
  printf("quinto bias\n");  fflush(stdout);
	
  for(b=0;b<i;b++){
    B[b]=MPb5[b+j];		//	<sqrt(dv*dv/dm*dm)>
    S[b]=EPb5[b+j];}
  printf("antes de interpolar\n");  fflush(stdout);
  Fitab linea6(K,B,S);
  fprintf(BI,"%e %e\n",linea6.a,linea6.siga);
  printf("sexto bias\n");  fflush(stdout);
*/
  fprintf(BI,"\n");
  fclose(BI);




  printf("ahora imprime los ajustes\n");
  fflush(stdout);

  //////////////////////		imprimo ajuste
  sprintf(nombre,"ajuste_%02d_%s",nbin,NomArch);
  ES=fopen(nombre,"w+");
  double err1,err2,err3,err4,err5,err6;
  double val=0.01;		//	imprime el ajuste lineal
  while(val<k_lineal){
    err1=linea1.a+linea1.b*val*val;
    err2=linea2.a+linea2.b*val*val;
//    err3=linea3.a+linea3.b*val;
//    err4=linea4.a+linea4.b*val;
//    err5=linea5.a+linea5.b*val;
//    err6=linea6.a+linea6.b*val;
    fprintf(ES,"%e %e %e\n",val,err1,err2);
    val*=1.15;
  }
  fclose(ES);

  printf("calcula el bias...echo\n");
  fflush(stdout);
}














double poly(int order, double *coeff, double x){
  double valor=coeff[0];
  for(int o=1;o<order+1;o++)
    valor+=coeff[o]*pow(x,o);
  return valor;
}



void mcmc(int len_data,VecDoub x,VecDoub y, VecDoub sx, int order,double *coeff, double *sigma_coeff){
Crandom ran(92830);
double delta_coeff[order+1];
int iterations=0,n_salta=0;
int max_iterations=100000;
int burn_in=100000;
double chi2,chi_new,chi_old;
int i,j;
FILE *NOM;
char nomAr[200];
sprintf(nomAr,"lista.txt");
NOM=fopen(nomAr,"w+");

chi2=0.;
for(i=0;i<len_data;i++){
  chi2+=pow((y[i]-poly(order,coeff,x[i]))/sx[i],2);
}
chi_old=chi2;
chi_new=chi2;
while (iterations<max_iterations+burn_in){
  for(i=0;i<order+1;i++){
    delta_coeff[i]=ran.gauss(0.,sigma_coeff[i]);
    coeff[i]+=delta_coeff[i];}
  chi2=0.;
  for(i=0;i<len_data;i++){
    chi2+=pow((y[i]-poly(order,coeff,x[i]))/sx[i],2);
  }
  chi_old=chi_new;
  chi_new=chi2;
  if(chi_new<=chi_old){
    if(iterations>burn_in){
      fprintf(NOM,"%d	%le ",n_salta,chi_new);
      for(int k=0;k<order+1;k++)
        fprintf(NOM,"%le ",coeff[k]);
      fprintf(NOM,"\n");
    }
    iterations++;
    n_salta=0;
  }
  else {
    if((exp(0.5*(chi_old-chi_new)))>=ran.r()) {
      if(iterations>burn_in){
        fprintf(NOM,"%d %le ",n_salta,chi_new);
        for(int k=0;k<order+1;k++)
          fprintf(NOM,"%le	",coeff[k]);
      fprintf(NOM,"\n");
    }
      iterations++;
      n_salta=0;
    }
    else{
      for(int k=0;k<order+1;k++)
        coeff[k]-=delta_coeff[k];
      chi_new=chi_old;
    }
  }
}



//media_error
  FILE * ES;
  double a[order+1],pMin[order+1],c2,chi_min;
  int peso,contador;
  double mean[order+1],stdev[order+1];

  for(i=0;i<order+1;i++){			//	analisis
    mean[i]=0.;
    stdev[i]=0.;
  }
  contador=0;

  sprintf(nomAr,"lista.txt");
  ES=fopen(nomAr,"r");

  for(j=0;j<max_iterations;j++){
    fscanf(ES,"%d %le",&peso,&chi_new);
    contador++;//=peso;

    for(i=0;i<order+1;i++){
      fscanf(ES," %le",&a[i]);
      mean[i]+=a[i];//*peso;
    }
    fscanf(ES,"\n");

    if(j==0) chi_min=chi_new;

    if(chi_new<chi_min){
      chi_min=chi_new;				//	chi_minimo
      for(i=0;i<order+1;i++){
        pMin[i]=a[i];			//	parametros
       }
    }
  }

  for(i=0;i<order+1;i++)
    mean[i]/=contador;

  rewind(ES);
  for(j=0;j<max_iterations;j++){
    fscanf(ES,"%d %le",&peso,&c2);
    for(i=0;i<order+1;i++){
      fscanf(ES," %le",&a[i]);
      stdev[i]+=pow(a[i]-mean[i],2);//*peso;
    }
    fscanf(ES,"\n");
  }
  for(i=0;i<order+1;i++)
    stdev[i]/=contador;
  fclose(ES);

corte=pMin[0];
scorte=sqrt(stdev[0]);
pendiente=mean[1];
if(order==2)
segunda=mean[2];

}














//	solo calcula el bias dados los espectros
void calcula_bias(const char *NomArch,int nbin){
  printf("calcula el bias para el bin %d...\n",nbin);
  fflush(stdout);

  double bias,ebias;
  int i,j;

  //////////////////////		imprimo medidas
  char nombre[400];
  FILE *ES;
  if(bin_lineal==0)
    sprintf(nombre,"medidas_%02d_%s",nbin,NomArch);
  else
    sprintf(nombre,"medidas_lin_%02d_%s",nbin,NomArch);
  ES=fopen(nombre,"r");

  for(i=0;i<bin;i++){
    fscanf(ES,"%lf %lf ",&bias,&ebias);    //	K Pmm Phm Phh
    Mk[i]=bias;	Ek[i]=ebias;
    fscanf(ES,"%le %le ",&bias,&ebias);	  //	<dm*dv+dv*dm>/2   /  <dm*dm>
    fscanf(ES,"%le %le ",&bias,&ebias);    //	sqrt(  <dv*dv>   /  <dm*dm>  )
//    fprintf(ES,"%e %e ",bias,ebias);    //	<sqrt(dv*dv)>   /  <sqrt(dm*dm)>
    fscanf(ES,"%le %le ",&bias,&ebias);    //	<(dm*dv+dv*dm)/(2 dm*dm)>
    fscanf(ES,"%le %le ",&bias,&ebias);    //	sqrt( < dv*dv / dm*dm > )
//    fprintf(ES,"%e %e ",bias,ebias);    //	<sqrt(dv*dv/dm*dm)>
  }


  ////////////////////			calculando el bias lineal
  j=0;
  int b;
  for(b=0;b<bin;b++)
    if(Mk[b]<k_lineal)
      j++;
  i=0;
  for(b=0;b<bin;b++)
    if((Mk[b]>0.0)&&(Mk[b]<k_lineal))
      i++;
  j-=i;

//  printf("j=%d i=%d Mk[j]=%e Mk[i]=%e\n",j,i,Mk[j],Mk[i]);  fflush(stdout);

  VecDoub K(i);
  VecDoub K2(i);
  VecDoub B(i);
  VecDoub S(i);

//  printf("antes de interpolar\n");  fflush(stdout);

  rewind(ES);
  for(b=0;b<i;b++){
    fscanf(ES,"%lf %lf ",&bias,&ebias);    //	K Pmm Phm Phh
    K[b]=bias;
    K2[b]=bias*bias;
    fscanf(ES,"%le %le ",&bias,&ebias);	  //	<dm*dv+dv*dm>/2   /  <dm*dm>
    B[b]=bias;	//	<dm*dv+dv*dm>/2   /  <dm*dm>

///////////////////////////////////////////////////////////////////////////////
////////////////////////////	error del espectro
///////////////////////////////////////////////////////////////////////////////
    S[b]=ebias;
//    S[b]=B[b]*sqrt(2./cantidad[b]);
    
    fscanf(ES,"%le %le ",&bias,&ebias);    //	sqrt(  <dv*dv>   /  <dm*dm>  )
//    fprintf(ES,"%e %e ",bias,ebias);    //	<sqrt(dv*dv)>   /  <sqrt(dm*dm)>
    fscanf(ES,"%le %le ",&bias,&ebias);    //	<(dm*dv+dv*dm)/(2 dm*dm)>
    fscanf(ES,"%le %le ",&bias,&ebias);    //	sqrt( < dv*dv / dm*dm > )
//    fprintf(ES,"%e %e ",bias,ebias);    //	<sqrt(dv*dv/dm*dm)>
  }
  fclose(ES);

  //	imprime ajuste
  FILE *BI;
  BI=fopen(NomArch,"a");
  fprintf(BI,"%f %f ",carga[nbin],Ecarga[nbin]);
//  printf("antes de interpolar\n");  fflush(stdout);

///////////////////////////////////////////////////////////////////////////////
////////////////////////////	procedimiento para calcular el bias lineal
///////////////////////////////////////////////////////////////////////////////

//	usando un ajuste lineal en k^2
  Fitab linea1(K2,B,S);
  pendiente=0.;
  corte=linea1.a,segunda=linea1.b,scorte=linea1.siga;
  printf("a=%e sa=%e c=%e\n",corte,scorte,segunda);  fflush(stdout);
  double       coeff[3]={corte,0.*pendiente,pendiente};
  double sigma_coeff[3]={0.5,0.*0.5,1.5};
//  mcmc(i, K, B, S, 2, coeff, sigma_coeff);
//  printf("a=%e sa=%e c=%e\n",corte,scorte,segunda);  fflush(stdout);


  fprintf(BI,"%e %e ",corte,scorte);
//  printf("primer bias\n");  fflush(stdout);
  fprintf(BI,"\n");
  fclose(BI);



  //////////////////////		imprimo ajuste
  sprintf(nombre,"ajuste_%02d_%s",nbin,NomArch);
  ES=fopen(nombre,"w+");
  double err1,err2,err3,err4,err5,err6;
  double val=0.001;		//	imprime el ajuste lineal
  while(val<k_lineal*1.2){
    err1=corte+pendiente*val+segunda*val*val;
    fprintf(ES,"%e %e\n",val,err1);
    val*=1.15;
  }
  fclose(ES);

  printf("calcula el bias...echo\n");
  fflush(stdout);
}















int main(int argc, char **argv){

////////////////////////////            FFTW en paralelo
int fftw_init_threads (void);
printf("0 debe ser direfrente de %d\n",fftw_init_threads ());
void fftw_plan_with_nthreads (int nthreads);
fftw_plan_with_nthreads (nt);
// -lfftw3_threads -lfftw3 -lm -lpthread
// -lfftw3_omp -lfftw3 -lm -fopenmp
omp_set_num_threads(nt);


  allocar();
  char ArchEn[500];
  char ArchSa[500];
  int cuantos;

  if((espectro==1)&&(solo_analiza_bias_voids==0)){
    int lado_caja=(int)Lbox;
    sprintf(ArchEn,"%s%s",prefijo,caso);
    sprintf(ArchSa,"P%s",caso);

    if(distorcion==0){
      densidad(ArchEn,0.,6);
      fftw_execute(dft);
      densidad(ArchEn,0.5,6);
      fftw_execute(dftc);
      media(ArchSa);
    }

    else{
      allocard();
      for(int eje=0;eje<3;eje++){                     //      ciclo sobre los ejes
        cout<<"eje "<<ejes[eje]<<endl;
        densidad_distorcion(ArchEn,eje);
        fftw_execute(dft);
        media(ArchSa);
      }
    }
  }





  FILE * BIA;
  if((bias_halos==1)||(bias_voids==1))
    allocarb();

  if(bias_halos==1){
    sprintf(ArchEn,"halos%s",caso);
    max_min_c(ArchEn);

    if (puntual==1)
      sprintf(ArchEn,"halos%s",caso);			//	puntual
    else
      sprintf(ArchEn,"particulas_halos%s",caso);	//	continuo
    sprintf(ArchSa,"bh%s",caso);

    BIA=fopen(ArchSa,"w+");
    fclose(BIA);

    for(int b=0;b<binb;b++){
      if(puntual==1)
        densidad_bias(ArchEn,b,0.,5);			//	puntual
      else
        densidad_bias_halos(ArchEn,b,0.0);		//	contínuo
      fftw_execute(dftb);
      if(puntual==1)
        densidad_bias(ArchEn,b,0.5,5);			//	puntual
      else
        densidad_bias_halos(ArchEn,b,0.5);		//	contínuo
      fftw_execute(dftc);
      media_bias(ArchSa,b);
    }
  }




  sprintf(ArchEn,"esferas_%s%s",prefijo,caso);
//  sprintf(ArchEn,"voids%s",caso);
  if((bias_voids==1)||(espectro_voids==1))
    max_min_cv(ArchEn);		//	también alloca


  if(bias_voids==1){
    if((puntual==1)&&(solo_analiza_bias_voids==0))
      sprintf(ArchEn,"esferas_%s%s",prefijo,caso);
//      sprintf(ArchEn,"voids%s",caso);			//	puntual
    else if (solo_analiza_bias_voids==0)
      sprintf(ArchEn,"particulas_voids%s",caso);	//	continuo
    sprintf(ArchSa,"bv%s",caso);

    BIA=fopen(ArchSa,"w+");
    fclose(BIA);


    if (solo_analiza_bias_voids==0){
    for(int b=0;b<binb;b++){
      if(puntual==1)
        densidad_bias(ArchEn,b,0.,4);			//	puntual
      else
        densidad_bias_voids(ArchEn,b,0.);		//	continuo
      fftw_execute(dftb);

      if(puntual==1)
        densidad_bias(ArchEn,b,0.5,4);			//	puntual
      else
        densidad_bias_voids(ArchEn,b,0.5);		//	continuo
      fftw_execute(dftc);
      media_bias(ArchSa,b);
    }
    }
    else{
    for(int b=0;b<binb;b++)
      calcula_bias(ArchSa,b);
    }
  }






  
  if(espectro_halos==1){
    if(puntual==1)
      sprintf(ArchEn,"halos%s",caso);			//	puntual
    else
      sprintf(ArchEn,"particulas_halos%s",caso);	//	continuo
    sprintf(ArchSa,"Ph%s",caso);

    if(puntual==1)
      densidad(ArchEn,0.,5);				//	puntual
    else
      densidad_halos(ArchEn,0.);			//	contínuo
    fftw_execute(dft);
    if(puntual==1)
      densidad(ArchEn,0.5,5);				//	puntual
    else
      densidad_halos(ArchEn,0.5);			//	contínuo
    fftw_execute(dftc);
    media(ArchSa);
  }





  if(espectro_voids==1){
    if(puntual==1)
      sprintf(ArchEn,"voids%s",caso);			//	puntual
    else
      sprintf(ArchEn,"particulas_voids%s",caso);	//	continuo
    sprintf(ArchSa,"Pv%s",caso);

    if(puntual==1)
      densidad(ArchEn,0.,4);				//	puntual
    else
      densidad_voids(ArchEn,0.);			//	continuo
    fftw_execute(dft);

    if(puntual==1)
      densidad(ArchEn,0.5,4);				//	puntual
    else
      densidad_voids(ArchEn,0.5);			//	continuo
    fftw_execute(dftc);
    media(ArchSa);
  }


printf("radios de los voids\n");
for(int i=0;i<binb;i++)
printf("%e\n",carga[i]);



  
return 0;

}

