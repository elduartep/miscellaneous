//	hace una submuestra de halos aislados, halos rockstar
//	tambien genera subcatalogos con baja y alta concentracion

#include <iostream>
#include <fstream>
#include <cmath>



#include "../l/nr3.h"

#include "../l/parametros.h"

double minmin=0.,maxmax=1.e23;
double Lbox2=128.;

int ID, DescID;
double mvir, Vmax, Vrms, rvir, Rs;
int Npart;
double x,y,z,VX, VY, VZ, JX, JY, JZ, Spin, rs_klypin, Mvir_all, M200b, M200c, m500c, M2500c;
double Xoff, Voff, spin_bullock, b_to_a, c_to_a, Ax, Ay, Az, b_to_a_500c, c_to_a_500c;
double Ax_500c, Ay_500c, Az_500c, TU, M_pe_Behroozi, M_pe_Diemer;

double *X, *Y, *Z , *R, *Rvir, *Mvir, *M500c;

double min_mass_low[7]={3.e18};
double min_mass_high[7]={3.e18};
double max_mass_low[7]={0.};
double max_mass_high[7]={0.};

int NumObjetos;

double c,mean;


//////////////////////////////////////////////////////////////////////////////////
void alloca_xyzr(void){
//////////////////////////////////////////////////////////////////////////////////
  printf("alloca halos...");fflush(stdout);
  if(!(X = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory X.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Y = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Y.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Z = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Z.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Rvir = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Rvir.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(R = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory R.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(Mvir = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory Mvir.\n");
      fflush(stderr);
      exit(0);
    }
  if(!(M500c = (double*) malloc ( NumObjetos * sizeof(double))))
    {
      fprintf(stderr, "failed to allocate memory M500c.\n");
      fflush(stderr);
      exit(0);
    }
  printf("hecho\n");fflush(stdout);
}


//////////////////////////////////////////////////////////////////////////////////
void read_halos(int id_theory){
//////////////////////////////////////////////////////////////////////////////////
  printf("lee halos...\n");fflush(stdout);
  char NomArch[300];
  FILE * HA;
  sprintf(NomArch,"../%c/out_0.list",NameTheories[id_theory]);
  HA=fopen(NomArch,"r");
  char * line = NULL;
  size_t len = 0;
  int i;
  for(i=0;i<16;i++)
    getline(&line, &len, HA);
  NumObjetos=0;
  while(getline(&line, &len, HA)!=-1) NumObjetos++;
  printf("NumObjetos=%d\n",NumObjetos);fflush(stdout);
  rewind(HA);
  alloca_xyzr();
  double min=1.e18,max=0.;
  for(i=0;i<16;i++)
    getline(&line, &len, HA);
  for(i=0;i<NumObjetos;i++){
    fscanf(HA,"%d %d %le %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %le %le %le %lf %lf %le %le %le %le %le %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %le %le\n",&ID, &DescID, &mvir, &Vmax, &Vrms, &rvir, &Rs, &Npart, &x, &y, &z, &VX, &VY, &VZ, &JX, &JY, &JZ, &Spin, &rs_klypin, &Mvir_all, &M200b, &M200c, &m500c, &M2500c, &Xoff, &Voff, &spin_bullock, &b_to_a, &c_to_a, &Ax, &Ay, &Az, &b_to_a_500c, &c_to_a_500c, &Ax_500c, &Ay_500c, &Az_500c, &TU, &M_pe_Behroozi, &M_pe_Diemer);
    X[i]=x;
    Y[i]=y;
    Z[i]=z;
    Rvir[i]=rvir;
    R[i]=rvir;
    Mvir[i]=mvir;
    M500c[i]=m500c;
    if(min>mvir) min=mvir;
    if(max<mvir) max=mvir;
  }
  fclose(HA);

  printf("max=%le min=%le\n",max,min);
  if(minmin<min) minmin=min;
  if(maxmax>max) maxmax=max;
  printf("lee halos...hecho\n");fflush(stdout);
}






//////////////////////////////////////////////////////////////////////////////////
void free_xyzr(void){
//////////////////////////////////////////////////////////////////////////////////
free(Mvir);
free(M500c);
free(R);
free(Rvir);
free(Z);
free(Y);
free(X);

}


//////////////////////////////////////////////////////////////////////////////////
void subsample_halos(int id_theory){
//////////////////////////////////////////////////////////////////////////////////
  printf("subsample halos...");fflush(stdout);
  int t=id_theory;
  int i,j,k;
  double dx,dy,dz,rvjr,rad;
  for(i=0;i<NumObjetos;i++){
    x=X[i];
    y=Y[i];
    z=Z[i];
    rvir=Rvir[i];
    for(j=i+1;j<NumObjetos;j++){
      rvjr=Rvir[j];
      if(rvir<rvjr){
        k=i;
        rad=rvjr;}
      else{
        k=j;
        rad=rvir;}
      double dist2=0.;
      dx=abs(x-X[j]);	if(dx>Lbox2)	dx-=Lbox2;
      dy=abs(y-Y[j]);	if(dy>Lbox2)	dy-=Lbox2;
      dz=abs(z-Z[j]);	if(dz>Lbox2)	dz-=Lbox2;
      dist2=dx*dx+dy*dy+dz*dz;
      if(sqrt(dist2)<(rad*0.001))
        R[k]=-1.;	//	declaro a j no aislado
    }
  }
 
  //	computes the average of concentration
  int NumHalos=0;
  mean=0.;
  for(i=0;i<NumObjetos;i++){
    if((R[i]>0.)&&(M500c[i]>0.)){		//	isolated
      c=pow(Mvir[i]/M500c[i],1./3);
      mean+=c;
      NumHalos++;
    }
  }

  mean/=NumHalos;
  printf("\nmean=%le\n",mean);
  mean=1.23;
  printf("\nmean=%le\n",mean);

  char NomArch[300];
  FILE * HA;
  sprintf(NomArch,"../%c/out_isolated.list",NameTheories[id_theory]);
  HA=fopen(NomArch,"w+");
  FILE * HB;
  sprintf(NomArch,"../%c/out_all.list",NameTheories[id_theory]);
  HB=fopen(NomArch,"w+");



  char catalogo_halos_low[500];
  char catalogo_halos_high[500];
  sprintf(catalogo_halos_low,"../%c/halos_lowc.list",NameTheories[t]);
  sprintf(catalogo_halos_high,"../%c/halos_highc.list",NameTheories[t]);
  FILE *LO;
  FILE *HI;
  LO=fopen(catalogo_halos_low,"w+");
  HI=fopen(catalogo_halos_high,"w+");

  //	split the catalog
  for(i=0;i<NumObjetos;i++){
    if(R[i]>0.){		//	isolated
      fprintf(HA,"%f %f %f %f %f\n",X[i],Y[i],Z[i],Mvir[i]*1.e-11,Rvir[i]*0.001);
      if(M500c[i]>0.){
        c=pow(Mvir[i]/M500c[i],1./3);
        if(c<mean){	//	high concentration
          fprintf(HI,"%lf %lf %lf %lf %lf\n",X[i],Y[i],Z[i],Mvir[i]*1.e-11,Rvir[i]*0.001);
          if(min_mass_high[t]>Mvir[i])	min_mass_high[t]=Mvir[i];
          if(max_mass_high[t]<Mvir[i])	max_mass_high[t]=Mvir[i];
        }
        else{		//	low concentration
          fprintf(LO,"%lf %lf %lf %lf %lf\n",X[i],Y[i],Z[i],Mvir[i]*1.e-11,Rvir[i]*0.001);
          if(min_mass_low[t]>Mvir[i])	min_mass_low[t]=Mvir[i];
          if(max_mass_low[t]<Mvir[i])	max_mass_low[t]=Mvir[i];
        }
      }
    }
    fprintf(HB,"%f %f %f %f %f\n",X[i],Y[i],Z[i],Mvir[i]*1.e-11,Rvir[i]*0.001);
  }
  fclose(HA);
  fclose(HB);
  fclose(LO);
  fclose(HI);

  free_xyzr();
  printf("subsample halos...hecho\n");fflush(stdout);
}








//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(){
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
  for (int id_theory=0;id_theory<=6;id_theory++){
    printf("\nid_theory=%d\n",id_theory);fflush(stdout);
    read_halos(id_theory);
    subsample_halos(id_theory);
  }
  printf("max=%le min=%le\n",maxmax,minmin);

double min_mass=0.;
double max_mass=3.e18;
for(int t=0;t<7;t++){
  if(min_mass<min_mass_low[t])	min_mass=min_mass_low[t];
  if(min_mass<min_mass_high[t])	min_mass=min_mass_high[t];
  if(max_mass>max_mass_low[t])	max_mass=max_mass_low[t];
  if(max_mass>max_mass_high[t])	max_mass=max_mass_high[t];
}
printf("sup mass = %le\n",min_mass*1.e-11);
printf("inf mass = %le\n",max_mass*1.e-11);

  return 0;
}










