
#define PI 3.141592
#define Nphi 500
#define Ndel 200

//	basado en el codigo de Rodrigo
//	devuelve valores de sigma tanto para los perfiles (un valor único para todas las teorias por cada parámetro)

//	basado en FuncionRadioPredice
//	calcula y devuelve sigma y d ln sdigma/d ln R  para comparar con la abundancia medida



/*Structure with the information about the cosmolgy*/
typedef struct Cosmology {
	double H0;
	double w;
	double Ob;
	double Odm;
	double Om;
	double Ol;
	double Ok;
	double Onu;
	double As;
	double ns;
	double Yhe;
	double T;
} COSMOLOGY;

/*Structure with the information about Hu-Sawicki theory parameters*/
typedef struct HuSawicki {
	double n;
	double fr0;
} HS;

/*Structure with the information about Symmetron theory parameters*/
typedef struct Symmetron {
	double LL;
	double zssb;
	double beta;
} SYM;

/*Global variables contaning the information about the cosmoloy and modified gravity*/
double *u;
COSMOLOGY cosmo;
HS hs;
SYM sym;

/*The Hubble parameters over H0 (E(a))*/
double E(double a){
	double resp;
	resp = sqrt(cosmo.Om*pow(a,-3) + cosmo.Ok*pow(a,-2) + cosmo.Ol*pow(a,-3.0*(1.0 + cosmo.w)));
	return resp;
}

/*Derivative of equation above*/
double dE(double a){
	double resp;
	resp = -0.5*(3.0*cosmo.Om*pow(a,-4) + 2.0*cosmo.Ok*pow(a,-3))/E(a);
	return resp;
}

/*Equation for \dot{v} = \ddot{u}*/
double ddu(double a, double u, double v, double mu2, double assb){
	double resp;
	resp = -(4.0/a + dE(a)/E(a))*v - mu2/(a*a*E(a)*E(a))*((pow(assb/a,3) - 1.0)*u + pow(u,3));
	return resp;
}

/*Equation for \dot{u}*/
double du(double a, double u, double v){
	return v;
}

/*Solve the equation for the symmetron field*/
void phi(double aini, double assb, double mu2, double u[], int Ns){
	double *v, h, a, hh;
	double q1, q2, q3, q4, k1, k2, k3, k4;
	int i, j;
	
	/*Allocate the velocity vector*/
	v = (double*)malloc(Ns*sizeof(double));	

	/*Initial conditions*/
	u[0] = 1e-3;	
	v[0] = 1e-3;

	/*Solve the EDO*/
	h = (1.0 - aini)/(Ns-1);
        hh = h*0.5;

	/*Solve the EDO for each k using Runge-Kutta 4*/
	a = aini;
	for(i=0;i<Ns-1;i++){
		q1 = du(a,    u[i],       v[i]);	k1 = ddu(a,    u[i],       v[i],       mu2, assb);
		q2 = du(a+hh, u[i]+hh*q1, v[i]+hh*k1);	k2 = ddu(a+hh, u[i]+hh*q1, v[i]+hh*k1, mu2, assb);
		q3 = du(a+hh, u[i]+hh*q2, v[i]+hh*k2);	k3 = ddu(a+hh, u[i]+hh*q2, v[i]+hh*k2, mu2, assb);
		q4 = du(a+h,  u[i]+h*q3,  v[i]+h*k3);	k4 = ddu(a+h,  u[i]+h*q3,  v[i]+h*k3,  mu2, assb);

		u[i+1] = u[i] + h/6.0*(q1 + 2.0*q2 + 2.0*q3 + q4);	
		v[i+1] = v[i] + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		a = a + h;
	}
}

/*Mass of the symmetron field*/
double Mass_sym(double a, double assb, double mu2){
	double m2;	
	if(a < assb)
		m2 = mu2*(pow(assb/a, 3.0) - 1.0);
	else
		m2 = 2.0*mu2*(1.0 - pow(assb/a, 3.0));
	return m2;
}

/*Wave lenght of the symmetron field*/
double Lamb_sym(double a, double assb){
	double lamb2;	
	if(a < assb)
		lamb2 = 2.0*pow(sym.LL,2)/(pow(assb/a, 3.0) - 1.0);
	else
		lamb2 = pow(sym.LL,2)/(1.0 - pow(assb/a, 3.0));
	return lamb2;
}

/*The value of the potential minimum*/
double min(double a, double assb){
	double resp;
	if(a<assb)
		resp = 0.0;
	else
		resp = sqrt(1.0 - pow(assb/a, 3.0));	
	return resp;
}

/*Function that sumarize the modied Poisson's equation in f(R)*/
double Mu(double k, double a){
	double resp;
	if(MG == 1){
		double m0, m;
		m0 = (100.0/299792.0)*sqrt((cosmo.Om + 4.0*cosmo.Ol)/((hs.n + 1.0)*hs.fr0));
		m = m0*pow((cosmo.Om*pow(a,-3)+4.0*cosmo.Ol)/(cosmo.Om + 4.0*cosmo.Ol),(hs.n+2.0)/2.0);
		resp = ((4.0/3.0)*pow(k,2) + pow(m*a,2))/(pow(k,2) + pow(m*a,2));
	}
	if(MG == 2){
		double lamb2, assb;
	//	int i;
		assb = 1.0/(1.0 + sym.zssb);
		lamb2 = Lamb_sym(a,assb);
	//	i = (int)floor((a - 0.01)/(1.0 - 0.01)*(Nphi*Ndel - 1.0));
		/*resp = 1.0 + 2.0*pow(min(a,assb)*sym.beta,2)/(1.0 + pow(a/k,2)/lamb2);*/
		if(a<=assb)	resp = 1.0;
		else resp = 1.0 + (2.0*pow(sym.beta*sym.LL*k, 2.0)*(1.0 - pow(assb/a,3.0)))/(pow(k*sym.LL,2.0) + a*a - pow(assb,3)/a);
	}
	return resp;
}

/*Equation for \dot{v} = \ddot{\delta}*/
double dv(double k, double a, double d, double v){
	double resp;
	resp = -(3.0/a + dE(a)/E(a))*v + 3.0/2.0*cosmo.Om*Mu(k, a)/(pow(E(a),2)*pow(a,5))*d;
	return resp;
}

/*Equation for \dot{\delta}*/
double dd(double a, double d, double v){
	return v;
}




/*Define the window function for the sigma's calculation*/
double W(double k, double R){
	double resp;
	resp = 3.0/(pow(k*R,2))*(sin(k*R)/(k*R) - cos(k*R));
	return resp;
}

/*Evaluate the square root of matter variance*/
double calc_sigma(int nk,double rad){
	int i;
	double resp;
	resp = 0.0;
	for(i=0;i<nk-2;i++)
		resp += (K[i+1] - K[i])*(P[i]*pow(K[i],2)*pow(W(K[i],rad),2) + P[i+1]*pow(K[i+1],2)*pow(W(K[i+1],rad),2));
	return resp*0.25/(PI*PI);
}






/*Define the window function for the xi's calculation*/
long double W_xi(double k, double R){
	double resp;
        double kr=k*R;
	resp = sin(kr)/(kr);
        if(kr<1e-8)
          return 1.;
        else
  	  return resp;
}

/*Compute the linear matter correlation fuction, xi*/
double calc_xi(double rad){
	int i;
	double resp;
	resp = 0.0;
	for(i=0;i<Nk-2;i++)
		resp += ( KK[i+1] - KK[i] ) * ( PP[i] * pow(KK[i],2) * W_xi( KK[i] , rad ) + PP[i+1] * pow( KK[i+1] , 2 ) * W_xi( KK[i+1] ,rad));
	return resp*0.25/(PI*PI);
}





/*Evaluate the square root of matter variance*/
inline double W_dW(double x, double &w, double &dw){//	funcion ventana al cuadrado
  double xc=x*cos(x);
  double s=sin(x);
  double x2=pow(x,-2);
  double x4=3.*x2*x2;
  w=x*x4*(s-xc);
  dw=x4*(3.*xc+(x*x-3.)*s);}




double aini, aaux, h, a;


void carga_espectro_99(void)
{
printf("carga_espectro_99...");fflush(stdout);
FILE *camb;
double n, fr0, zini;
int i, j;

//	parametros cosmologicos planck/isis
Rmin=0.5;Rmax=280.;
cosmo.H0=H0;//71.9;
cosmo.w=-1.;
cosmo.Odm=Odm;//0.222;
cosmo.Ob=Ob;//0.045;
cosmo.Ol=Ol;//0.733;
cosmo.Ok=Ok;
cosmo.Onu=Onu;//0.;
cosmo.Yhe=Yhe;
cosmo.As=As;//1.085e-09;
cosmo.ns=ns;//1.;
cosmo.T=T;
cosmo.Om = cosmo.Ob + cosmo.Odm;
zini = 99.0;





/*Read the CAMB output (z=99)*/
//camb = fopen("/home/cosmousp/Descargas/CAMB/planck_folder/planck_2018_matterpower_99.dat", "r");
camb = fopen("/home/cosmousp/Descargas/CAMB/isis_folder/isis_matterpower_99.dat", "r");
if (camb == NULL) {
	printf("Unable to open %s\n", "test_matterpower99.dat");
	exit(0);
}

nk = 0;
while(fscanf(camb,"%lf %lf", &K[nk], &Pini[nk])!=EOF){
	nk ++;
}
fclose(camb);


/*Read the aux power spectrum (z=100)*/
//camb = fopen("/home/cosmousp/Descargas/CAMB/planck_folder/planck_2018_matterpower_100.dat", "r");
camb = fopen("/home/cosmousp/Descargas/CAMB/isis_folder/isis_matterpower_100.dat", "r");
if (camb == NULL) {
	printf("Unable to open %s\n", "test_matterpower100.dat");
	exit(0);
}

nk = 0;
while(fscanf(camb,"%lf %lf", &K[nk], &Paux[nk])!=EOF){
        if(kmin>K[nk])	kmin=K[nk];
	if(kmax<K[nk])	kmax=K[nk];
	nk ++;
}
fclose(camb);

/*Construct the initial condictions*/
aini = 1.0/(1.0 + zini);
aaux = 1.0/(2.0 + zini);

printf("echo\n");fflush(stdout);
}









void main_rodrigo(int id_theory){
int i,j;
double q1, q2, q3, q4, k1, k2, k3, k4, trash;


if(MG == 1){
  hs.n=1.;
  hs.fr0=pow(10.,param_cosmo[0]);
//        printf("%d %f %le\n", MG, hs.n, hs.fr0);
}
if(MG == 2){
  sym.zssb=param_cosmo[0];
  sym.beta=1.;
  sym.LL=1.;
//        printf("%d %le %f %f\n", MG, sym.zssb, sym.beta, sym.LL);
}

for(i=0;i<nk;i++){
	d[i] = sqrt(Pini[i]);
	v[i] = (d[i] - sqrt(Paux[i]))/(aini - aaux);
}

/*Solve the EDO*/
int Ns = Ndel;	/*number of steps*/
h = (1.0 - aini)/(Ns-1);
double hh=0.5*h;


for(j=0;j<nk;j++){
/*Solve the EDO for each k using Runge-Kutta 4*/
	a = aini;
	for(i=0;i<Ns;i++){ 
		q1 = dd(a,    d[j],       v[j]);	k1 = dv(K[j], a,    d[j],       v[j]);
		q2 = dd(a+hh, d[j]+hh*q1, v[j]+hh*k1);	k2 = dv(K[j], a+hh, d[j]+hh*q1, v[j]+hh*k1);
		q3 = dd(a+hh, d[j]+hh*q2, v[j]+hh*k2);	k3 = dv(K[j], a+hh, d[j]+hh*q2, v[j]+hh*k2);
		q4 = dd(a+h,  d[j]+ h*q3, v[j]+ h*k3);	k4 = dv(K[j], a+h,  d[j]+ h*q3, v[j]+ h*k3);

		d[j] = d[j] + h/6.0*(q1 + 2.0*q2 + 2.0*q3 + q4);	
		v[j] = v[j] + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		a = a + h;
	}
}

  char nomAr[210];
//  FILE * MPS;
//  sprintf(nomAr,"../%c/P.txt",NameTheories[id_theory]);
//  MPS=fopen(nomAr,"w+");
  for(i=0;i<nk;i++){
	P[i] = d[i]*d[i]/1.0285;
//    fprintf(MPS,"%le %le\n",K[i],P[i]);
  }
//  fclose(MPS);





//	calcula sigma, lo necesito para el biasy para la densidad
if((si_bias==1)||(si_densidad==1)){
  /*Evaluating the square root of variance (\sigma)*/

  for(int s=first_stack;s<=last_stack;s++){
    double rad=radio_bias[id_theory][s];
    sigma_bias[id_theory][s] = sqrt(calc_sigma(nk,rad));
//printf("%le %le %le\n",radio_stack[id_theory][s],radio_bias[id_theory][s],sigma_bias[id_theory][s]);
  }
}


//	es una fucion externa, solo sirve para calcular la funcion de correlacion y sigma de cada teoria e imprimirla
//	es llamada por cosmo7_matter_correlation_function.c
if(imprime_correlation_fuction==1){

  FILE * ES;
  sprintf(nomAr,"../%c/lin_matter_power_spectrum.txt",NameTheories[id_theory]);
  ES=fopen(nomAr,"w+");
  for(i=0;i<nk;i++){
    fprintf(ES,"%le %le\n",K[i],P[i]);
  }
  fclose(ES);

  //spline
  Spline_interp pk(K,P);
  //write
  double k_left=K[0], k_right=K[nk-1], c=pow(k_right/k_left,1./(Nk-1)), k_new=k_left;printf("kleft=%le k_right=%le k_new=%le c=%le\n",k_left,k_right,k_new,c);
  for(i=0;i<Nk;i++){
    KK[i]=k_new;
    PP[i]=pk.interp(k_new);
    k_new*=c;//printf("%le %le\n",KK[i],PP[i]);
  }

  FILE * XI;
  FILE * SI;
  sprintf(nomAr,"../%c/lin_matter_correlation_fuction.txt",NameTheories[id_theory]);
  XI=fopen(nomAr,"w+");
  sprintf(nomAr,"../%c/lin_matter_sigma.txt",NameTheories[id_theory]);
  SI=fopen(nomAr,"w+");
  for(double rad=5.e-3; rad<2.; rad+=1.e-1){
    fprintf(XI,"%le %le\n",rad,calc_xi(rad));
    fprintf(SI,"%le %le\n",rad,sqrt(calc_sigma(nk,rad)));}
  for(double rad=2.; rad<11.; rad+=2.e-1){
    fprintf(XI,"%le %le\n",rad,calc_xi(rad));
    fprintf(SI,"%le %le\n",rad,sqrt(calc_sigma(nk,rad)));}
  for(double rad=11.; rad<20.; rad+=5.e-1){
    fprintf(XI,"%le %le\n",rad,calc_xi(rad));
    fprintf(SI,"%le %le\n",rad,sqrt(calc_sigma(nk,rad)));}
  for(double rad=20.; rad<120.; rad+=1.e-0){
    fprintf(XI,"%le %le\n",rad,calc_xi(rad));
    fprintf(SI,"%le %le\n",rad,sqrt(calc_sigma(nk,rad)));}
  fclose(XI);
  fclose(SI);
}




//	basada en el parrafo anterior
//	calcula y guarda la funcion de correlacion de la materia
//	la necesito para la parte externa de los perfiles de densidad
if(calcula_correlation_fuction==1){
  //spline
  Spline_interp pk(K,P);
  //write
  double k_left=K[0], k_right=K[nk-1], c=pow(k_right/k_left,1./(Nk-1)), k_new=k_left;//printf("kleft=%le k_right=%le k_new=%le c=%le\n",k_left,k_right,k_new,c);
  for(i=0;i<Nk;i++){
    KK[i]=k_new;
    PP[i]=pk.interp(k_new);
    k_new*=c;//printf("%le %le\n",KK[i],PP[i]);
  }
  int id_xi=0;
  for(double rad=5.e-3; rad<2.; rad+=1.e-1){
    radio_xi[id_xi]=rad;
    xi[id_theory][id_xi]=calc_xi(rad);
    id_xi++;}
  for(double rad=2.; rad<11.; rad+=2.e-1){
    radio_xi[id_xi]=rad;
    xi[id_theory][id_xi]=calc_xi(rad);
    id_xi++;}
  for(double rad=11.; rad<20.; rad+=5.e-1){
    radio_xi[id_xi]=rad;
    xi[id_theory][id_xi]=calc_xi(rad);
    id_xi++;}
  for(double rad=20.; rad<120.; rad+=1.e-0){
    radio_xi[id_xi]=rad;
    xi[id_theory][id_xi]=calc_xi(rad);
    id_xi++;}//printf("id_xi=%d\n",id_xi);fflush(stdout);


}


}


























