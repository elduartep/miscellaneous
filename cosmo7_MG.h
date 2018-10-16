
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
	double A;
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
	double *v, h, a;
	double q1, q2, q3, q4, k1, k2, k3, k4;
	int i, j;
	
	/*Allocate the velocity vector*/
	v = (double*)malloc(Ns*sizeof(double));	

	/*Initial conditions*/
	u[0] = 1e-3;	
	v[0] = 1e-3;

	/*Solve the EDO*/
	h = (1.0 - aini)/(Ns-1);

	/*Solve the EDO for each k using Runge-Kutta 4*/
	a = aini;
	for(i=0;i<Ns-1;i++){
		q1 = du(a, u[i], v[i]);
		k1 = ddu(a, u[i], v[i], mu2, assb);
		q2 = du(a+h/2.0, u[i]+h/2.0*q1, v[i]+h/2.0*k1);
		k2 = ddu(a+h/2.0, u[i]+h/2.0*q1, v[i]+h/2.0*k1, mu2, assb);
		q3 = du(a+h/2.0, u[i]+h/2.0*q2, v[i]+h/2.0*k2);
		k3 = ddu(a+h/2.0, u[i]+h/2.0*q2, v[i]+h/2.0*k2, mu2, assb);
		q4 = du(a+h, u[i]+h*q3, v[i]+h*k3);		
		k4 = ddu(a+h, u[i]+h*q3, v[i]+h*k3, mu2, assb);

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
		int i;
		assb = 1.0/(1.0 + sym.zssb);
		lamb2 = Lamb_sym(a,assb);
		i = (int)floor((a - 0.01)/(1.0 - 0.01)*(Nphi*Ndel - 1.0));
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
		resp += (K[i+1] - K[i])/2.0*(P[i]*pow(K[i],2)*pow(W(K[i],rad),2) + P[i+1]*pow(K[i+1],2)*pow(W(K[i+1],rad),2));
	return resp/(2.0*PI*PI);
}





/*Evaluate the square root of matter variance*/
inline double W_dW(double x, double &w, double &dw){//	funcion ventana al cuadrado
  double xc=x*cos(x);
  double s=sin(x);
  double x2=pow(x,-2);
  double x4=3.*x2*x2;
  w=x*x4*(s-xc);
  dw=x4*(3.*xc+(x*x-3.)*s);}




double Pini[647], Paux[647], d[647], v[647], aini, aaux, h, a;


void carga_espectro_99(void)
{
printf("carga_espectro_99...");fflush(stdout);
FILE *camb;
double n, fr0, zini;
int i, j;

//	parametros cosmologicos
Rmin=0.5;Rmax=280.;
cosmo.H0=71.9;
cosmo.w=-1.;
cosmo.Odm=0.222;
cosmo.Ob=0.045;
cosmo.Ol=0.733;
cosmo.Ok=0.;
cosmo.Onu=0.;
cosmo.Yhe=0.24;
cosmo.A=1.085e-09;
cosmo.ns=1.;
cosmo.T=2.726;
cosmo.Om = cosmo.Ob + cosmo.Odm;
zini = 99.0;





/*Read the CAMB output (z=99)*/
camb = fopen("../l/test_matterpower99.dat", "r");
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
camb = fopen("../l/test_matterpower100.dat", "r");
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
//        printf("%d %f %le\n", model, hs.n, hs.fr0);
}
if(MG == 2){
  sym.zssb=param_cosmo[0];
  sym.beta=1.;
  sym.LL=1.;
//        printf("%d %le %f %f\n", model, sym.zssb, sym.beta, sym.LL);
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

for(i=0;i<nk;i++)
	P[i] = d[i]*d[i]/1.0285;




//	calcula sigma, lo necesito para el perfil de densidad
if(si_densidad==1){
  /*Evaluating the square root of variance (\sigma)*/
  for(int s=first_stack;s<=last_stack;s++){
    double rad=radio_stack[id_theory][s];
    sigma_stack[id_theory][s] = sqrt(calc_sigma(nk,rad));
  }
}


}
