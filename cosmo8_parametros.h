// configuración para el isis lcdm
//const double Om=0.267;
//const double Ol=0.733;
//const double H0=71.9;


  const int first_theory=0;
  const int  last_theory=6;


//	from planck_2018.ini
double ombh2          = 0.0223828;
double omch2          = 0.1201075;
double omnuh2         = 0.6451439E-03;
double omk            = 0.;
double hubble         = 67.32117;
double temp_cmb           = 2.7255;
double helium_fraction    = 0.2454006;
double scalar_amp                = 2.100549e-9;
double scalar_spectral_index     = 0.9660499;

const double h2=pow(hubble*0.01,2);
const double H0=hubble;
const double Odm=omch2/h2;
const double Ob=ombh2/h2;
const double Ok=omk;
const double Onu=omnuh2/h2;
const double Yhe=helium_fraction;
const double As=scalar_amp;
const double ns=scalar_spectral_index;
const double T=temp_cmb;
const double Om=Ob+Odm+Onu;
const double Ol=1.-Om;

const double Lbox = 256.;               //	longitud de la caja en [Mpc/h]
const int    np   = 512;                //	numero de particulas en cada dimension

char caso[]="_03.txt";                  //	ultimo snapshot z=0
char prefijo[]="512x256";               //	prefijo para el nombre de los archivos


const double radio_min_todos=1.155162;  //	radio minimo para los bines stack de todos los catalogos
const double radio_max_todos=18.28353;  //	radio maximo para los bines stack de todos los catalogos

const double masa_min_todos=3.8;
const double masa_max_todos=13095.;

double Rmin=radio_min_todos*1.01, Rmax=radio_max_todos*0.99;





//cruzado
//char NameTheories[]={'D','B','A','l','6','5','4'};
//normal
char NameTheories[]={'4','5','6','l','A','B','D'};
  const double parametro[]={1.000000e-04,1.000000e-05,1.000000e-06,1.000000e-16,1.000000e-00,2.000000e-00,3.000000e-00};





  #define Np 1		//	fR0 o zSSB
  int seed;// 3255
  const double dcc=1.686;
  const double pi=4.*atan(1.);
  const int theories=7;
  int MG;
  double mean[Np],stdev[Np]; int contador;	//	analisis
  double ChiMax, ChiMin, param_cosmo[Np],delta[Np],sigma_param_cosmo[Np]; //	variables basicas
  double pMax[Np], pMin[Np], ChiOld, ChiNew;		//	variables auxiliares
  double cn, co;

  int si_bias;
  int si_densidad;
  int si_abundancia;
  int si_imprime_mejor_ajuste_abundancia=0;









  //	abundancia
  #define Npca 6

  const double si_segunda_potencia=1.;

  //	bines para la abundancia
  const int       bin_abundancia=15;
  const int first_bin_abundancia=0;
  const int  last_bin_abundancia=12;
  //	para imprimir el mejor ajuste
  const int bin_print_abundance=200;

  //	coeficientes que vienen de mcmc_ajusta_abundancia4.c
  double param_abundancia[Npca];
  //	sigma evaluado en los R que hay medida de abundancia
  double  sigma_abundancia[theories][bin_abundancia];
  double dsigma_abundancia[theories][bin_abundancia];
  //	prediccion para la abundancia
  double teoria_abundancia[theories][bin_abundancia];

  //	medida abundancia
  double radio_abundancia[theories][bin_abundancia];
  double       abundancia[theories][bin_abundancia];
  double error_abundancia[theories][bin_abundancia];

  //	usados por rodrigo y por mi
  Int nk=694, Nk=nk*10;		//	numero de datos en el espectro camb 99 y 100
  VecDoub P(nk), K(nk), PP(Nk), KK(Nk);
  Doub kmin=1.e8,kmax=-1.;
  double gama;










  //	densidad
  #define uniformiza_densidad 1
  #define densidad_integrada 0

  double A,B;

  //	bines stack a analizar
  const int   bin_stack=10;
  const int first_stack=2;
  const int  last_stack=8;

  //	bines radiales para el perfil de densidad
  const int   bin_radio_densidad=200;			//	maximo número de bines
  const int bines_radio_densidad[bin_stack]={200,200,200,200,150,150,150,150,100,100};
  double radio_densidad[bin_radio_densidad];	//	va entre 0 y 5,7.5,10,		es la distancia al centro del void normalizada por el radio del void

  //	medida perfiles
  double       densidad[theories][bin_stack][bin_radio_densidad];	//	medida	
  double error_densidad[theories][bin_stack][bin_radio_densidad];	//	medida

  double radio_stack[theories][bin_stack];	//	es el radio del void
  double dc=0.867554;
  double xi[theories][124];		//	funcion de correlacion, todas las teorias
  VecDoub radio_xi(124),xi_actual(124);	//	funcion de correlacion actual para interpolar

  int imprime_correlation_fuction;
  int calcula_correlation_fuction;








  //	bias
  double b0,b1,b2;
  double  sigma_bias[theories][bin_stack];	//	rodrigo
  double  radio_bias[theories][bin_stack];	//	medida
  double medida_bias[theories][bin_stack];	//	medida
  double  error_bias[theories][bin_stack];	//	medida


  #define Npd 7	//	densidad: parametros geometricos, c,a,b,G; d1,d2,l1,l2		Number of Parameters for Density
  #define Npbb 3		//	bias: coeficientes, potencias, cortes, constantes
  #define Npbd 10	//	juntando los dos anteriores
  //	parametros libres del modelo conjunto bias perfil
  double param_bd[Npbd];








