// configuración para el isis lcdm

const double Om=0.267;
const double Ol=0.733;
const double H0=71.9;

const double Lbox = 256.;               //	longitud de la caja en [Mpc/h]
const int    np   = 512;                //	numero de particulas en cada dimension

char caso[]="_03.txt";                  //	ultimo snapshot z=0
char prefijo[]="512x256";               //	prefijo para el nombre de los archivos


const double radio_min_todos=1.155162;  //	radio minimo para los bines stack de todos los catalogos
const double radio_max_todos=18.28353;  //	radio maximo para los bines stack de todos los catalogos

double Rmin=radio_min_todos*1.01, Rmax=radio_max_todos*0.99;

char NameTheories[]={'4','5','6','l','A','B','D'};
const double parametro[]= {1.000000e-04,1.000000e-05,1.000000e-06,1.000000e-18,1.000000e-00,2.000000e-00,3.000000e-00};




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

  int si_densidad;
  int si_abundancia;









  //	abundancia
  #define Npca 6

  const double si_segunda_potencia=1.;

  //	bines para la abundancia
  const int       bin_abundancia=15;
  const int first_bin_abundancia=0;
  const int  last_bin_abundancia=12;

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
  VecDoub P(647), K(647);
  Int nk;
  Doub kmin=1.e8,kmax=-1.;
  double gama;










  //	densidad
  #define uniformiza_densidad 1
  #define densidad_integrada 0
  #define Npd 3	//	densidad: parametros geometricos, c,a,b,G(r.2;sigma) (más la constante del G(sigma))
  #define Ncd 5		//	densidad: coeficientes, potencias, cortes, constantes

  double alpha,beta,A,B;

  //	bines stack a analizar
  const int   bin_stack=10;
  const int first_stack=1;
  const int  last_stack=7;

  //	bines radiales para el perfil de densidad
  const int   bin_radio_densidad=200;			//	maximo número de bines
  const int bines_radio_densidad[bin_stack]={200,200,200,200,150,150,150,150,100,100};

  //	coeficientes que bienen de mcmc_todo_libre5nuevofast_2.c
  double param_densidad[Npd][Ncd];

  //	medida perfiles
  double radio_densidad[bin_radio_densidad];	//	va entre 0 y 5,7.5,10
  double       densidad[theories][bin_stack][bin_radio_densidad];
  double error_densidad[theories][bin_stack][bin_radio_densidad];

  double radio_stack[theories][bin_stack];	//	es el radio del void
  double sigma_stack[theories][bin_stack];
  double dc=0.867554;






