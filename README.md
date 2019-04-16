# voids
spherical void finding algorithm and analysis routines


## spherical void finder
VoidFinderPrueba2.c

requires: a tracer catalog

returns: an spherical void catalog






## Desity Profile (DP)

### compute density profile
Perfiles_2019_octantes.c

requires: a void catalog and the original tracers catalog

returns: the spherical averaged density profile for stacked voids


### trend of the density profile parameters
mcmc_perfill_lineal_4.c ..... .c

requires: spherical averaged density profile for stacked voids for different MG cases of the same theory

returns: the best fit parameters (associated a given expression describing the voids density profil) for the different staks and MG cases

### best density profile
cosmo7_desciende_bias_densidad.c


requires: spherical averaged density profile for stacked voids for different MG cases of the same theory

returns: the best fit coefficients (given an expression for the DP parameters dependece, which should be infered fron the last result) for the different common for all the staks and MG cases






## Abundance (AB)

### compute abundance or AB
cosmo7_FuncionRadioMide.c


requires: a void catalog

returns: the void abundance


### trend of the abundance parameters
mcmc_ajusta_abundancia4.c
cosmo7_mcmc_abundancia.c

requires: void abundance for different MG cases of the same theory

returns: the best fit parameters (associated a given expression describing the voids abundance) for the different MG cases

### best abundance
cosmo7_desciende_ajusta_abundancia.c
cosmo7_desciende_abundacia.c

requires: voids abundance for different MG cases of the same theory

returns: the best fit coefficients (given an expression for the AB parameters dependece, which should be infered fron the last result) for the different common for all the MG cases




## Linear Void Bias (LB)


### Measure
Bias_pragma_lb.c

requires: matter and void catalogs

return: linear void-matter bias

### Trend
cosmo8_mcmc_bias_densidad.c

### best parameters
cosmo7_desciende_bias.c






## applay AB and DP learned expressions

### infer cosmology
cosmo8.c

requires: coefficients of the abundace and density profile fits, AB and DP measurements for the different MG cases

return: given the coefficients that describe the parameters of the AB and DP fits, returns the value of the MG parameter for each MG case









# Void and halo Profiles
Perfiles_2018_octantes.c

computes halo density profiles and voids density profiles


# Power Spectrum and bias
Bias_pragma_lb.c

computes matter power spectrum, halo bias and voids bias


# Spherical Halo Finder
HaloFinder.c


# Correlation Function
correl_pix_nuevo.cpp  measure

cosmo7_matter_correlation_function.c  linear-theory




# MG P(k)
MG_Power_Spectrum.c ./a.out params_mg.dat





# planck cosmology analysis

UiO

Bias_pragma_ln.c measures the linear bias using fourier transform

funcionRadio.c counts the voids

Perfiles_2019.c measures the void-matter correlation fuction

DELL

cosmo8_covariance.c computes the covariance matrix for the bias, abundance and correlation

invert_ma.py inverts the covariance matrix




# Halos

HaloFinder2019.c using spherical overdensity

SplitHalosByConcentration.c split so200 halo catalog into low and high concentration estimated by the ratio r200/r400

SplitHalosByConcentrationRockstar.c split rockstar halo catalog into low and high concentration estimated by the ratio (Mvir/M500c)**(1/3)

Perfiles_2019_median.c estimates halo profiles with the median over 48 angular patchs given by healpix

cosmo8_halos_subsample.c subsample rockstar haloes by selecting those which are not inside others

savitzki3.ipynb savitzki4.ipynb : locates the splash-back radius of halos


