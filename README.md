# voids
spherical void finding algorithm and analysis routines


## spherical void finder
VoidFinderPrueba3.c

requires: a tracer catalog

returns: an spherical void catalog






## desity profile DP

### compute density profile
Perfiles_2018_octantes.c

requires: a void catalog and the original tracers catalog

returns: the spherical averaged density profile for stacked voids


### trend of the density profile parameters
mcmc_perfill_lineal_4.c ..... .c

requires: spherical averaged density profile for stacked voids for different MG cases of the same theory

returns: the best fit parameters (associated a given expression describing the voids density profil) for the different staks and MG cases

### best density profile
cosmo7_desciende_perfil.c


requires: spherical averaged density profile for stacked voids for different MG cases of the same theory

returns: the best fit coefficients (given an expression for the DP parameters dependece, which should be infered fron the last result) for the different common for all the staks and MG cases






## abundance AB

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




## Bias

### Trend
cosmo7_mcmc_bias.c

### best parameters
cosmo7_desciende_bias.c






## applay AB and DP learned expressions

### infer cosmology
cosmo7.c

requires: coefficients of the abundace and density profile fits, AB and DP measurements for the different MG cases

return: given the coefficients that describe the parameters of the AB and DP fits, returns the value of the MG parameter for each MG case









# Void and halo Profiles
Perfiles_2018_octantes_velocidad.c

computes halo density profiles and voids density profiles


# Power Spectrum and bias
Bias_pragma_lb.c

computes matter power spectrum, halo bias and voids bias


# Spherical Halo Finder
HaloFinder.c


# Correlation Function
correl_pix_nuevo.cpp
