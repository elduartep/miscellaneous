# voids
spherical void finding algorithm and analyzing routines


## spherical void finder
requires: a tracer catalog
returns: an spherical void catalog






## desity profile DP

### compute density profile
requires: a void catalog and the original tracers catalog
returns: the spherical averaged density profile for stacked voids

trend of the density profile parameters
requires: spherical averaged density profile for stacked voids for different MG cases of the same theory
returns: the best fit parameters (associated a given expression describing the voids density profil) for the different staks and MG cases

best density profile
requires: spherical averaged density profile for stacked voids for different MG cases of the same theory
returns: the best fit coefficients (given an expression for the DP parameters dependece, which should be infered fron the last result) for the different common for all the staks and MG cases






## abundance AB

compute abundance or AB
requires: a void catalog
returns: the void abundance


trend of the abundance parameters
requires: void abundance for different MG cases of the same theory
returns: the best fit parameters (associated a given expression describing the voids abundance) for the different MG cases

best abundance
requires: voids abundance for different MG cases of the same theory
returns: the best fit coefficients (given an expression for the AB parameters dependece, which should be infered fron the last result) for the different common for all the MG cases




## applay AB and DP learned expressions

infer cosmology
requires: coefficients of the abundace and density profile fits, AB and DP measurements for the different MG cases
return: given the coefficients that describe the parameters of the AB and DP fits, returns the value of the MG parameter for each MG case

