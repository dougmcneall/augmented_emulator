# augmented_emulator

Code developing testing a Gaussian Process emulator, augmented with model outputs as inputs.
This code augments a Gaussian process emulator using temperature and precipitation (usually considered model outputs) as inputs. It allows biases in those fields to be accounted for when choosing plausible values of input parameters for the land surface component of the reduced-resolution climate model, FAMOUS.

This is the main analysis code for:
McNeall, D.J., Williams J., Betts, R.A. , Booth, B.B.B., Challenor, P.G. , Good, P. & Wiltshire A. (2019) 
Submitted to Geoscientific Model Development
Contact Doug McNeall dougmcneall@gmail.com @dougmcneall

Main R code for analysis is augmented_emulator.R
