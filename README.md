# micstoich
R package for microbial reaction stoichiometry calculations following Rittmann and McCarty's "Environmental Biotechnology"

# Maintainer
Sasha D. Hafner (https://au.dk/sasha.hafner@bce.au.dk)

# Description
The micstoich package calculates stoichiometry of microbial reactions using the half-reaction approach of Rittmann and McCarty (2001, 2020).
It supports aerobic respiration, anaerobic respiration (nitrate, sulfate, CO2, and other acceptors), and fermentation.
The user specifies the electron donor, acceptor, and optionally a synthesis fraction `fs` to include biomass production.

# Installation

Presently the package is only available on GitHub, so can be installed using the remotes (or devtools) package.

```r
remotes::install_github("AU-BCE-EE/micstoich")
```

# Quick example
```r
library(micstoich)

# Aerobic oxidation of glucose with biomass synthesis
micstoich(donor = "C6H12O6", acceptor = "O2", fs = 0.2)

#     C6H12O6          O2        NH4+       HCO3-     C5H7O2N         CO2 
# -0.04166667 -0.20000000 -0.01000000 -0.01000000  0.01000000  0.21000000 
#         H2O 
#  0.24000000 


# Methanogenesis from acetic acid
micstoich(donor = "CH3COOH", acceptor = "CO2")
# CH3COOH     CO2     CH4 
#  -0.125   0.125   0.125 

# Denitrification
micstoich(donor = "CH3COOH", acceptor = "NO3-", product = "N2")
# CH3COOH    NO3-      H+      N2     CO2     H2O 
#  -0.125  -0.200  -0.200   0.100   0.250   0.350 

# Empirical substrate, mixed fermentation products, with biomass synthesis 
# and specified biomass formula (default is C5H7O2N)
micstoich(
  donor = "C226.36 H404.24 O187.69 N1", 
  product = "(C3H6O3)3 CH3COOH (H2)5", 
  fs = 0.18, 
  bioform = "C3.82H6.45O1.62N0.95"
)

# C226.36 H404.24 O187.69 N1                       NH4+ 
#               -0.001073768               -0.009859736 
#                      HCO3-                        H2O 
#               -0.009859736               -0.038400033 
#    (C3H6O3)3 CH3COOH (H2)5       C3.82H6.45O1.62N0.95 
#                0.015185185                0.011508951 
#                        CO2 
#                0.041916595 

```

Other functions include: `molmass()`, `calcCOD()`, `rxnbal()`.

For more, see the vignette: `vignette("micstoich-start")`.

# References
Rittmann, B.E. and McCarty, P.L. (2001) *Environmental Biotechnology: Principles and Applications*. McGraw-Hill, New York.

Rittmann, B.E. and McCarty, P.L. (2020) *Environmental Biotechnology: Principles and Applications*. 2nd ed. McGraw-Hill, New York.
