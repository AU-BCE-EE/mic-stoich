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

# Methanogenesis from acetic acid
micstoich(donor = "CH3COOH", acceptor = "CO2")

# Denitrification
micstoich(donor = "CH3COOH", acceptor = "NO3-", product = "N2")
```

For more, see the vignette: `vignette("micstoich-start")`.

# References
Rittmann, B.E. and McCarty, P.L. (2001) *Environmental Biotechnology: Principles and Applications*. McGraw-Hill, New York.

Rittmann, B.E. and McCarty, P.L. (2020) *Environmental Biotechnology: Principles and Applications*. 2nd ed. McGraw-Hill, New York.
