---
title: "Getting started with micstoich"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started with micstoich}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Testing


``` r
micstoich('H2', acceptor = 'CO2', product = 'CH3COOH')
```

```
##      H2     CO2 CH3COOH     H2O 
##  -0.500  -0.250   0.125   0.250
```

``` r
micstoich('Fe+2', acceptor = 'O2')
```

```
##  Fe+2    O2    H+  Fe+3   H2O 
## -1.00 -0.25 -1.00  1.00  0.50
```

# Some simple examples

Let's use glucose oxidation as an example.


``` r
micstoich(donor = 'C6H12O6', acceptor = 'O2')
```

```
##     C6H12O6          O2         CO2         H2O 
## -0.04166667 -0.25000000  0.25000000  0.25000000
```

Reactants have negative values and products positive.
The compounds themselves are given in the element names.

There is no biomass production in that example.
We can add it by providing a value for `fs`, perhaps 0.5 for this aerobic reaction.


``` r
micstoich(donor = 'C6H12O6', acceptor = 'O2', fs = 0.5)
```

```
##     C6H12O6          O2        NH4+       HCO3-     C5H7O2N         CO2 
## -0.04166667 -0.12500000 -0.02500000 -0.02500000  0.02500000  0.15000000 
##         H2O 
##  0.22500000
```

The default biomass composition is `C5H7O2N`, but it can be changed with `bioform`.

Is that reaction really balanced?
It *should* be, but we can check it or any reaction with `rxnbal()`.



``` r
r1 <- micstoich(donor = 'C6H12O6', acceptor = 'O2', fs = 0.5)
r1
```

```
##     C6H12O6          O2        NH4+       HCO3-     C5H7O2N         CO2 
## -0.04166667 -0.12500000 -0.02500000 -0.02500000  0.02500000  0.15000000 
##         H2O 
##  0.22500000
```

``` r
rxnbal(r1)
```

It returns `TRUE` invisibly if everything is OK.


``` r
isTRUE(rxnbal(r1))
```

```
## [1] TRUE
```

Here is an example that is not balanced.


``` r
rxnbal(c(CH3COOH = -2, CH4 = 1, CO2 = 1))
```

```
## Warning in rxnbal(c(CH3COOH = -2, CH4 = 1, CO2 = 1)): Elemental balance is off
## in reaction: -2 CH3COOH, 1 CH4, 1 CO2.
```

```
##  C  H  O 
## -2 -4 -2
```

We get a warning plus the overall balance, which shows 2 more moles of C and O and 4 of H on the reactant side than product side.
This is indicated by the negative values.

# Other electron acceptors

Let's use acetic acid as a donor with some other acceptors.

First nitrate for complete denitrification to N2.


``` r
micstoich(donor = 'CH3COOH', acceptor = 'NO3-')
```

```
## Error in `halfrxn[[names(halfrxn)[reactant == rnm & product == pnm]]]`:
## ! attempt to select less than one element in get1index
```

Here we see we need to specify a product as well, because there are multiple possibilities.
We can find the one we want in the error message, where they are listed in `acceptor product` pairs: `'NO3- N2'`.
Or more easily, view them all with the `halfrxns()` function.


``` r
halfrxns()
```

```
##    acceptor product                                         reaction
## 1        O2     H2O                              -0.25O2 -1H+ 0.5H2O
## 2       CO2     CH4                  -0.125CO2 -1H+ 0.125CH4 0.25H2O
## 3      NO3-      N2                     -0.2NO3- -1.2H+ 0.1N2 0.6H2O
## 4      NO3-     NH3            -0.125NO3- -1.25H+ 0.125NH4+ 0.375H2O
## 5      NO3-    NH4+            -0.125NO3- -1.25H+ 0.125NH4+ 0.375H2O
## 6      NO3-    NO2-                     -0.5NO3- -1H+ 0.5NO2- 0.5H2O
## 7      NO2-     NH3    -0.16667NO2- -1.3333H+ 0.16667NH4+ 0.33333H2O
## 8      NO2-    NH4+    -0.16667NO2- -1.3333H+ 0.16667NH4+ 0.33333H2O
## 9      NO2-      N2      -0.33333NO2- -1.3333H+ 0.16667N2 0.66667H2O
## 10       N2     NH3                 -0.16667N2 -1.3333H+ 0.33333NH4+
## 11       N2    NH4+                 -0.16667N2 -1.3333H+ 0.33333NH4+
## 12     Fe+3    Fe+2                                     -1Fe+3 1Fe+2
## 13       H+      H2                                       -1H+ 0.5H2
## 14    SO4-2     H2S -0.125SO4-2 -1.1875H+ 0.0625H2S 0.0625HS- 0.5H2O
## 15    SO4-2   SO3-2                   -0.5SO4-2 -1H+ 0.5SO3-2 0.5H2O
## 16    SO4-2       S      -0.16667SO4-2 -1.3333H+ 0.16667S 0.66667H2O
## 17    SO4-2  S2O3-2          -0.25SO4-2 -1.25H+ 0.125S2O3-2 0.625H2O
```

The one we want is in the third row.


``` r
micstoich(donor = 'CH3COOH', acceptor = 'NO3-', product = 'N2')
```

```
## CH3COOH    NO3-      H+      N2     CO2     H2O 
##  -0.125  -0.200  -0.200   0.100   0.250   0.350
```

How about sulfate?


``` r
micstoich(donor = 'CH3COOH', acceptor = 'SO4-2', product = 'H2S')
```

```
## CH3COOH   SO4-2      H+     H2S     CO2     H2O     HS- 
## -0.1250 -0.1250 -0.1875  0.0625  0.2500  0.2500  0.0625
```

But note that the package ignores elements other than C, H, O, and N in organic compounds, following Rittmann and McCarty (2001, 2020).
For example, if we have some empirical waste biomass formula like this one (spaces added for reading clarity, and could be removed),


``` r
micstoich(donor = 'C226.36 H404.24 O187.69 N1 S0.52', acceptor = 'SO4-2', product = 'H2S')
```

```
## Warning in rxnbal(rtot, tol = tol): Elemental balance is off in reaction:
## -0.00107 C226.36 H404.24 O187.69 N1 S0.52, -0.125 SO4-2, -0.188 H+, 0.0625 H2S,
## 0.242 CO2, 0.00107 NH4+, 0.00107 HCO3-, 0.214 H2O, 0.0625 HS-.
```

```
## C226.36 H404.24 O187.69 N1 S0.52                            SO4-2 
##                     -0.001073768                     -0.125000000 
##                               H+                              H2S 
##                     -0.187500000                      0.062500000 
##                              CO2                             NH4+ 
##                      0.241984323                      0.001073768 
##                            HCO3-                              H2O 
##                      0.001073768                      0.214345538 
##                              HS- 
##                      0.062500000
```

we get an overall reaction, but it is not balanced.
And where is the problem?
We can check.


``` r
r3 <- micstoich(donor = 'C226.36 H404.24 O187.69 N1 S0.52', acceptor = 'SO4-2', product = 'H2S')
```

```
## Warning in rxnbal(rtot, tol = tol): Elemental balance is off in reaction:
## -0.00107 C226.36 H404.24 O187.69 N1 S0.52, -0.125 SO4-2, -0.188 H+, 0.0625 H2S,
## 0.242 CO2, 0.00107 NH4+, 0.00107 HCO3-, 0.214 H2O, 0.0625 HS-.
```

``` r
rxnbal(r3)
```

```
## Warning in rxnbal(r3): Elemental balance is off in reaction: -0.00107 C226.36
## H404.24 O187.69 N1 S0.52, -0.125 SO4-2, -0.188 H+, 0.0625 H2S, 0.242 CO2,
## 0.00107 NH4+, 0.00107 HCO3-, 0.214 H2O, 0.0625 HS-.
```

```
##             C             H             O             N             S 
##  0.0000000000  0.0000000000  0.0000000000  0.0000000000 -0.0005583593
```

The problem is extra S on the reactant side, in the substrate of course.
The negative sign in the results indicates that we have excess reactants, not products.
To use the package as intended, don't include S in the substrate!


``` r
micstoich(donor = 'C226.36 H404.24 O187.69 N1', acceptor = 'SO4-2', product = 'H2S')
```

```
## C226.36 H404.24 O187.69 N1                      SO4-2 
##               -0.001073768               -0.125000000 
##                         H+                        H2S 
##               -0.187500000                0.062500000 
##                        CO2                       NH4+ 
##                0.241984323                0.001073768 
##                      HCO3-                        H2O 
##                0.001073768                0.214345538 
##                        HS- 
##                0.062500000
```

# Methanogenesis

Either use CO2 as the acceptor, but don't specify a product,


``` r
micstoich(donor = 'CH3COOH', acceptor = 'CO2')
```

```
## CH3COOH     CO2     CH4 
##  -0.125   0.125   0.125
```

or specify the product, like this,


``` r
micstoich(donor = 'CH3COOH', product = 'CH4')
```

```
## CH3COOH     CH4     CO2 
##  -0.125   0.125   0.125
```

Or both.


``` r
micstoich(donor = 'CH3COOH', acceptor = 'CO2', product = 'CH4')
```

```
## CH3COOH     CO2     CH4 
##  -0.125   0.125   0.125
```

Hydrogenotrophic methanogenesis has similar flexibility.


``` r
micstoich(donor = 'H2', acceptor = 'CO2')
```

```
##     H2    CO2    H2O    CH4 
## -0.500 -0.125  0.250  0.125
```

``` r
micstoich(donor = 'H2', product = 'CH4')
```

```
##     H2    CO2    CH4    H2O 
## -0.500 -0.125  0.125  0.250
```

``` r
micstoich(donor = 'H2', acceptor = 'CO2', product = 'CH4')
```

```
##     H2    CO2    CH4    H2O 
## -0.500 -0.125  0.125  0.250
```

The different calls do result in different element order in the output.
We could use some of the other ordering options to align them.
The calls below demonstrate this, and have some cell biomass added to make the effects more clear.


``` r
micstoich(donor = 'H2', acceptor = 'CO2', fs = 0.1, arrange = 'decreasing')
```

```
##      H2     CO2    NH4+   HCO3-     H2O     CH4 C5H7O2N 
## -0.5000 -0.1325 -0.0050 -0.0050  0.2700  0.1125  0.0050
```

``` r
micstoich(donor = 'H2', acceptor = 'CO2', fs = 0.1, arrange = 'increasing')
```

```
##    NH4+   HCO3-     CO2      H2 C5H7O2N     CH4     H2O 
## -0.0050 -0.0050 -0.1325 -0.5000  0.0050  0.1125  0.2700
```

``` r
micstoich(donor = 'H2', acceptor = 'CO2', fs = 0.1, arrange = 'alphanum')
```

```
##     CO2      H2   HCO3-    NH4+ C5H7O2N     CH4     H2O 
## -0.1325 -0.5000 -0.0050 -0.0050  0.0050  0.1125  0.2700
```

All options list reactants first and products later.

# Fermentation

For fermentation, specify the organic product.

Glucose to acetic acid.


``` r
micstoich(donor = 'C6H12O6', product = 'CH3COOH')
```

```
##     C6H12O6     CH3COOH 
## -0.04166667  0.12500000
```

Or to ethanol.


``` r
micstoich(donor = 'C6H12O6', product = 'CH3CH2OH')
```

```
##     C6H12O6    CH3CH2OH         CO2 
## -0.04166667  0.08333333  0.08333333
```

Here we get carbon dioxide too.
This is why bread rises.

For production of only H2, use one of these approaches.


``` r
micstoich(donor = 'C6H12O6', acceptor = 'H+')
```

```
##     C6H12O6         H2O         CO2          H2 
## -0.04166667 -0.25000000  0.25000000  0.50000000
```

``` r
micstoich(donor = 'C6H12O6', product = 'H2')
```

```
##     C6H12O6         H2O          H2         CO2 
## -0.04166667 -0.25000000  0.50000000  0.25000000
```

Of course, we have been ignoring cell biomass, but it could be added.


``` r
micstoich(donor = 'C6H12O6', product = 'H2', fs = 0.2)
```

```
##     C6H12O6        NH4+       HCO3-         H2O          H2     C5H7O2N 
## -0.04166667 -0.01000000 -0.01000000 -0.16000000  0.40000000  0.01000000 
##         CO2 
##  0.21000000
```

For mixed fermentations, just combine the products into one formula.
This example follows example 5.6 in Rittmann and McCarty (2020): fermentation of citrate to formate and acetate in a 1:2 molar ratio.
We need uncharged substrates, so will use citric acid.


``` r
micstoich(donor = 'COOHCH2COHCOOHCH2COOH', product = 'HCOOH (CH3COOH)2')
```

```
## COOHCH2COHCOOHCH2COOH                   H2O      HCOOH (CH3COOH)2 
##           -0.05555556           -0.05555556            0.05555556 
##                   CO2 
##            0.05555556
```

# Mixed substrates
Mixed substrates follow the mixed product example.

For example, a mix of lactose and protein (skim milk?) fermented to a mix of lactic and acetic acids plus hydrogen:


``` r
micstoich(donor = '(C12H22O11)5 (C5H7O2N)', product = '(C3H6O3)3 CH3COOH (H2)5')
```

```
##  (C12H22O11)5 (C5H7O2N)                     H2O (C3H6O3)3 CH3COOH (H2)5 
##            -0.003846154            -0.080911681             0.018518519 
##                     CO2                    NH4+                   HCO3- 
##             0.042450142             0.003846154             0.003846154
```

The mixed donor and product are no problem, but we have to manually separate the results.

# Errors in `micstoich()`

The following calls are missing something, and the error messages are meant to explain the problem.


``` r
micstoich(donor = 'C6H12O6')
```

```
## Error in `micstoich()`:
## ! acceptor and product arguments cannot both be missing.
```

``` r
micstoich(donor = 'C6H12O6', acceptor = 'Fe')
```

```
## Error in `if (product %in% pnm) ...`:
## ! argument is of length zero
```

``` r
micstoich(donor = 'C100P', acceptor = 'O2')
```

```
## Warning in rxnbal(rtot, tol = tol): Elemental balance is off in reaction:
## -0.0025 C100P, -0.25 O2, 0.25 CO2.
```

```
##   C100P      O2     CO2 
## -0.0025 -0.2500  0.2500
```

Presently mistakes in element names can give confusing errors.


``` r
micstoich(donor = 'X', acceptor = 'O2')
```

```
## Error in `hrlookup()`:
## ! Problem with reactant  or product X: Not found. Extra space? Choices are: O2 H2O, CO2 CH4, NO3- N2, NO3- NH3, NO3- NH4+, NO3- NO2-, NO2- NH3, NO2- NH4+, NO2- N2, N2 NH3, N2 NH4+, Fe+3 Fe+2, H+ H2, SO4-2 H2S, SO4-2 SO3-2, SO4-2 S, SO4-2 S2O3-2
```

A good way to proceed in case of a strange error message is to try the lower-level functions described below.
For example:


``` r
molmass('X')
```

```
## Error in `molmass()`:
## ! One or more elements in "form" is not in the database.
```

``` r
calcCOD('X')
```

```
## Error in `molmass()`:
## ! One or more elements in "form" is not in the database.
```

The next error is related to how all compounds are identified internally--simply by the character strings provided by formulas (or calculated).


``` r
micstoich(donor = 'C5H6O2N', product = 'C5H6O2N')
```

```
## Error in `micstoich()`:
## ! There are duplicates in input formulas (donor, acceptor, product, or bioform)!
## Change at least one (change in element order is OK)
```

``` r
micstoich(donor = 'C5H7O2N', acceptor = 'O2', bioform = 'C5H7O2N', fs = 0.2)
```

```
## Error in `micstoich()`:
## ! The bioform formula is the donor! Change one (change in element order is OK).
```

They can be fixed by simply changing element order,


``` r
micstoich(donor = 'C5H6O2N', product = 'H6C5O2N')
```

```
##     C5H6O2N     H6C5O2N 
## -0.05263158  0.05263158
```

``` r
micstoich(donor = 'H7C5O2N', acceptor = 'O2', bioform = 'C5H7O2N', fs = 0.2)
```

```
## H7C5O2N      O2 C5H7O2N     CO2    NH4+   HCO3-     H2O 
##   -0.05   -0.20    0.01    0.16    0.04    0.04    0.04
```

or even just adding a space.


``` r
micstoich(donor = 'C5 H6O2N', product = 'C5H6O2N')
```

```
##    C5 H6O2N     C5H6O2N 
## -0.05263158  0.05263158
```

``` r
micstoich(donor = 'C5 H7O2N', acceptor = 'O2', bioform = 'C5H7O2N', fs = 0.2)
```

```
## C5 H7O2N       O2  C5H7O2N      CO2     NH4+    HCO3-      H2O 
##    -0.05    -0.20     0.01     0.16     0.04     0.04     0.04
```

# Molar mass and COD'

For convenience, the package includes a function for calculating molar mass and "calculated" or "theoretical" chemical oxygen demand, COD'.
These are closely related to the functions with similar names in the biogas package.
Any of the chemical formulas used above work here.


``` r
molmass('C6H12O6')
```

```
## [1] 180.156
```

``` r
molmass('C226.36 H404.24 O187.69 N1 S0.52')
```

```
## [1] 6159.797
```

And so do others.


``` r
molmass('S')
```

```
## [1] 32.1
```

``` r
molmass('NaCl')
```

```
## [1] 58.38977
```

``` r
molmass('XePtF6')
```

```
## [1] 440.3694
```

COD' is returned as g oxygen per g of the given compound by default.


``` r
calcCOD('C6H12O6')
```

```
## [1] 1.065743
```

``` r
calcCOD('C226.36 H404.24 O187.69 N1 S0.52')
```

```
## [1] 1.20952
```

But this could be changed to per mole.


``` r
calcCOD('C6H12O6', per = 'mol')
```

```
## [1] 192
```

``` r
calcCOD('C226.36 H404.24 O187.69 N1 S0.52', per = 'mol')
```

```
## [1] 7450.4
```

As expected, COD' will be zero for inorganic compounds.


``` r
calcCOD('NaCl')
```

```
## [1] 0
```

# More information
For more on the micstoich package, see <https://github.com/AU-BCE-EE/micstoich>.

# References
Rittmann, B.E. and McCarty, P.L. (2001) *Environmental Biotechnology:
Principles and Applications*. McGraw-Hill, New York.

Rittmann, B.E. and McCarty, P.L. (2020) *Environmental Biotechnology:
Principles and Applications*. 2nd ed. McGraw-Hill, New York.

