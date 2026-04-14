# Names are 'acceptor product' (space between)
halfrxn <- list(

  'O2 H2O'       = c(O2 = -1/4,    `H+` = -1,     H2O = 1/2),                              # I-14

  'CO2 CH4'      = c(CO2 = -1/8,   `H+` = -1,     CH4 = 1/8,    H2O = 1/4),                # O-12

  'NO3- N2'      = c(`NO3-` = -1/5,  `H+` = -6/5,   N2 = 1/10,   H2O = 3/5),               # I-7
  'NO3- NH3'     = c(`NO3-` = -1/8,  `H+` = -5/4,   `NH4+` = 1/8, H2O = 3/8),              # I-1
  'NO3- NH4+'    = c(`NO3-` = -1/8,  `H+` = -5/4,   `NH4+` = 1/8, H2O = 3/8),              # I-1
  'NO3- NO2-'    = c(`NO3-` = -1/2,  `H+` = -1,     `NO2-` = 1/2, H2O = 1/2),              # I-6

  'NO2- NH3'     = c(`NO2-` = -1/6,  `H+` = -4/3,   `NH4+` = 1/6, H2O = 1/3),              # I-2
  'NO2- NH4+'    = c(`NO2-` = -1/6,  `H+` = -4/3,   `NH4+` = 1/6, H2O = 1/3),              # I-2
  'NO2- N2'      = c(`NO2-` = -1/3,  `H+` = -4/3,   N2 = 1/6,    H2O = 2/3),               # I-8

  'N2 NH3'       = c(`N2` = -1/6,  `H+` = -4/3,   `NH4+` = 1/3),                           # I-3
  'N2 NH4+'      = c(`N2` = -1/6,  `H+` = -4/3,   `NH4+` = 1/3),                           # I-3

  'Fe+3 Fe+2'    = c(`Fe+3` = -1, `Fe+2` = 1),                                             # I-4

  'MnO2 Mn+2'    = c(MnO2 = -1/2, `H+` = -2, `Mn+2` = 1/2, H2O = 1),                      # Mn(IV) reduction

  'H+ H2'        = c(`H+` = -1, `H2` = 1/2),                                               # I-5

  'SO4-2 H2S'    = c(`SO4-2` = -1/8, `H+` = -19/16, H2S = 1/16,  `HS-` = 1/16, H2O = 1/2), # I-9
  'SO4-2 SO3-2'  = c(`SO4-2` = -1/2, `H+` = -1, `SO3-2` = 1/2,  H2O = 1/2),                # I-11
  'SO4-2 S'      = c(`SO4-2` = -1/6, `H+` = -4/3, `S` = 1/6,  H2O = 2/3),                  # I-12
  'SO4-2 S2O3-2' = c(`SO4-2` = -1/4, `H+` = -5/4, `S2O3-2` = 1/8,  H2O = 5/8),            # I-13

  'S H2S'        = c(S = -1/2, `H+` = -1, H2S = 1/2),                                      # Elemental S reduction

  # Bioremediation
  'HAsO4-2 H3AsO3' = c(`HAsO4-2` = -1/2, `H+` = -2,   H3AsO3 = 1/2,  H2O = 1/2),          # As(V) reduction
  'ClO4- Cl-'      = c(`ClO4-` = -1/8,   `H+` = -1,   `Cl-` = 1/8,   H2O = 1/2),           # Perchlorate reduction
  'CrO4-2 Cr+3'    = c(`CrO4-2` = -1/3,  `H+` = -8/3, `Cr+3` = 1/3,  H2O = 4/3),           # Chromate reduction
  'UO2+2 U+4'      = c(`UO2+2` = -1/2,   `H+` = -2,   `U+4` = 1/2,   H2O = 1),             # Uranium reduction
  'SeO4-2 Se'      = c(`SeO4-2` = -1/6,  `H+` = -4/3, Se = 1/6,      H2O = 2/3)             # Selenate reduction

)

halfrxns <- function() {
  data.frame(
    acceptor = sapply(strsplit(names(halfrxn), ' '), `[[`, 1),
    product  = sapply(strsplit(names(halfrxn), ' '), `[[`, 2),
    reaction = sapply(halfrxn, function(rxn) paste0(paste0(signif(rxn, 5), names(rxn)), collapse = ' ')),
    row.names = NULL
  )
}
