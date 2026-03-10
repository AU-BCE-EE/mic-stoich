# Names are 'acceptor product' (space between)
halfrxn <- list(

  'O2 H2O'     = c(O2 = -1/4,    `H+` = -1,     H2O = 1/2),

  'CO2 CH4'    = c(CO2 = -1/8,   `H+` = -1,     CH4 = 1/8,    H2O = 1/4),

  'NO3- N2'    = c(`NO3-` = -1/5,  `H+` = -6/5,   N2 = 1/10,   H2O = 3/5),
  'NO3- NH3'   = c(`NO3-` = -1/8,  `H+` = -5/4,   `NH4+` = 1/8, H2O = 3/8),
  'NO3- NH4+'  = c(`NO3-` = -1/8,  `H+` = -5/4,   `NH4+` = 1/8, H2O = 3/8),              # I-3
  'NO3- NO2-'  = c(`NO3-` = -1/2,  `H+` = -1,     `NO2-` = 1/2, H2O = 1/2),

  'NO2- NH3'   = c(`NO2-` = -1/6,  `H+` = -4/3,   `NH4+` = 1/6, H2O = 1/3),              # I-2
  'NO2- NH4+'  = c(`NO2-` = -1/6,  `H+` = -4/3,   `NH4+` = 1/6, H2O = 1/3),              # I-2
  'NO2- N2'    = c(`NO2-` = -1/3,  `H+` = -4/3,   N2 = 1/6,    H2O = 2/3),

  'N2 NH3'     = c(`N2` = -1/6,  `H+` = -4/3,   `NH4+` = 1/3),                           # I-3
  'N2 NH4+'    = c(`N2` = -1/6,  `H+` = -4/3,   `NH4+` = 1/3),                           # I-3

  'SO4-2 H2S'  = c(`SO4-2` = -1/8, `H+` = -19/16, H2S = 1/16,  `HS-` = 1/16, H2O = 1/2)

)

halfrxns <- function() {
  data.frame(
    acceptor = sapply(strsplit(names(halfrxn), ' '), `[[`, 1),
    product  = sapply(strsplit(names(halfrxn), ' '), `[[`, 2),
    reaction = sapply(halfrxn, function(rxn) paste0(paste0(signif(rxn, 5), names(rxn)), collapse = ' ')),
    row.names = NULL
  )
}
