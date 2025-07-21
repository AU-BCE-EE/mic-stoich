# Fermentation stoichiometry from Rittmann and McCarty

# Function to get stoichiometry for custom organic reaction (O-19)
customOrgStoich <- function(
  form, 
  elements =  c('C', 'H', 'O', 'N'),
  dover = FALSE
  ) {
  
  fc <- readFormula(form, elements)

  # Use symbols from O-19 in R&M
  n <- as.numeric(fc['C'])
  a <- as.numeric(fc['H'])
  b <- as.numeric(fc['O'])
  cc <- as.numeric(fc['N'])
  
  # dover for d override deals with case of no available electrons
  if (dover) {
    d <- 1
    H.c <- 0
  } else {
    d <- 4 * n + a - 2 * b - 3 * cc
    H.c <- -1
  }

  if (d == 0) {
    return(form)
  }
  
  # Put together
  #rr <- c(CO2 = - (n - cc) / d, NH4. = - cc / d, HCO3. = - cc / d, H. = -1, H2O = (2*n - b + cc) /d)
  rr <- c(CO2 = - n / d, NH3 = - cc / d, H. = H.c, H2O = (2*n - b) / d)
  rr[form] <-  1/d

  return(rr)
  
}

# General Rittmann and McCarty stype stoichiometry calculations
RMStoich <- function(donor, rd, ra, rc, fs, dropzero, dropsub, order, tol) {

  ii <- unique(names(c(rd, rc, ra)))

  # Blanks
  rd[ii[!ii %in% names(rd)]] <- 0
  rc[ii[!ii %in% names(rc)]] <- 0
  ra[ii[!ii %in% names(ra)]] <- 0

  # Order
  rd <- rd[ii]
  rc <- rc[ii]
  ra <- ra[ii]

  fe <- 1 - fs
  
  # Combine
  rtot <- fe * ra + fs * rc  - rd

  # Drop substrate if requested
  if (isTRUE(dropsub)) {
    # Adjust coefficients to 1 mol substrate
    rtot <- - rtot / rtot[donor]
    rtot <- rtot[names(rtot) != donor] 
  }

  rtot[abs(rtot) < tol] <- 0 
  
  # Drop empty elements
  if (dropzero) {
    rtot <- rtot[rtot != 0]
  }

  if (!is.na(order[1]) && tolower(order[1]) == 'sort') {
    rtot <- rtot[order(rtot < 0, abs(rtot), decreasing = TRUE)]
  } else if (!is.na(order[1]) && all(sort(order) == sort(names(rtot)))) {
    rtot <- rtot[order]
  } else if (inherits(order, 'logical') && isFALSE(order)) {
    # Skips order
  } else if (!is.na(order[1])) {
    warning('order argument ignored')
  }

  return(rtot)

}

predFerm <- function(
  donor,                     # Character chemical formula of electron donor or substrate
  acceptor,                  # Character chemical formula of electron acceptor or fermentation product
  biomassform = 'C5H7O2N',   # Biomass empirical formula
  fs = 0,                    # Fraction substrate going to cell synthesis, f_s in Rittmann and McCarty
  elements = c('C', 'H', 'O', 'N'),
  order = 'sort',
  dropzero = TRUE,
  dropsub = FALSE,
  tol = 1E-10
  ) {

  if (biomassform %in% donor) {
    stop('The biomassform formula is the/a substrate! Change one (change in just order is OK).')
  }

  # Donor half reaction
  rd <- customOrgStoich(donor, elements = elements)

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reactions
  ra <- customOrgStoich(acceptor, elements = elements)

  rtot <- RMStoich(donor = donor, rd = rd, ra = ra, rc = rc, fs = fs, 
                   dropzero = dropzero, dropsub = dropsub, order = order, 
                   tol = tol)

  return(rtot)

}


