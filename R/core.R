micstoich <- function(
  donor,
  acceptor = NULL,
  product = NULL,
  bioform = 'C5H7O2N',
  fs = 0,
  elements = c('C', 'H', 'O', 'N'),
  order = 'rxn',
  tol = 1E-10
) {

  if (is.null(acceptor) && is.null(product)) {
    stop('acceptor and product arguments cannot both be missing.')
  }

  if (bioform %in% donor) {
    stop('The bioform formula is the donor! Change one (change in just order is OK).')
  }

  acceptor_orig <- acceptor
  if (!is.null(product)) {
    acceptor <- paste(acceptor, product)
  }

  # Half reactions
  # Donor
  rd <- orgrxn(donor, elements = elements)
  # Acceptor
  if (is.null(product) || !grepl('C|H', product)) {
    # Trim half reaction names to length of acceptor
    substr(names(halfrxn), 1, nchar(acceptor))
    anm <- as.character(lapply(strsplit(names(halfrxn), ' '), `[[`, 1))
    pnm <- as.character(lapply(strsplit(names(halfrxn), ' '), `[[`, 2))
    fnm <- names(halfrxn)
    if (acceptor %in% anm || acceptor %in% fnm) {
      # If in only one anm, use that
      if (sum(acceptor == anm) == 1) {
        ra <- halfrxn[[names(halfrxn)[acceptor == anm]]]
      } else if (sum(acceptor == fnm) == 1) {
        ra <- halfrxn[[names(halfrxn)[acceptor == fnm]]]
      } else {
        stop(paste0('product needed. More than one for this acceptor. Pairs are: ', paste(names(halfrxn), collapse = ', ')))
      }
    } else {
      stop(paste0('Problem with acceptor argument: Not found. Extra space? Choices are: ', paste(names(halfrxn), collapse = ', ')))
    }
  } else {
    ra <- orgrxn(product, elements = elements)
  }

  # Synthesis
  rc <- orgrxn(bioform, elements = elements)

  # N source adjustment
  # NTS: need to figure out!

  # Align
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

  rtot[abs(rtot) < tol] <- 0
  rtot <- rtot[rtot != 0]

  if (!is.null(order[1]) && !is.na(order[1])) {
    if (tolower(order[1]) == 'sort') {
      rtot <- rtot[order(rtot < 0, abs(rtot), decreasing = TRUE)]
    } else if (tolower(order[1]) == 'rxn') {
      # First donor and acceptor
      it <- which(names(rtot) %in% c(donor, acceptor_orig))
      rord <- rtot[it]
      rtot <- rtot[-it]

      # Then other reactants
      it <- which(rtot < 0)
      if (length(it) > 0) {
        rord <- c(rord, rtot[it])
        rtot <- rtot[-it]
      }

      # Then product
      if (!is.null(product)) {
        it <- which(names(rtot) == product)
        rord <- c(rord, rtot[it])
        rtot <- rtot[-it]
      }

      # And biomass
      if (!is.null(bioform) && fs > 0) {
        it <- which(names(rtot) == bioform)
        rord <- c(rord, rtot[it])
        rtot <- rtot[-it]
      }

      # Then remaining products
      rtot <- c(rord, rtot)

    } else {
      warning('order argument ignored. Options are: "sort", "rxn".')
    }
  }

  return(rtot)

}

# Function to get stoichiometry for custom organic reaction (O-19)
orgrxn <- function(
  form, 
  elements =  c('C', 'H', 'O', 'N'),
  dover = FALSE
  ) {
  
  fc <- readform(form, elements)

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
  rr <- c(CO2 = - (n - cc)/d, `NH4+` = -cc/d, `HCO3-` = -cc/d, `H+` = H.c, H2O = (2*n - b + cc)/d)
  #rr <- c(CO2 = - n / d, NH3 = - cc / d, `H+` = H.c, H2O = (2*n - b) / d)
  rr[form] <-  1/d

  return(rr)
  
}

readform <- function(
  form,
  elements = NULL,        # Set of elements returned, all others ignored, e.g., c('C', 'H', 'N', 'O')
  min.elements = NULL,    # Minimum set of elements, will return error if these at least are not included
  cdigits = 6,
  value = 'numeric'       # Type of output, 'numeric' for named vector, 'shortform' for shortened formula
) {

  form.orig <- form

  # Remove spaces
  form <- gsub(' ', '', form)

  # Add implied coefficients of 1 (also after ")")
  form <- gsub('([a-zA-Z\\)])([A-Z\\)\\(])', '\\11\\2', form)
  form <- gsub('([a-zA-Z\\)])([A-Z\\)\\(])', '\\11\\2', form) # Repeated for e.g., COOH
  form <- gsub('([a-zA-Z\\)])$', '\\11', form)

  # Find parentheses and remove them, multipying coefficients inside by coefficient at end
  # So (CH2)2 ---> C2H4
  # First add ( after N), e.g., (CH2)2CH3 ---> (CH2)2(CH3 for separation below
  form <- gsub('(\\)[0-9\\.]+)', '\\1(', form)
  # Drop extra (
  form <- gsub('^\\(', '', form)
  form <- gsub('\\($', '', form)
  form <- gsub('\\(\\(', '(', form)
  s1  <- strsplit(form, '\\(')[[1]]

  # Build up elementwise formula piecewise
  formpw <- NULL
  for(i in seq_along(s1)) {
    xx <- s1[i]
    if(grepl('\\)', xx)) {
      nn <- as.numeric(gsub('.+\\)','',xx))
      ff <- gsub('\\).+', '', xx)
      cc <- nn*as.numeric(strsplit(ff, '[A-Za-z]+')[[1]][-1])
      ee <- strsplit(ff, '[0-9.]+')[[1]]
      formpw <- paste0(formpw, paste0(ee, cc, collapse = ''))
    } else {
      formpw <- paste0(formpw, xx)
    }
  }

  form <- formpw

  # Extract integer coefficients 
  cc <- as.numeric(strsplit(form, '[A-Za-z]+')[[1]][-1])
  names(cc) <- strsplit(form, '[0-9.]+')[[1]]

  # Sort out elements to return
  if(is.null(elements)) elements <- unique(names(cc))
  fc <- numeric(length(elements))
  names(fc) <- elements

  # Fill in fc, summing elements of cc if required (if elements are repeated)
  for(i in elements) {
    for(j in 1:length(cc)) {
      if(names(cc)[j]==i) fc[i] <- fc[i] + cc[j]
    }
  }

  # Simplify form based on fc (for output only)
  # format() is for the rare case with something like C0.00001, to avoid C1e-5 which will result in an error
  form <- paste0(names(fc), format(signif(fc/min(fc), cdigits), scientific = FALSE), collapse = '')
  # Drop spaces that come in with format(...scientific = FALSE) (scipen fix)
  form <- gsub(' ', '', form)
  # And drop coefficients of 1
  form <- gsub('([a-zA-Z])1([a-zA-Z])', '\\1\\2', form)
  form <- gsub('([a-zA-Z])1$', '\\1', form)
  
  # Check for minimum set of elements
  if(!is.null(min.elements)) if(any(!min.elements %in% names(fc)) | any(fc[min.elements] == 0)) stop('Minimum elements required are ', min.elements, ' (from min.elements argument), but form is ', form.orig, ', interpreted as ', form)

  if(value == 'numeric') return(fc)
  if(value == 'shortform') return(as.vector(form))

}

molmass <- function(form, elements = NULL) {

  # Loop through all elements in form
  mmass <- NULL
  for(f in form) {
    # If and only if first letter of form is lowercase, entire string is capitalized
    if(grepl('^[a-z]', f)) f <- toupper(f) 

    # Get coefficients of formula
    fc <- readform(f)

    if (!is.null(elements)) {
      fc <- fc[intersect(names(fc), elements)]
    }

    # Check for unidentified element
    if(any(!names(fc) %in% names(atomweights))) stop('One or more elements in \"form\" is not in the database.')

    # Calculate molar mass, using names of fc for indexing
    mmass <- c(mmass, sum(atomweights[names(fc)]*fc))
  }

  return(mmass)
}

# COD' as g O per mol substrate
calcCOD <- function(form, per = 'g') {

  # If and only if first letter of form is lowercase, entire string is capitalized
  if(grepl('^[a-z]', form)) {
    form <- toupper(form)
  }

  # Read formula (function not vectorized)
  fc <- readform(form, elements = c('C', 'H', 'O', 'N'))
  # Calculate COD based on Rittmann and McCarty
  COD <- as.vector((2*fc['C'] + 0.5*fc['H'] - 1.5*fc['N'] - fc['O']) * molmass('O'))

  if (per == 'g') {
    COD <- COD / molmass(form)
  } else if (per != 'mol') {
    stop(paste('Argument per must be \"g\" or \"mol\" but is', per))
  }

  return(COD)

}

# Get mass conversion factor to go from moles of component to mass COD, N, C, S, or total, in that order
massconv <- function(form) {

  # Remove p and m (+/-)
  form <- gsub('p$|m$', '', form)
  
  cod <- calcCOD(form)
  fn <- readform(form)
  
  if (cod > 0) {
    cf <- cod
  } else if ('N' %in% names(fn)) {
    cf <- molmass(form, elements = 'N')
  } else if ('C' %in% names(fn)) {
    cf <- molmass(form, elements = 'C')
  } else if ('S' %in% names(fn)) {
    cf <- molmass(form, elements = 'S')
  } else {
    cf <- molmass(form)
  }

  return(cf)

}
