\name{micstoich}
\alias{micstoich}
\title{Microbial reaction stoichiometry}
\description{
  Determines stoichiometry of a microbial reaction using the half-reaction
  approach of Rittmann and McCarty (2001, 2020). Combines electron donor, acceptor,
  and biomass synthesis half-reactions weighted by the synthesis fraction
  \code{fs}. 
}
\usage{
micstoich(donor, acceptor = NULL, product = NULL,
          bioform = "C5H7O2N", fs = 0,
          order = "rxn", tol = 1e-10)
}
\arguments{
  \item{donor}{Chemical formula of the electron donor as a character string,
    e.g., \code{"C6H12O6"} for glucose. Can handle decimal coefficients
    (e.g., for empirical formulas) and mixed donors. Only C, H, O, and N are
    tracked; other elements in the formula are ignored.}
  \item{acceptor}{Chemical formula of the electron acceptor. For inorganic
    acceptors (O2, NO3-, SO4-2, CO2), use the built-in half-reactions by
    passing the acceptor formula, e.g., \code{"O2"} or \code{"NO3-"}.
    Omit for fermentation and provide \code{product} instead.}
  \item{product}{Chemical formula of the primary reaction product. Required
    for fermentation (no acceptor) or when an acceptor has multiple possible
    products (e.g., NO3- can yield N2, NH3, or NO2-).}
  \item{bioform}{Chemical formula of cell biomass, default \code{"C5H7O2N"}.}
  \item{fs}{Fraction of electrons directed to biomass synthesis (0 to 1),
    default \code{0} (no synthesis). The energy fraction is \code{fe = 1 - fs}.}
  \item{order}{Controls ordering of compounds in the output vector.
    \code{"rxn"} (default) orders as reactants then products, with \code{donor} and
    \code{acceptor} (if any) first among reactants and \code{product} and \code{bioform} first
    among products. \code{"sort"} sorts by descending absolute coefficient,
    reactants first. \code{NULL} or \code{NA} returns compounds in order they
    are returned by internal functions.}
  \item{tol}{Numerical tolerance: coefficients with absolute value below
    \code{tol} are set to zero and dropped from the output. Default
    \code{1e-10}.}
}
\value{
  A named numeric vector of stoichiometric coefficients. Names are reactants 
  and products. Negative values are reactants, positive values are products. 
  Coefficients are normalised to one electron equivalent.
}
\references{
  Rittmann, B.E. and McCarty, P.L. (2001) \emph{Environmental Biotechnology:
  Principles and Applications}. McGraw-Hill, New York.

  Rittmann, B.E. and McCarty, P.L. (2020) \emph{Environmental Biotechnology:
  Principles and Applications}. 2nd ed. McGraw-Hill, New York.
}

\examples{
# Aerobic respiration of glucose, no biomass
micstoich(donor = "C6H12O6", acceptor = "O2")

# With biomass synthesis (10\% of electrons to synthesis)
micstoich(donor = "C6H12O6", acceptor = "O2", fs = 0.1)

# Methanogenesis
micstoich(donor = "C6H12O6", acceptor = "CO2")

# Fermentation to ethanol
micstoich(donor = "C6H12O6", product = "CH3CH2OH")

# Denitrification: specify product because NO3- has multiple product options
micstoich(donor = "C6H12O6", acceptor = "NO3-", product = "N2")
}
