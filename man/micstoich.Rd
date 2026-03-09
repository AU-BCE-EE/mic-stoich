\name{micstoich}
\alias{micstoich}
\title{Microbial reaction stoichiometry}
\description{
  Calculates stoichiometry of a microbial reaction using the half-reaction
  approach of Rittmann and McCarty (2001). Combines electron donor, acceptor,
  and biomass synthesis half-reactions weighted by the synthesis fraction
  \code{fs}.
}
\usage{
micstoich(donor, acceptor = NULL, product = NULL,
          bioform = "C5H7O2N", fs = 0,
          elements = c("C", "H", "O", "N"),
          order = "rxn", tol = 1e-10)
}
\arguments{
  \item{donor}{Chemical formula of the electron donor as a character string,
    e.g., \code{"C6H12O6"}. For mixed donors, provide a single
    space-separated string, e.g., \code{"CH3COOH H2"}.}
  \item{acceptor}{Chemical formula of the electron acceptor. For inorganic
    acceptors (O2, NO3-, SO4-2, CO2), use the built-in half-reactions by
    passing the acceptor formula alone, e.g., \code{"O2"} or \code{"NO3-"}.
    Can be \code{NULL} for fermentation (requires \code{product}).}
  \item{product}{Chemical formula of the primary reaction product. Required
    for fermentation (no acceptor) or when an acceptor has multiple possible
    products (e.g., NO3- can yield N2, NH3, or NO2-). Also accepts a
    space-separated string for multiple products.}
  \item{bioform}{Chemical formula of the biomass, default \code{"C5H7O2N"}.}
  \item{fs}{Fraction of electrons directed to biomass synthesis (0 to 1),
    default \code{0} (no synthesis). The energy fraction is \code{fe = 1 - fs}.}
  \item{elements}{Character vector of elements to track in the output,
    default \code{c("C", "H", "O", "N")}.}
  \item{order}{Controls ordering of compounds in the output vector.
    \code{"rxn"} (default) orders as reactants then products, with donor and
    acceptor first among reactants and \code{product} and \code{bioform} first
    among products. \code{"sort"} sorts by descending absolute coefficient,
    reactants first. \code{NULL} or \code{NA} returns compounds in natural
    output order. A character vector of compound names gives a custom order.}
  \item{tol}{Numerical tolerance: coefficients with absolute value below
    \code{tol} are set to zero and dropped from the output. Default
    \code{1e-10}.}
}
\value{
  A named numeric vector of stoichiometric coefficients. Negative values are
  reactants, positive values are products. Coefficients are normalised per
  electron equivalent.
}
\references{
  Rittmann, B.E. and McCarty, P.L. (2001) \emph{Environmental Biotechnology:
  Principles and Applications}. McGraw-Hill, New York.
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

# Denitrification: specify product because NO3- has multiple options
micstoich(donor = "C6H12O6", acceptor = "NO3-", product = "N2")
}
