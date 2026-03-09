\name{rxnbal}
\alias{rxnbal}
\title{Check elemental balance of a reaction}
\description{
  Checks whether a reaction is elementally balanced by summing element counts
  weighted by stoichiometric coefficients.
}
\usage{
rxnbal(rxn, tol = 1e-10)
}
\arguments{
  \item{rxn}{A named numeric vector of stoichiometric coefficients, where
    names are compound formulas. Negative values are reactants, positive
    values are products, as returned by \code{\link{micstoich}}.}
  \item{tol}{Numerical tolerance for balance check. Default \code{1e-10}.}
}
\value{
  Returns \code{TRUE} invisibly if the reaction is balanced. If not, issues
  a warning and returns a named numeric vector of elemental imbalances (the
  sum of element counts weighted by stoichiometric coefficients; zero means
  balanced for that element).
}
\seealso{
  \code{\link{micstoich}}
}
\examples{
# A balanced reaction should return TRUE invisibly
rxn <- micstoich(donor = "C6H12O6", acceptor = "O2")
rxnbal(rxn)

# Check the return value explicitly
isTRUE(rxnbal(rxn))
}
