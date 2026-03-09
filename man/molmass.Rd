\name{molmass}
\alias{molmass}
\title{Molar mass of a compound}
\description{
  Calculates the molar mass of one or more compounds from their chemical
  formulas.
}
\usage{
molmass(form, elements = NULL)
}
\arguments{
  \item{form}{Character vector of chemical formulas, e.g.,
    \code{c("H2O", "C6H12O6")}.}
  \item{elements}{Optional character vector of element symbols. If provided,
    only the mass contribution from those elements is included, e.g.,
    \code{"N"} to get g N per mol compound.}
}
\value{
  Numeric vector of molar masses in g/mol, one value per formula in
  \code{form}.
}
\seealso{
  \code{biogas::molMass()}, which provides similar functionality.
}
\examples{
molmass("H2O")
molmass(c("C6H12O6", "CH4", "CO2"))

# g N per mol
molmass("C5H7O2N", elements = "N")
}
