\name{calcCOD}
\alias{calcCOD}
\title{Chemical oxygen demand of a compound}
\description{
  Calculates the theoretical chemical oxygen demand (COD) of a compound from
  its chemical formula, using the approach of Rittmann and McCarty (2001).
}
\usage{
calcCOD(form, per = "g")
}
\arguments{
  \item{form}{Chemical formula as a character string, e.g., \code{"C6H12O6"}.}
  \item{per}{Units for the result: \code{"g"} (default) for g O2 per g
    compound, or \code{"mol"} for g O2 per mol compound.}
}
\details{
  Similar to \code{biogas::calcCOD()} but also supports mole-based output
  via the \code{per} argument. Note that \code{biogas::calcCOD()} returns
  g O2 per g compound; the default here is the same.
}
\value{
  Numeric value of COD in the units specified by \code{per}.
}
\seealso{
  \code{biogas::calcCOD()}
}
\references{
  Rittmann, B.E. and McCarty, P.L. (2001) \emph{Environmental Biotechnology:
  Principles and Applications}. McGraw-Hill, New York.
}
\examples{
# g O2 per g glucose
calcCOD("C6H12O6")

# g O2 per mol glucose
calcCOD("C6H12O6", per = "mol")

# g O2 per g methane
calcCOD("CH4")
}
