\name{halfrxns}
\alias{halfrxns}
\title{List available built-in half-reactions}
\description{
  Returns a data frame of the inorganic acceptor half-reactions built into
  \code{\link{micstoich}}. Use this to check the exact acceptor and product
  names required by the \code{acceptor} and \code{product} arguments.
}
\usage{
halfrxns()
}
\value{
  A data frame with one row per built-in half-reaction and three columns:
  \describe{
    \item{acceptor}{Electron acceptor formula, as used in the \code{acceptor}
      argument of \code{\link{micstoich}}.}
    \item{product}{Reduction product formula, as used in the \code{product}
      argument of \code{\link{micstoich}}.}
    \item{reaction}{Half-reaction stoichiometry as a character string, with
      coefficients per electron equivalent (rounded decimals).}
  }
}
\seealso{
  \code{\link{micstoich}}
}
\examples{
halfrxns()
}
