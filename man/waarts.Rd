\name{waarts}
\docType{data}
\alias{waarts}
\title{A limit-state-function defined by Waarts}
\description{
  The limit-state function is defined by:
    \deqn{ b1 = 3 + (u_1-u_2)^2/10 - sign(u_1 + u_2)*(u_1+u_2)/sqrt(2) }
    \deqn{ b2 = sign(u_2-u_1)*(u_1-u_2)+7/sqrt(2) }
    \deqn{ f(u) = min(b1, b2) }
}
\usage{waarts}
\format{The function can handle a vector or matrix with column vectors.}
\references{
  Waarts, PH:\cr
  \emph{An appraisal of DARS: directional adaptive response surface sampling}\cr
  Delft University Press, The Netherlands, 2000.\cr
}
\keyword{datasets}