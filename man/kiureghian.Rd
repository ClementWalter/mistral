\name{kiureghian}
\docType{data}
\alias{kiureghian}
\title{A limit-state-function defined by Der Kiureghian}
\description{
  The limit-state function is defined by:
    \deqn{
      f(x) = b - x_2 - \kappa*(x_1-e)^2
    }
  with \eqn{b = 5}, \eqn{\kappa = 0.5} and \eqn{e = 0.1}.
}
\usage{kiureghian}
\format{The function can handle a vector or matrix with column vectors.}
\references{
  Der Kiureghian, A and Dakessian, T:\cr
  \emph{Multiple design points in first and second-order reliability}\cr
  Structural Safety, 20, 1, 37-49, 1998.\cr
}
\keyword{datasets}