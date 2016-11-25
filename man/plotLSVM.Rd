\name{plotLSVM}
\alias{plotLSVM}
\title{plot of LSVM}

\description{
  Make a plot of the data and the LSVM classifier
}

\usage{
  plotLSVM(X,
          Y,
          A.model.lsvm,
          hyperplanes = FALSE,
          limit.state.estimate = TRUE,
          convexity)
}

\arguments{
  \item{X}{a matrix containing the data sets}
  \item{Y}{a vector containing -1 or +1 that reprensents the class of each elements of X.}
  \item{A.model.lsvm}{a matrix containing the parameters of all hyperplanes.}
  \item{hyperplanes}{A boolean. If TRUE, plot the hyperplanes obtained.}
  \item{limit.state.estimate}{A boolean. If TRUE, plot the estimate of the limit state.}
  \item{convexity}{Either -1 if the set of data associated to the label "-1" is convex or +1 otherwise.}
}

\details{
  plotLSVM makes a plot of the data as well as the estimate limit state and the hyperplanes involved in this construction.
}

\note{
  This function is useful only in dimension 2.
}

\references{
  \itemize{
    \item
      R.T. Rockafellar:\cr
      \emph{Convex analysis}\cr
      Princeton university press, 2015.\cr
    \item
      N. Bousquet, T. Klein and V. Moutoussamy :\cr
      \emph{Approximation of limit state surfaces in monotonic Monte Carlo settings}\cr
      Submitted .\cr
  }
}

\author{
  Vincent Moutoussamy\cr
}

\seealso{
  \code{\link{LSVM}}
  \code{\link{modelLSVM}}
}

\examples{

# A limit state function
f <- function(x){  sqrt(sum(x^2)) - sqrt(2)/2 }

# Creation of the data sets

n <- 200
X <- matrix(runif(2*n), nrow = n)
Y <- apply(X, MARGIN = 1, function(w){sign(f(w))})

\dontrun{
  model.A <- modelLSVM(X,Y, convexity = -1)
  plotLSVM(X, Y, model.A, hyperplanes = FALSE, limit.state.estimate = TRUE, convexity = -1)
}

}

