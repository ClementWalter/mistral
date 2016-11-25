\name{updateLSVM}
\alias{updateLSVM}
\title{Update LSVM classifier}

\description{
  Update the existing classifier LSVM with a new set of data.
}

\usage{
  updateLSVM(X.new,
             Y.new,
             X,
             Y,
             A.model.lsvm,
             convexity,
             PLOTSVM = FALSE,
             step.plot.LSVM = 1,
             hyperplanes = FALSE,
             limit.state.estimate = TRUE)
}

\arguments{
  \item{X.new}{a matrix containing a new data sets}
  \item{Y.new}{a vector containing -1 or +1 that reprensents the class of each elements of X.new.}
  \item{X}{a matrix containing the data sets}
  \item{Y}{a vector containing -1 or +1 that reprensents the class of each elements of X.}
  \item{A.model.lsvm}{a matrix containing the parameters of all hyperplanes.}
  \item{convexity}{Either -1 if the set of data associated to the label "-1" is convex or +1 otherwise.}
  \item{PLOTSVM}{A boolean. If TRUE, plot the data.}
  \item{step.plot.LSVM}{A plot is made each \code{step.plot.LSVM} steps.}
  \item{hyperplanes}{A boolean. If TRUE, plot the hyperplanes obtained.}
  \item{limit.state.estimate}{A boolean. If TRUE, plot the estimate of the limit state.}
}

\value{
  An object of class \code{matrix} containing the parameters of a set of hyperplanes
}

\details{
  updateLSVM allows to make an update of the classifier LSVM.
}

\note{
  The argument PLOTSVM is useful only in dimension 2.
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
  M <- 20
  X.new <- matrix(runif(2*M), nrow = M)
  Y.new <- apply(X.new, MARGIN = 1, function(w){ sign(f(w))})

  X.new.S <- X.new[which(Y.new > 0), ]
  Y.new.S <- Y.new[which(Y.new > 0)]
  model.A.new <- updateLSVM(X.new.S, Y.new.S, X, Y,
                            model.A, convexity = -1, PLOTSVM = TRUE, step.plot.LSVM = 5)
}

}











