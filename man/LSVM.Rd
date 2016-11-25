\name{LSVM}
\alias{LSVM}
\title{Linear Support Vector Machine under monotonicity constraints}

\description{
  Produce a globally increasing binary classifier built from linear monotonic SVM
}

\usage{
  LSVM(x, A.model.lsvm, convexity)
}

\arguments{
  \item{x}{a set of points where the class must be estimated.}
  \item{A.model.lsvm}{a matrix containing the parameters of all hyperplanes.}
  \item{convexity}{Either -1 if the set of data associated to the label "-1" is convex or +1 otherwise.}
}

\value{
  An object of class \code{integer} representing the class of x

  \item{res}{A vector of -1 or +1.}
}

\details{
  LSVM is a monotonic binary classifier built from linear SVM under the constraint that one of the two classes of data is convex.
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
  \code{ \link{modelLSVM} }
}

\examples{
# A limit state function
f <- function(x){  sqrt(sum(x^2)) - sqrt(2)/2 }

# Creation of the data sets
n <- 200
X <- matrix(runif(2*n), nrow = n)
Y <- apply(X, MARGIN = 1, function(w){sign(f(w))})

#The convexity is known
\dontrun{
  model.A <- modelLSVM(X, Y, convexity = -1)
  m <- 10
  X.test <- matrix(runif(2*m), nrow = m)
  classOf.X.test <- LSVM(X.test, model.A, convexity = -1)
}

}


