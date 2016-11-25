\name{modelLSVM}
\alias{modelLSVM}
\title{Estimation of the parameters of the LSVM}

\description{
  Produce a matrix containing the parameters of a set of hyperplanes separating the two classes of data
}

\usage{
  modelLSVM(X, Y, convexity)
}

\arguments{
  \item{X}{a matrix containing the data sets}
  \item{Y}{a vector containing -1 or +1 that reprensents the class of each elements of X.}
  \item{convexity}{Either -1 if the set of data associated to the label "-1" is convex or +1 otherwise.}
}

\value{
  An object of class \code{matrix} containing the parameters of a set of hyperplanes

  \item{res}{A matrix where each lines contains the parameters of a hyperplane.}
}

\details{
  modelLSVM evaluate the classifier on a set of points.
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
}

}


