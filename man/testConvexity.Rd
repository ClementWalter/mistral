\name{testConvexity}
\alias{testConvexity}
\title{Test the convexity of set of data}

\description{
  Provides the 
}

\usage{
  testConvexity(X,Y)
}

\arguments{
  \item{X}{a matrix containing the data sets}
  \item{Y}{a vector containing -1 or +1 that reprensents the class of each elements of X.}
}

\value{
  An object of class \code{list} containing the number of the class which is convex and the parameters of a set of hyperplanes separating the two classes
}

\details{
  testConvexity test if one of the two data set is potentially convex.
}


\references{
  \itemize{
    \item
      R.T. Rockafellar:\cr
      \emph{Convex analysis}\cr
      Princeton university press, 2015.\cr
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
  TEST.Convexity <- testConvexity(X, Y)
  if(length(TEST.Convexity) == 2){
    Convexity <- TEST.Convexity[[1]] 
    model.A   <- TEST.Convexity[[2]]
  }
  if(length(TEST.Convexity) == 1){
    # The problem is not convex
    Convexity <- 0 #the problem is not convex
  }
}

}











