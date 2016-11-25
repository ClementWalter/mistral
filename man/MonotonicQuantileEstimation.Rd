\name{MonotonicQuantileEstimation}
\alias{MonotonicQuantileEstimation}
\title{Quantile estimation under monotonicity constraints}

\description{
  Estimate a quantile with the constraints that the function is monotone
}

\usage{
  MonotonicQuantileEstimation(f, 
                            inputDimension, 
                            inputDistribution,
                            dir.monot,
                            N.calls,
                            p,
                            method,
                            X.input = NULL,
                            Y.input = NULL)
}

\arguments{
 \item{f}{
    a failure fonction
  }
  \item{inputDimension}{
    dimension of the inputs
  }
  \item{inputDistribution}{
    a list of length \sQuote{inputDimension} which contains the name of the input distribution and their parameters. For the input "i",
    inputDistribution[[i]] = list("name_law",c(parameters1,..., parametersN))
  }
  \item{dir.monot}{
    vector of size \code{inputDimension} which represents the monotonicity of the failure function. dir.monot[i] = -1 (resp. 1) if the failure function f is decreasing (resp. increasing) according with direction i.
  }
  \item{N.calls}{
     Number of calls to f allowed
  }
  \item{method}{
    there are four methods available. "MonteCarloWB" provides the empirical quantile estimator, "MonteCarloWB" provides the empirical quantile estimator as well as two bounds for the searched quantile, "Bounds" provides two bounds for a quantile from a set of points and "MonteCarloIS" provides an estimate of a quantile based on a sequential framework of simulation.
  }
  \item{p}{
    the probability associated to the quantile
  }
  \item{X.input}{
    a set of points
  }
  \item{Y.input}{
    value of f on X.input
  }
}

\value{
  An object of class \code{list} containing the quantile as well as:

  \item{qm}{A lower bound of the quantile.}
  \item{qM}{A upperer bound of the quantile.}
  \item{q.hat}{An estimate of the quantile.}
  \item{Um}{A lower bounds of the probability obtained from the desing of experiments.}
  \item{UM}{An upper bounds of the probability obtained from the desing of experiments.}
  \item{XX}{Design of experiments}
  \item{YY}{Values of on XX}
}

\details{
  MonotonicQuantileEstimation provides many methods to estimate a quantile under monotonicity constraints.
}

\note{
 Inputs X.input and Y.input are useful only for method = "Bounds" 
}

\references{
Bousquet, N. (2012) Accelerated monte carlo estimation of exceedance probabilities under monotonicity constraints. Annales de la Faculte des Sciences de Toulouse. XXI(3), 557-592.
}

\author{
  Vincent Moutoussamy \cr
}

\examples{
\dontrun{
 inputDistribution <- list()
 inputDistribution[[1]] <- list("norm",c(4,1))
 inputDistribution[[2]] <- list("norm",c(0,1))

 inputDimension <- length(inputDistribution)
 dir.monot <- c(1, -1)
 N.calls <- 80

 f <- function(x){
   return(x[1] - x[2])
 }

 probability <- 1e-2

 trueQuantile <- qnorm(probability,
                     inputDistribution[[1]][[2]][1] - inputDistribution[[2]][[2]][1],
                     sqrt(inputDistribution[[1]][[2]][2] + inputDistribution[[1]][[2]][2]))

 resQuantile <- MonotonicQuantileEstimation(f, inputDimension, inputDistribution,
                                      dir.monot, N.calls, p = probability, method = "MonteCarloIS")

 quantileEstimate <- resQuantile[[1]][N.calls, 3]

}
}