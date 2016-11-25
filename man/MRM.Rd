\name{MRM}
\alias{MRM}
\title{MRM method}
\description{
  Estimate a failure probability by MRM method.
}
\usage{
MRM(f, inputDimension, inputDistribution, dir.monot, N.calls, Method, silent = FALSE)
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
  \item{Method}{
    there is two methods available. "MC" is an adapation of the Monte Carlo method under constraints of monotony. "MRM" is based on a sequential sampling.
  }
  \item{silent}{
    if silent = TRUE, print curent number of call to f. Default: FALSE.
  }
}

\details{
    These methods compute the probability that the output of the failure function is negative
}

\value{
    \item{Um}{Exact lower bounds of the failure probability }
    \item{UM}{Exact upper bounds of the failure probability }
    \item{MLE}{Maximum likelihood estimator of the failure probability }
    \item{IC.inf}{Lower bound of the confidence interval of the failure probability based on MLE }
    \item{IC.sup}{Upper bound of the confidence interval of the failure probability based on MLE }
    \item{CV.MLE}{Coefficient of variation of the MLE }
    \item{X}{design of experiments }
    \item{Y}{value of f on X }
    \item{N.tot}{Total number of simulation (only for "MC_monotone") }
}

\author{
  Vincent Moutoussamy and Nicolas Bousquet
}

\references{
    Bousquet, N. (2012) Accelerated monte carlo estimation of exceedance probabilities under monotonicity constraints. Annales de la Faculte des Sciences de Toulouse. XXI(3), 557-592.
}


\examples{
\dontrun{

  inputDistribution <- list()
  inputDistribution[[1]] <- list("norm",c(4,1))
  inputDistribution[[2]] <- list("norm",c(0,1))
  inputDistribution[[3]] <- list("norm",c(-1,3))

  inputDimension <- length(inputDistribution)

  p <- 1e-5

  threshold <- qnorm(p, 3, sqrt(11))

  f <- function(Input){
    sum(Input) - threshold
  }

  dir.monot <- c(1, 1, 1)

  N.calls <- 300
 
  res.MRM <- MRM(f, inputDimension, inputDistribution,
                  dir.monot, N.calls, Method = "MRM", silent = FALSE)

  N <- 1:dim(res.MRM[[1]])[1]
  
  plot(N, res.MRM[[1]][, 1],
        col = "black", lwd=2, type='l', ylim=c(0, 50*p),
        xlab="Number of runs to the failure function",
        ylab="")
  lines(N, res.MRM[[1]][, 2], col = "black", lwd = 2)
  lines(N, res.MRM[[1]][, 3], col = "red", lwd = 2)
  lines(N, res.MRM[[1]][, 7], col = "blue", lwd = 2, lty = 2)
  lines(N, rep(p, length(N)), lwd= 2, col= "orange", lty=3 )
  legend("topright",
          c("Exact Bounds", "MLE","p.hat", "p"), 
          col = c("black", "red", "blue", "orange"),
          text.col = c("black", "red", "blue", "orange"),
          lty = c(1, 1, 2, 3),
          merge = TRUE)

}

}


