\name{UtoX}
\alias{UtoX}
\title{Iso-probabilistic transformation from U space to X space}

\description{
  \code{UtoX} performs as iso-probabilistic transformation from standardized space (U) to physical space (X) according to the NATAF transformation, which requires only to know the means, the standard deviations, the correlation matrix \eqn{\rho(Xi,Xj) = \rho_{ij}} and the marginal distributions of Xi.
  In standard space, all random variables are uncorrelated standard normal distributed variables whereas they are correlated and defined using the following distribution functions: Normal (or Gaussian), Lognormal, Uniform, Gumbel, Weibull and Gamma.
}

\usage{ UtoX(U, input.margin, L0) }

\arguments{
  \item{U}{a matrix containing the realisation of all random variables in U-space}
  \item{input.margin}{ A list containing one or more list defining the marginal distribution functions of all random variables to be used }
  \item{L0}{the lower matrix of the Cholesky decomposition of correlation matrix R0 (result of \code{\link{ModifCorrMatrix}}) }
}

\value{ 
  \item{X}{a matrix containing the realisation of all random variables in X-space}
}

\details{
  Supported distributions are :
  \itemize{
  \item{NORMAL: }{distribution, defined by its mean and standard deviation
  \deqn{distX <- list(type="Norm", MEAN=0.0, STD=1.0, NAME="X1")}
}
  \item{LOGNORMAL: }{distribution, defined by its internal parameters P1=meanlog and P2=sdlog (\code{\link{plnorm}})
  \deqn{distX <- list(type="Lnorm", P1=10.0, P2=2.0, NAME="X2")}
}
  \item{UNIFORM: }{distribution, defined by its internal parameters P1=min and P2=max (\code{\link{punif}})
  \deqn{distX <- list(type="Unif",P1=2.0, P2=6.0, NAME="X3")}
}
  \item{GUMBEL: }{distribution, defined by its internal parameters P1 and P2
  \deqn{distX <- list(type='Gumbel',P1=6.0, P2=2.0, NAME='X4')}
}
  \item{WEIBULL: }{distribution, defined by its internal parameters P1=shape and P2=scale (\code{\link{pweibull}})
  \deqn{distX <- list(type='Weibull', P1=NULL, P2=NULL, NAME='X5')}
}
  \item{GAMMA: }{distribution, defined by its internal parameters P1=shape and P2=scale (\code{\link{pgamma}})
  \deqn{distX <- list(type='Gamma', P1=6.0, P2=6.0, NAME='X6')}
}
  \item{BETA: }{distribution, defined by its internal parameters P1=shape1 and P2=shapze2 (\code{\link{pbeta}})
  \deqn{distX <- list(type='Beta', P1=6.0, P2=6.0, NAME='X7')}
}
}
}

\seealso{
  \code{\link{ModifCorrMatrix}}, \code{\link{ComputeDistributionParameter}}
}

\author{ gilles DEFAUX, \email{gilles.defaux@cea.fr} }

\examples{ 
Dim = 2

distX1 <- list(type='Norm',  MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X1')
distX2 <- list(type='Norm',  MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X2')

input.margin <- list(distX1,distX2)
input.Rho    <- matrix( c(1.0, 0.5,
                          0.5, 1.0),nrow=Dim)
input.R0     <- ModifCorrMatrix(input.Rho)
L0           <- t(chol(input.R0))

lsf = function(U) {   
    X <- UtoX(U, input.margin, L0)
    G <- 5.0 - 0.2*(X[1,]-X[2,])^2.0 - (X[1,]+X[2,])/sqrt(2.0)
    return(G)
}

u0 <- as.matrix(c(1.0,-0.5))
lsf(u0)
}

\references{
  \itemize{
    \item M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009
    \item V. Dubourg, Meta-modeles adaptatifs pour l'analyse de fiabilite et l'optimisation sous containte fiabiliste, PhD Thesis, Universite Blaise Pascal - Clermont II,2011
  }
}
