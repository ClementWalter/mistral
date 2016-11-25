\name{ComputeDistributionParameter}
\alias{ComputeDistributionParameter}
\title{Compute internal parameters and moments for univariate distribution functions}

\description{
  Compute the internal parameters needed in the definition of several distribution functions when unknown
}

\usage{
ComputeDistributionParameter(margin)
}

\arguments{
  \item{margin}{ A list containing the definition of the marginal distribution function}
}
\value{
  \item{margin}{The updated list}
}

\author{ gilles DEFAUX, \email{gilles.defaux@cea.fr} }

\examples{ 
distX1 <- list(type='Lnorm', MEAN=120.0, STD=12.0, P1=NULL, P2=NULL, NAME='X1')
distX1 <- ComputeDistributionParameter(distX1)
print(distX1)
}