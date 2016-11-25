\name{ModifCorrMatrix}
\alias{ModifCorrMatrix}
\title{Modification of a correlation matrix to use in UtoX}

\description{
  \code{ModifCorrMatrix} modifies a correlation matrix originally defined using SPEARMAN correlation coefficients to the correlation matrix to be used in the NATAF transformation performed in \code{UtoX}.
}

\usage{ ModifCorrMatrix(Rs) }

\arguments{
  \item{Rs}{ Original correlation matrix defined using SPEARMAN correlation coefficient : \deqn{R_s=[\rho_{ij}^s]}{Rs=[\rho ij^s]} }
}

\value{
  \item{R0}{Modified correlation matrix}
}

\note{
  The NATAF distribution is reviewed from the (normal) copula viewpoint as a particular and convenient means to describe a joint probabilistic model assuming that the normal copula fits to the description of the input X.
  The normal copula is defined by a symmetric positive definite matrix R0. Even though the off-diagonal terms in this matrix are comprised in ]-1; 1[ and its diagonal terms are equal to 1, it shall not be confused with the more usual correlation matrix.
  Lebrun and Dutfoy point out that the SPEARMAN (or rank) correlation coefficient is better suited to parametrize a copula because it leads to a simpler closed-form expression for \eqn{\rho_{ij}}.
  \deqn{}{ \rho ij^0 = 2 sin(\rho ij^s pi/6)}
}

\seealso{ \code{\link{UtoX}} }

\author{ Gilles DEFAUX, \email{gilles.defaux@cea.fr} }


\examples{ 
Dim <- 2
input.Rho    <- matrix( c(1.0, 0.5,
                          0.5, 1.0),nrow=Dim)
input.R0     <- ModifCorrMatrix(input.Rho)
print(input.R0)
}

\references{
  \itemize{
    \item M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009
    \item Lebrun, R. and A. Dutfoy. A generalization of the Nataf transformation to distributions with elliptical copula. Prob. Eng. Mech., 24(2), 172-178.
    \item
      V. Dubourg, Meta-modeles adaptatifs pour l'analyse de fiabilite et l'optimisation sous containte fiabiliste, PhD Thesis, Universite Blaise Pascal - Clermont II,2011
  }
}
