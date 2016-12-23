\name{FORMv0}
\alias{FORMv0}
\title{FORM method (old version)}

\description{
  Calculate failure probability by FORM method and important sampling.
}

\usage{
FORMv0(f, u.dep, inputDist, N.calls, eps = 1e-7,
     Method = "HLRF", IS = FALSE, q = 0.5, copula = "unif")
}

\arguments{
  \item{f}{
     A failure fonction
   }
  \item{u.dep}{
     A vector, starting point to the research of the design point  
   }
  \item{inputDist}{
     A list which contains the name of the input distribution and their parameters.
     For the input "i", inputDistribution[[i]] = list("name_law",c(parameters1,..., parametersN))
   }
  \item{N.calls}{
     Number of calls to f allowed
   }
  \item{eps}{
     Stop criterion : distance of two points between two iterations
   }
  \item{Method}{
     Choice of the method to research the design point: "AR" for Abdo-Rackwitz and "HLRF" for Hasofer-Lindt-Rackwitz-Fiessler
   }
  \item{IS}{
    "TRUE" for using importance Sampling method (applied after FORM which provides the importance density). Default = "FALSE".
   }
  \item{q}{
     Ratio of N.calls for the research of the design point by FORM. Default = 0.5. 1-q = the remaining ratio to use importance sampling. 
   }
  \item{copula}{
     Choice of the copula. Default = "unif" (uniform copula)
   }
}

\value{
  \item{pf}{Failure probability}
  \item{beta}{Reliability index (beta)}
  \item{compt.f}{Number of calls to f}
  \item{design.point}{Coordinates of the design point}
  \item{fact.imp}{Importance factors}
  \item{variance}{Standard error of the probability estimator (if IS = TRUE)}
  \item{conf}{Confidence interval of the estimator at 0.95 (if IS = TRUE)}
  \item{x}{A data frame containing the input design of experiments}
  \item{y}{A vector of model responses (corresponding to x)}
  \item{dy}{A data frame of model response derivatives (wrt each input and corresponding to x); for the IS sample, the derivatives are not computed}
}

\author{
Vincent Moutoussamy and Bertrand Iooss
}

\details{
    This function estimate the probability that the output of the failure function is negative using FORM algorithm. The importance sampling procedure estimate a probability using a Gaussian distribution centered in the design point with a covariance matrix equal to the indentity.
}

\references{
O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.
}

\examples{ 
\dontrun{
 distribution = list()
 distribution[[1]] = list("gamma",c(2,1))
 distribution[[2]] = list("gamma",c(3,1))
 
 f <- function(X){
    X[1]/sum(X) - qbeta((1e-5),2,3)
 }

 res <- mistral:::FORMv0(f, u.dep = c(0,0.1), inputDist = distribution, 
     N.calls = 1000, eps = 1e-7, Method = "HLRF", IS = "TRUE", 
     q = 0.1, copula = "unif")
                  
names(res)
print(res)
print(res$pf)
}
}


