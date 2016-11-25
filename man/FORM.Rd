\name{FORM}
\alias{FORM}
\title{ FORM method}

\description{
  Calculate failure probability by FORM method and important sampling.
}

\usage{
FORM(dimension, lsf, u.dep, N.calls, eps = 1e-7,
     Method = "HLRF", IS = FALSE, IS.ratio = 0.5)
}

\arguments{
  \item{dimension}{
     The dimension of the input space
   }
  \item{lsf}{
     A failure fonction
   }
  \item{u.dep}{
     A vector, starting point to the research of the design point  
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
  \item{IS.ratio}{
     Ratio of N.calls for the reserach of the design point by FORM. Default = 0.5. 1-\code{IS.ratio} = the remaining ratio to use importance sampling. 
   }
}

\value{
  \item{p}{Failure probability}
  \item{indice.reliab}{Reliability index}
  \item{Ncall}{Number of calls to f}
  \item{Design.Point}{Coordinates of the design point}
  \item{fact.imp}{Importance factors}
  \item{variance}{Standard error of the probability estimator (if IS = TRUE)}
  \item{Interval.conf}{Confidence interval of the estimator at 0.95 (if IS = TRUE)}
  \item{DOE}{List which contains the design of experiments}
}

\author{
  Vincent Moutoussamy, Clement Walter
}

\details{
    This function estimate the probability that the output of the failure function is negative using FORM algorithm.
    The importance sampling procedure estimate a probability using a gaussian distribution centered in the design point with a covariance matrix equal to the indendity.
}

\references{
O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.
}

\examples{ 
\dontrun{
# u.dep is a starting point for the research of the Most Probable Failing Point
# N.calls is a total number of calls

form <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000, u.dep = c(0,0))
form$p

# use IS=TRUE to use an Importance Sampling scheme with a Gaussian standard
# proposal distribution centred at the MPFP
form.IS <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000,
                         u.dep = c(0,0),
                         IS = TRUE)
form.IS$p}
}