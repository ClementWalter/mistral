\name{WilksFormula}
\alias{WilksFormula}
\title{Sample size by Wilks formula}

\description{
  Compute Wilks formula for setting size of a i.i.d. sample for quantile estimation with confidence level or for tolerance intervals
}

\usage{
WilksFormula(alpha=0.95,beta=0.95,bilateral=FALSE,order=1)
}

\arguments{
  \item{alpha}{
     order of the quantile (default = 0.95)
   }
  \item{beta}{
     level of the confidence interval (default = 0.95)
   }
  \item{bilateral}{
     TRUE for bilateral quantile (default = unilateral = FALSE)
   }
  \item{order}{
     order of the Wilks formula (default = 1)
   }
}

\value{
  \item{N}{The minimal sample size to apply Wilks formula}
}

\author{
  Paul Lemaitre and Bertrand Iooss
}

\examples{ 
N <- WilksFormula(0.95,0.95,order=1)
print(N)
}

\references{
H.A. David and H.N. Nagaraja. Order statistics, Wiley, 2003.

W.T. Nutt and G.B. Wallis. Evaluation of nuclear safety from the outputs of computer codes in the presence of uncertainties. Reliability Engineering and System Safety, 83:57-77, 2004.

S.S. Wilks. Determination of Sample Sizes for Setting Tolerance Limits. Annals Mathematical Statistics, 12:91-96, 1941.
}
