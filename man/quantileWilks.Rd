\name{quantileWilks}
\alias{quantileWilks}
\title{Computing quantiles with the Wilks formula}

\description{
From the Wilks formula, compute a quantile (or a tolerance interval) with a given confidence level from a i.i.d. sample, or compute the minimal sample size to estimate a quantile (or a tolerance interval) with a given confidence level.
}

\usage{
quantileWilks(alpha=0.95,beta=0.95,data=NULL,bilateral=FALSE)
}

\arguments{
  \item{alpha}{
     level of the unilateral or bilateral quantile (default = 0.95)
   }
  \item{beta}{
     level of the confidence interval on quantile value(s) (default = 0.95)
   }
  \item{data}{
     the data sample (vector format) to compute the quantile(s);
     if data=NULL (by default), the function returns the minimal sample size to compute the required quantile 
   }
  \item{bilateral}{
     TRUE for bilateral quantile (default = unilateral = FALSE)
   }
}

\value{
4 output values if 'data' is specified; 1 output value (nmin) if 'data' is not specified
  \item{lower}{
  lower bound of the bilateral tolerance interval; if bilateral=FALSE, no value
  }
  \item{upper}{
  upper bound of the tolerance interval (bilateral case) or quantile value (unilateral case)
  }
  \item{nmin}{
  minimal size of the required i.i.d. sample for given alpha and beta:
  - bilateral case: tolerance interval will be composed with the min and max of the sample;
  - unilateral case: the quantile will correspond to max of the sample.
  }
  \item{ind}{
  the index (unilateral case) or indices (bilateral case) of the quantiles in the ordered sample (increasing order) 
  }
}
  
\author{
  Claire Cannamela and Bertrand Iooss
}

\examples{ 
N <- quantileWilks(alpha=0.95,beta=0.95)
print(N)
}

\references{
H.A. David and H.N. Nagaraja. Order statistics, Wiley, 2003.

W.T. Nutt and G.B. Wallis. Evaluation of nuclear safety from the outputs of computer codes in the presence of uncertainties. Reliability Engineering and System Safety, 83:57-77, 2004.

S.S. Wilks. Determination of Sample Sizes for Setting Tolerance Limits. Annals Mathematical Statistics, 12:91-96, 1941.
}
