\name{mistral-package}
\alias{mistral-package}
\alias{mistral}
\docType{package}
\title{Methods In Structural Reliability Analysis}

\description{
Provide tools for structural reliability analysis (failure probability, quantile).
}

\details{
 \tabular{lll}{
   Package: \tab mistral\cr
   Type: \tab Package\cr
   Version: \tab 2.2.1
   Date: \tab 2016-12-21\cr
   License: \tab CeCILL
  }
This package provides tools for structural reliability analysis:
\itemize{
    \item Calculate failure probability with FORM method and importance sampling.
    \item Calculate failure probability with crude Monte Carlo method
    \item Calculate failure probability with Subset Simulation algorithm
    \item Calculate failure probability with Monotonic Reliability Methods (MRM)
    \item Calculate failure probability with metamodel based algorithms : AKMCS, SMART and MetaIS
    \item Calculate failure probability with a metamodel based Subset Simulation : S2MART
    \item Wilks formula: Compute a quantile (or tolerance interval) with a given confidence level from a i.i.d. sample,
    \item Wilks formula: Compute the minimal sample size to estimate a quantile with a given confidence level,
    \item Calculate a quantile under monotonicity constraints
}
}

\author{Clement Walter, Gilles Defaux, Bertrand Iooss, Vincent Moutoussamy, with contributions from Nicolas Bousquet, Claire Cannamela and Paul Lemaitre
(maintainer: Bertrand Iooss \email{biooss@yahoo.fr})}

\references{
S.-K. Au, J. L. Beck. Estimation of small failure probabilities in high dimensions by Subset Simulation. Probabilistic Engineering Mechanics, 2001.

J.-M. Bourinet, F. Deheeger, M. Lemaire. Assessing small failure probabilities by combined Subset Simulation and Support Vector Machines. Structural Safety, 2011.

N. Bousquet. Accelerated monte carlo estimation of exceedance probabilities under monotonicity constraints. Annales de la Faculte des Sciences de Toulouse. XXI(3), 557-592, 2012

H.A. David and H.N. Nagaraja. Order statistics, Wiley, 2003.

F. Deheeger. Couplage mecano-fiabiliste : 2SMART - methodologie d'apprentissage stochastique en fiabilite. PhD. Thesis, Universite Blaise Pascal - Clermont II, 2008

A. Der Kiureghian, T. Dakessian. Multiple design points in first and second-order reliability. Structural Safety, vol.20, 1998.

O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996.

V. Dubourg. Meta-modeles adaptatifs pour l'analyse de fiabilite et l'optimisation sous containte fiabiliste. PhD. Thesis, Universite Blaise Pascal - Clermont II, 2011.

B. Echard, N. Gayton, M. Lemaire. AK-MCS : an Active learning reliability method combining Kriging and Monte Carlo Simulation

M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.

V. Moutoussamy. Contributions to structural reliability analysis: accounting for monotonicity constraints in numerical models, PhD Thesis of Universite de Toulouse, France, 2015.

W.T. Nutt and G.B. Wallis. Evaluation of nuclear safety from the outputs of computer codes in the presence of uncertainties. Reliability Engineering and System Safety, 83:57-77, 2004.

P.-H. Waarts. Structural reliability using finite element methods: an appraisal of DARS, Directional Adaptive Response Surface Sampling. PhD. Thesis, Technical University of Delft, The Netherlands, 2000.

C. Walter. Using Poisson processes for rare event simulation, PhD Thesis of Universite Paris Diderot, France, 2016.
 
S.S. Wilks. Determination of Sample Sizes for Setting Tolerance Limits. Annals Mathematical Statistics, 12:91-96, 1941.

}

\examples{
########## FORM ########### 
# u.dep is a starting point for the research of the Most Probable Failing Point
# N.calls is a total number of calls
form <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000,
                    u.dep = c(0,0))
form$p

# use IS=TRUE to use an Importance Sampling scheme with a Gaussian standard
# proposal distribution centred at the MPFP
form.IS <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000,
                       u.dep = c(0,0),
                       IS = TRUE)
form.IS$p

########### Wilks ##########

N <- WilksFormula(0.95,0.95,order=1)
print(N)
}
\keyword{ package }
