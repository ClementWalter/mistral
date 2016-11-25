\name{SMART}
\alias{SMART}
\title{Support-vector Margin Algoritm for Reliability esTimation}

\description{
Calculate a failure probability with SMART method. This should not be used by itself but only through S2MART.
}

\usage{
SMART(dimension,
      lsf,
      N1 		          = 10000,
      N2 	          	= 50000,
      N3 	          	= 200000,
      Nu 	           	= 50,
      lambda1     		= 7,
      lambda2     		= 3.5,
      lambda3 	    	= 1,
      tune_cost     	= c(1,10,100,1000),
      tune_gamma     	= c(0.5,0.2,0.1,0.05,0.02,0.01),
      clusterInMargin = TRUE,
      alpha_margin 	  = 1,
      k1 	          	= round(6*(dimension/2)^(0.2)),
      k2 	          	= round(12*(dimension/2)^(0.2)),
      k3 	          	= k2 + 16,
      learn_db  	    = NULL,
      lsf_value 	    = NULL,
      failure   	    = 0,
      limit_fun_MH 	  = NULL,
      sampling_strategy = "MH",
      seeds 		      = NULL,
      seeds_eval 	    = NULL,
      burnin 		      = 30,
      thinning 		    = 4,
      plot 		        = FALSE,
      limited_plot 	  = FALSE,
      add 		        = FALSE,
      output_dir 	    = NULL,
      z_MH 		        = NULL,
      z_lsf 		      = NULL,
      verbose 		    = 0)
}

\arguments{
  \item{dimension}{an integer giving the dimension of the input space.}
  \item{lsf}{the failure fonction.}
  \item{N1}{an integer defining the number of uniform samples for (L)ocalisation stage.}
  \item{N2}{an integer defining the number of uniform samples for (S)tabilisation stage.}
  \item{N3}{an integer defining the number of gaussian standard samples for (C)onvergence stage, and so Monte-Carlo population size.}
  \item{Nu}{an integer defining the size of the first Design Of Experiment got by uniforme sampling in a sphere of radius the maximum norm of \code{N3} standard samples.}
  \item{lambda1}{a real defining the relaxing paramater in the Metropolis-Hastings algorithm for stage L.}
  \item{lambda2}{a real defining the relaxing paramater in the Metropolis-Hastings algorithm for stage S.}
  \item{lambda3}{a real defining the relaxing paramater in the Metropolis-Hastings algorithm for stage C. This shouldn't be modified as Convergence stage population is used to estimate failure probability.}
  \item{tune_cost}{a vector containing proposed values for the cost parameter of the SVM.}
  \item{tune_gamma}{a vector containing proposed values for the gamma parameter of the SVM.}
  \item{clusterInMargin}{margin points to be evaluated during refinements steps are got by mean of clustering of the \code{N1}, \code{N2} or \code{N3} points lying in the margin. Thus, they are not necessarily located into the margin. This boolean, if TRUE, enforces the selection of margin points by selecting points randomly in each cluster.}
  \item{alpha_margin}{a real value defining the margin. While 1 is the \sQuote{real} margin for a SVM, one can decide here to stretch it a bit.}
  \item{k1}{Rank of the first iteration of step S (ie stage L from 1 to \code{k1}-1).}
  \item{k2}{Rank of the first iteration of step C (ie stage S from \code{k1} to \code{k2}-1).}
  \item{k3}{Rank of the last iteration of step C (ie stage C from \code{k2} to \code{k3}).}
  \item{learn_db}{optional. A matrix of already known points, with dim : \code{dimension} x number_of_vector.}
  \item{lsf_value}{values of the limit state function on the vectors given in \code{learn_db}.}
  \item{failure}{the value defining the failure domain F = \{ x | \code{limit_state_function}(x) < \code{failure} \}.}
  \item{limit_fun_MH}{optional. If the working space is to be reduced to some subset defining by a function, eg. in case of use in a Subset Simulation algorithm. As for the \code{limit_state_function}, failure domain is defined by points whom values of \code{limit_fun_MH} are negative.}
  \item{sampling_strategy}{either "AR" or "MH", to specify which sampling strategy is to be used when generating Monte-Carlo population in a case of subset simulation : "AR" stands for \sQuote{accept-reject} while "MH" stands for Metropolis-Hastings.}
  \item{seeds}{optional. If \code{sampling_strategy}=="MH", seeds from which MH algorithm starts. This should be a matrix with \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.}
  \item{seeds_eval}{optional. The value of the \code{limit_fun_MH} on the \code{seeds}.}
  \item{burnin}{a burnin parameter for Metropolis-Hastings algorithm. This is used only for the last C step population while it is set to 0 elsewhere.}
  \item{thinning}{a thinning parameter for Metropolis-Hastings algorithm. This is used only for the last C step population while it is set to 0 elsewhere. \code{thinning} = 0 means no thinning.}
  \item{plot}{a boolean parameter specifying if function and samples should be plotted. The plot is refreshed at each iteration with the new data. Note that this option is only to be used when working on \sQuote{light} limit state functions as it requires the calculus of this function on a grid of size 161x161 (plot is done a -8:8 x -8:8 grid with 161 meshes.}
  \item{limited_plot}{only a final plot with \code{limit_state_function}, final DOE and metamodel. Should be used with \code{plot}==FALSE.}
  \item{add}{optional. "TRUE" if plots are to be added to the current active device.}
  \item{output_dir}{optional. If plots are to be saved in .jpeg in a given directory. This variable will be pasted with "_SMART.jpeg" to get the full output directory.}
  \item{z_MH}{optional. For plots, if metamodel has already been evaluated on the grid then \code{z_MH} (from outer function) can be provided to avoid extra computational time.}
  \item{z_lsf}{optional. For plots, if LSF has already been evaluated on the grid then \code{z_lsf} (from outer function) can be provided to avoid extra computational time.}
  \item{verbose}{Eiher 0 for an almost no output message, or 1 for medium size or 2 for full size}
}

\value{
  An object of class \code{list} containing the failure probability and some more outputs as described below:

  \item{proba}{The estimated failure probability.}
  \item{cov}{The coefficient of variation of the Monte-Carlo probability estimate.}
  \item{gamma}{The gamma value corresponding to the correlation between Monte-Carlo samples got from Metropolis-Hastings algorithm.}
  \item{Ncall}{The total number of calls to the \code{limit_state_function}.}
  \item{learn_db}{The final learning database, ie. all points where \code{limit_state_function} has been calculated.}
  \item{lsf_value}{The value of the \code{limit_state_function} on the learning database.}
  \item{meta_fun}{The metamodel approximation of the \code{limit_state_function}. A call output is a list containing the value and the standard deviation.}
  \item{meta_model}{The final metamodel.}
  \item{points}{Points in the failure domain according to the metamodel.}
  \item{meta_eval}{Evaluation of the metamodel on these points.}
  \item{z_meta}{If \code{plot}==TRUE, the evaluation of the metamodel on the plot grid.}
}

\details{

\code{SMART} is a reliability method proposed by J.-M. Bourinet et al. It makes uses of a SVM-based metamodel to approximate the limit state function and calculate the failure probability with a crude Monte-Carlo method using the metamodel-based limit state function. As SVM is a classification method, it makes use of limit state function values to create two classes : greater and lower than the failure threshold. Then the border is taken as a surogate of the limit state function.

Concerning the refinement strategy, it distinguishes 3 stages, known as Localisation, Stalibilsation and Convergence stages. The first one is proposed to reduce the margin as much as possible, the second one focuses on switching points while the last one works on the final Monte-Carlo population and is designed to insure a strong margin ; see F. Deheeger PhD thesis for more information.
}

\note{
Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}} to do so.
Furthermore, each time a set of vector is defined as a matrix, \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.
}

\references{
  \itemize{
    \item
      J.-M. Bourinet, F. Deheeger, M. Lemaire:\cr
      \emph{Assessing small failure probabilities by combined Subset Simulation and Support Vector Machines}\cr
      Structural Safety (2011)

    \item
      F. Deheeger:\cr
      \emph{Couplage mecano-fiabiliste : 2SMART - methodologie d'apprentissage stochastique en fiabilite}\cr
      PhD. Thesis, Universite Blaise Pascal - Clermont II, 2008
  }
}

\author{
  Clement Walter\cr
  \email{clement.walter@cea.fr}
}

\seealso{
  \code{\link{SubsetSimulation}}
  \code{\link{MonteCarlo}}
  \code{\link[e1071]{svm}} (in package \pkg{e1071})
  \code{\link{S2MART}}
}