#' @title The modified Metropolis-Hastings algorithm
#'
#' @description The function implements the specific modified Metropolis-Hastings algorithm
#' as described first by Au \& Beck and including another scaling parameter for an extended
#' search in initial steps of the \code{SMART} algorithm.
#' 
#' @details The modified Metropolis-Hastings algorithm is supposed to be used in the
#' Gaussian standard space. Instead of using a proposed point for the multidimensional
#' Gaussian random variable, it applies a Metropolis step to each coordinate. Then it
#' generates the multivariate candidate by checking if it lies in the right domain.
#' 
#' This version proposed by Bourinet et al. includes an scaling parameter \code{lambda}.
#' This parameter is multiplied with the likelihood ratio in order to increase the chance
#' of accepting the candidate. While it biases the output distribution of the Markov chain,
#' the authors of \code{SMART} suggest its use (\code{lambda > 1}) for the exploration phase.
#' Note such a value disable to possiblity to use the output population for Monte Carlo
#' estimation.
#'
#' @export
MetropolisHastings = function(x0,
                              #' @param x0 the starting point of the Markov chain
                              eval_x0=-1,
                              #' @param eval_x0 the value of the limit-state function on \code{x0}
                              chain_length,
                              #' @param chain_length the length of the Markov chain. At the end the chain
                              #' will be chain_length + 1 long
                              modified = TRUE,
                              #' @param  modified a boolean to use either the original Metropolis-Hastings
                              #' transition kernel or the coordinate-wise one
                              sigma = 0.3,
                              #' @param sigma a radius parameter for the Gaussian or Uniform proposal
                              proposal = "Uniform",
                              #' @param proposal either "Uniform" for a Uniform random variable in an interval
                              #' [-sigma, sigma] or "Gaussian" for a centred Gaussian random variable with
                              #' standard deviation sigma
                              lambda = 1,
                              #' @param lambda the coefficient to increase the likelihood ratio
                              limit_fun = function(x) {-1},
                              #' @param limit_fun the limite-state function delimiting the domain to sample in
                              burnin=20,
                              #' @param burnin a burnin parameter, ie a number of initial discards samples
                              thinning=4
                              #' @param thinning a thinning parameter, ie that one sample over \code{thinning}
                              #' samples is kept along the chain
) {
  
  # Set initial parameters
  niter = (thinning+1)*chain_length + burnin
  d = length(x0)
  U = matrix(nrow=d,ncol=(niter+1))
  U[,1] = x0
  res = list(points=U,
             eval=eval_x0,
             acceptation=NA,
             Ncall=NULL,
             samples=matrix(nrow=d),
             eval_samples=NA);
  tau = 0
  Ncall = 0;
  rand <- switch(proposal,
                 "Uniform" = function(x){x + runif(length(x), min = -sigma, max = sigma)},
                 "Gaussian" = function(x){x + rnorm(length(x), mean = 0, sd = sigma)}
  )
  
  for (i in 1:niter) {
    candidate = rand(res$points[,i])
    
    if(modified) {
      acceptance = pmin(lambda*stats::dnorm(candidate)/stats::dnorm(res$points[,i]),1)
    } else {
      acceptance = min(lambda*prod(stats::dnorm(candidate)/stats::dnorm(res$points[,i])),1)
    }
    test = rep(runif(acceptance) < acceptance, length.out = d)
    
    candidate[!test] <- res$points[!test,i]
    
    if (prod(test)==0){
      candidate = res$points[,i];
      eval = res$eval[i];
      tau = tau + 1
    }
    else {
      eval = limit_fun(candidate); Ncall = Ncall+1;
      res$samples = cbind(res$samples,candidate); res$eval_samples = c(res$eval_samples,eval);
      indicatrice = (1-sign(eval))/2 #1 if limit_fun(x) < 0, 0 otherwise
      if (indicatrice==0 | is.nan(indicatrice)) {
        candidate = res$points[,i];
        eval = res$eval[i]
        tau = tau+1
      }
    }
    res$points[,i+1] = candidate;
    res$eval[i+1] = eval
  }
  
  tau = tau/niter
  res$acceptation = 1-tau
  sel_samples = burnin+1+(thinning+1)*c(0:chain_length)
  res$points = res$points[,sel_samples]
  res$eval = res$eval[sel_samples]
  res$Ncall = Ncall
  
  return(res)
  #' @return A list containing the following entries:
  #' \item{points}{the generated Markov chain}
  #' \item{eval}{the value of the limit-state function on the generated samples}
  #' \item{acceptation}{the acceptation rate}
  #' \item{Ncall}{the total number of call to the limit-state function}
  #' \item{samples}{all the generated samples}
  #' \item{eval_samples}{the evaluation of the limit-state function on the
  #' \code{samples} samples}
}