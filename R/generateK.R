#' Generate Standard Gaussian samples with a Gaussian transiiton kernel
#' 
#' @details This function generates standard Gaussian samples with a Markov Chain
#' using a suitable transition kernel
#'
#' @author Clement WALTER \email{clementwalter@icloud.com}
#'
#' @examples
#' # Get a seed in dimension 2
#' X <- matrix(rnorm(2), nrow = 2)
#' X <- generateK(X, N = 1000)
#' 
#' library(ggplot2)
#' ggplot(as.data.frame(t(X)), aes(x_1,x_2)) + geom_point()
#' 
#' # One can also specify a limit-state function
#' lsf <- function(X){
#'      sqrt(colSums(X^2)) > 2
#' }
#' X <- matrix(c(2, 2), nrow = 2)
#' X <- generateK(X, N = 1000, lsf = lsf)
#' 
#' ggplot(as.data.frame(t(X)), aes(x_1,x_2)) + geom_point()
#' 
#' @export
generateK= function(X,
                    #' @param X the seeds for the Markov Chain. There are as many MC drawn as given
                    #' seeds
                    N = 100,
                    #' @param N the number of desired samples"'
                    thinning = 4,
                    #' @param thinning the proportion of kept samples, ie. 1 each \code{thinning} draw.
                    sigma = 1,
                    #' @param sigma the exploration parameter for the transition kernel
                    lsf,
                    #' @param lsf a boolean limit-state function for definig a subdomain of the input
                    #' space.
                    burnin=20
                    #' @param burnin the \code{burnin} parameter, ie. the number of discarded samples
                    #' before keeping one.
) {
  
  tic = proc.time()[3]
  
  # select only non duplicated seeds
  X = as.matrix(X)
  sel <- !duplicated(t(X))
  X <- as.matrix(X[,sel])
  n_seeds = ncol(X)
  
  # Define transition kernel
  K <- function(X, sigma){
    W <- matrix(rnorm(X), nrow = nrow(X))
    X <- (X + sigma*W)/sqrt(1 + sigma^2)
    X
  }
  
  # set variables
  # N is the number of desired samples
  # N_chain <- ceiling(N/n_seeds) is the number of samples per chain
  # chain_lenght <- burnin + thinning*N_chain is the number of iterations of the kernel
  N_chain <- ceiling(N/n_seeds)
  chain_length <- burnin + thinning*N_chain
  
  Ncall = 0;
  
  if(missing(lsf)) {
    lsf = function(x) {TRUE}
  }
  
  Xlist <- foreach(icount(chain_length)) %do% {
    Xstar <- K(X, sigma)
    ystar <- lsf(Xstar)
    Xstar[!ystar] <- X[!ystar]
    X <- Xstar
    as.matrix(X)
  }
  
  sel <- burnin + c(1:N_chain)*thinning
  X <- do.call(cbind, Xlist[sel])[,1:N]
  rownames(X) <- rep(c('x_1', 'x_2'), length.out = nrow(X))
  
  toc = proc.time()[3]-tic
  cat("   -",N,"points generated in",toc,"sec.; ",sum(duplicated(t(X))), "duplicated samples \n")

  #' @return A matrix \code{X} with the number of desired samples 
  return(X)
}