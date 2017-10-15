getPPP <- function(dimension,
                   N = 1000,
                   N.batch = foreach::getDoParWorkers(),
                   burnin = 20,
                   q,
                   p_max = 1e-20,
                   xi,
                   fast=TRUE,
                   Niter_max=function(N) ceiling(-N*log(p_max)),
                   keep_process = FALSE,
                   sorted = TRUE,
                   verbose = 0){

  # Fix foreach NOTE
  x_aug <- NULL
  
  # Particles initilisation
  X <- matrix(rnorm(dimension*N), dimension, dimnames = list(rep(c('x', 'y'), l = dimension)))
  xi_X = xi(X)
  U = rnorm(N, mean = xi_X$mean, sd = xi_X$sd)
  names(U) = 1:N
  X_aug = rbind(X,U)
  M = 0
  Niter = 0

  # Transition Kernel in standard Gaussian input space for MH
  # sigma_hist <- sigma <- krig@covariance@range.val/2
  sigma_hist <- sigma <- 0.3
  K = function(x){
    x = as.matrix(x);dimnames(x) <- NULL;
    w = matrix(rnorm(length(x)), nrow = dimension)
    x_star = (x + sigma*w)/(sqrt(1+sigma^2))
  }

  PPP <- foreach::foreach(
    x_aug = iterators::iter(X_aug, by ='col', chunksize = ceiling(N/N.batch)),
    .inorder = FALSE,
    .options.multicore = list(set.seed = TRUE),
    .packages = 'foreach',
    .verbose = (verbose>1)) %dopar% {
      X = head(x_aug, -1)
      U_tmp <- U <- c(tail(x_aug, 1))
      N = length(U)
      # Niter_max = ceiling(-N*log(p_max))
      names(U) = colnames(x_aug)
      x_min = which.min(U)

      # Sum Poisson process
      if(keep_process==TRUE){
        PPP_X <- X
        PPP_U <- U
      } else{
        PPP_X = NULL
        PPP_U = NULL
      }
      while((U[x_min]<q) && (Niter < Niter_max(N))){
        if(verbose>1){
          cat(" * level =",U[x_min],"&", Niter, "<", Niter_max(N), "\n")
          cat(' # Generate "big" population \n')
        }
        X_pop <- X
        U_pop <- U
        if(fast==FALSE){
          ind_star <- sample(which(U>U[x_min]), 1)
          X_pop = as.matrix(X_pop[,ind_star])
          U_pop = U_pop[ind_star]
        } else {
          ind_star = 1:length(U_pop)
        }
        acceptance = 0
        for(iter in 1:burnin) {
          X_star = K(X_pop)
          xi_star <- xi(X_star)
          U_star = rnorm(length(U_pop), mean = xi_star$mean, sd = xi_star$sd)
          accept <- U_star > U[x_min]
          X_pop[,accept] <- X_star[,accept]
          U_pop[accept] = U_star[accept]
          acceptance = acceptance + sum(accept)
        }
        X_pop = matrix(X_pop[,U_pop!=U[ind_star]], nrow = nrow(X_pop))
        U_pop = U_pop[U_pop!=U[ind_star]]
        #         acceptance = acceptance/burnin/N
        #         if(acceptance<0.2) sigma = sigma*0.9
        #         if(acceptance>0.5) sigma = sigma*1.1
        #         sigma_hist = c(sigma_hist, sigma)

        # Remove possible U[x_min] values
        keep <- U_pop>U[x_min]

        if(verbose>1) cat(' # Make move with the discretised distribution\n')
        Nmove_dis = M
        while(sum(keep)>0 && (U[x_min]<q) && (M-Nmove_dis<dimension)) {
          # get conditional population
          X_pop = matrix(X_pop[,keep], nrow = nrow(X_pop))
          U_pop = U_pop[keep]

          # sample new event with the discretised distribution
          ind = sample(1:length(U_pop), size = 1)
          X[,x_min] <- as.matrix(X_pop[,ind])
          U[x_min] <- U_pop[ind]

          # save state
          if(keep_process){
            PPP_X <- cbind(PPP_X, X_pop[,ind])
            PPP_U <- c(PPP_U, U_pop[ind])
          }
          M = M + 1
          U_pop <- U_pop[-ind]
          X_pop <- matrix(X_pop[,-ind], nrow = nrow(X_pop))

          # prepare next step
          x_min = which.min(U)
          keep <- U_pop>U[x_min]
        }
        Niter = Niter + 1
        Nmove_dis = M - Nmove_dis
      }
      if(Niter==Niter_max(N)) warning("In PPP generation: maximum number of iterations reached")
      ind_ppp <- PPP_U<q
      list(PPP_X = matrix(PPP_X[,ind_ppp], nrow = nrow(PPP_X)), PPP_U = PPP_U[ind_ppp], final_X = X, final_U = U, M = M, sigma = sigma_hist)
    }
  M = sum(sapply(PPP, function(l) {l$M}))
  sigma = sapply(PPP, function(l) {l$sigma})
  PPP_X <- do.call(cbind, lapply(PPP, function(l) {l$PPP_X}))
  PPP_U <- do.call(c, lapply(PPP, function(l) {l$PPP_U}))
  colnames(PPP_X) = names(PPP_U)
  final_X <- do.call(cbind, lapply(PPP, function(l) {l$final_X}))
  final_U <- do.call(c, lapply(PPP, function(l) {l$final_U}))

  if(sorted==TRUE){
    PPP_U <- sort(PPP_U, index.return = TRUE)
    PPP_X <- matrix(PPP_X[,PPP_U$ix], nrow = nrow(PPP_X))
    PPP_U <- PPP_U$x
    final_U <- sort(final_U, index.return = TRUE)
    final_X <- matrix(final_X[,final_U$ix], nrow = nrow(final_X))
    final_U <- final_U$x
  }

  return(list(X = PPP_X, U = PPP_U, final_X = final_X, final_U = final_U, M = M, sigma = sigma))
}
