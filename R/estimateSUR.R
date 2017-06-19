#' EstimateSUR
#'
#' A function for estimating a SUR criterion with a realisation of a PPP
# @importFrom Rcpp sourceCpp evalCpp
# @useDynLib BMP
estimateSUR <- function(PPP,
                        #' @param PPP the Poisson point process generated to get alpha.
                        xi_PPP_X,
                        #' @param xi_PPP_X the output of xi(cbind(PPP$X, PPP$final_X)).
                        integrated = TRUE,
                        #' @param integrated boolean to specify of SUR criterion is standard or
                        #' integrated.
                        N_ppp,
                        #' @param N_ppp the number of Poisson processes used for the SUR criterion estimation.
                        method = "discrete",
                        #' @param method eiter "genoud" for an optimisation using the package \code{rgenoud} or
                        #' "discrete" for a discrete search over \code{SUR_pop}.
                        SUR_pop,
                        #' @param SUR_pop if \code{optimcontrol$method=="discrete"}, SUR_pop is the population
                        #' onto which minimizer is sought. Should
                        #' be a matrix d x n.
                        r = N.batch,
                        #' @param r number of points to be added to the DoE.
                        optimcontrol = list(pop.size = 50*d,
                                            max.generations = 10*d),
                        #' @param optimcontrol a list of control parameters for the optimisation
                        #' of the SUR criterion using the \code{rgenoud} package.
                        approx.pnorm,
                        #' @param approx.pnorm (optional) an approximation of base pnorm function running faster.
                        J = 0,
                        #' @param J the center of an interval of size 8 for pnorm approximation.
                        N.batch = foreach::getDoParWorkers(),
                        #' @param N.batch Number of batchs for parallel computation.
                        verbose = 0,
                        #' @param verbose to control the print level of the algorithm
                        ...
                        #' @param ... further arguments to be passed to fSUR.
){

  # fix foreach NOTE
  X <- NULL
  
  N <- length(PPP$final_U)
  d <- dim(PPP$final_X)[1]
  n <- dim(xi_PPP_X$Kinv_kn)[1]

  if(method=='discrete' && missing(SUR_pop)){
    SUR_pop = cbind(PPP$X, PPP$final_X)
    u_SUR_pop <- c(PPP$U, PPP$final_U)
    N_pop = dim(SUR_pop)[2]
    xi_SUR_POP <- xi_PPP_X
    SUR_aug <- rbind(SUR_pop, matrix(xi_SUR_POP$mean, nrow = 1), matrix(xi_SUR_POP$sd, nrow = 1), xi_SUR_POP$kn, xi_SUR_POP$Kinv_kn)
  }

  if(!missing(N_ppp) && N_ppp<N) {
    ind <- sample(x = 1:N, size = N_ppp, replace = FALSE)
    sel1 <- names(PPP$U)%in%ind
    PPP$X <- PPP$X[,sel1]
    PPP$U <- PPP$U[sel1]
    sel2 <- names(PPP$final_U)%in%ind
    PPP$final_X <- PPP$final_X[,sel2]
    PPP$final_U <- PPP$final_U[sel2]
    xi_PPP_X$mean <- xi_PPP_X$mean[c(sel1, sel2)]
    xi_PPP_X$sd <- xi_PPP_X$sd[c(sel1, sel2)]
    xi_PPP_X$kn <- as.matrix(xi_PPP_X$kn[,c(sel1, sel2)])
    xi_PPP_X$Kinv_kn <- as.matrix(xi_PPP_X$Kinv_kn[,c(sel1, sel2)])
    N <- N_ppp
  }

  x = PPP$final_X
  u = PPP$final_U
  u_lpa <- NULL
  xi_x <- xi_PPP_X

  if(integrated == TRUE){
    x = cbind(PPP$X, x)
    u = c(PPP$U, u)
    u_lpa <- rep(PPP$U, each = N)
  } else {
    xi_x$mean = tail(xi_PPP_X$mean, N)
    xi_x$sd = tail(xi_PPP_X$sd, N)
    xi_x$kn = t(tail(t(xi_PPP_X$kn), N))
    xi_x$Kinv_kn = t(tail(t(xi_PPP_X$Kinv_kn), N))
  }


  args = list(...)

  args$x = x
  args$u = u
  args$u_lpa = u_lpa
  args$xi_x = xi_x
  args$integrated = integrated
  args$N = N

  if(approx.pnorm) {
    x_pnorm <- seq(from = -4, to = 4, by = 0.2)
    args$pnorm = stats::approxfun(x_pnorm, pnorm(x_pnorm), yleft = 0, yright = 1)
  }

  args$Xnew <- Xnew <- NULL
  args$xi_Xnew <- xi_Xnew <- NULL

  cat(" * Evaluation of SUR criterion: integrated = ",integrated, ", r = ", ifelse(is.null(args$r), 1, args$r), ", approx = ",ifelse(is.null(args$approx), FALSE, args$approx),", approx.pnorm = ",approx.pnorm,", optim = ", method,", N_ppp = ", ifelse(missing(N_ppp),N, N_ppp)," \n", sep = "")
  sur = list(x = NULL, u = NULL, t = NULL)
  for(k in 1:r){
    if(method=="discrete"){
      SUR <- foreach::foreach(X = iterators::iter(SUR_aug, by ='col', chunksize = ceiling(N_pop/N.batch)),
                              .combine = 'c',
                              .export = 'fSUR',
                              .options.multicore = list(set.seed = TRUE),
                              .errorhandling = "pass") %dopar% {
                                tmp <- c(apply(X, 2, function(xaug){
                                  args$xnew <- xaug[1:d]
                                  args$xi_xnew <- list(mean = xaug[d+1], sd = xaug[d+2], kn = xaug[(d+3):(d+2+n)], Kinv_kn = tail(xaug, n))
                                  do.call(fSUR, args)
                                }))
                                return(tmp)
                              }
      if(class(SUR)!="numeric"){
        message(' ! memory issue with parallel computing, approximated SUR used insteead !')
        args$approx = TRUE
        x_pnorm <- seq(from = -4, to = 4, by = 0.2)
        args$pnorm = stats::approxfun(x_pnorm, pnorm(x_pnorm), yleft = 0, yright = 1)
        SUR <- foreach::foreach(X = iterators::iter(SUR_aug, by ='col', chunksize = ceiling(N_pop/N.batch)),
                                .combine = 'c',
                                .export = 'fSUR',
                                .options.multicore = list(set.seed = TRUE),
                                .errorhandling = "pass") %dopar% {
                                  tmp <- c(apply(X, 2, function(xaug){
                                    args$xnew <- xaug[1:d]
                                    args$xi_xnew <- list(mean = xaug[d+1], sd = xaug[d+2], kn = xaug[(d+3):(d+2+n)], Kinv_kn = tail(xaug, n))
                                    do.call(fSUR, args)
                                  }))
                                  return(tmp)
                                }
      }
      if(class(SUR)!="numeric"){
        message(' ! memory issue with parallel computing, standard SUR used insteead !')
        args$integrated = FALSE
        SUR <- foreach::foreach(X = iterators::iter(SUR_aug, by ='col', chunksize = ceiling(N_pop/N.batch)),
                                .combine = 'c',
                                .export = 'fSUR',
                                .options.multicore = list(set.seed = TRUE),
                                .errorhandling = "pass") %dopar% {
                                  tmp <- c(apply(X, 2, function(xaug){
                                    args$xnew <- xaug[1:d]
                                    args$xi_xnew <- list(mean = xaug[d+1], sd = xaug[d+2], kn = xaug[(d+3):(d+2+n)], Kinv_kn = tail(xaug, n))
                                    do.call(fSUR, args)
                                  }))
                                  return(tmp)
                                }
      }
      if(class(SUR)!="numeric"){
        print(SUR)
      }
      ind_min = which.min(SUR)
      SUR_point <- SUR_pop[,ind_min]
      SUR_aug <- SUR_aug[,-ind_min]
      sur$u = c(sur$u, u_SUR_pop[ind_min])
      sur$t = c(sur$t, ind_min/length(SUR))

    } else {
      f_gen <- function(xnew) {
        args$xnew = xnew
        do.call(fSUR, args)
      }
      SUR_point <- do.call(rgenoud::genoud, c(optimcontrol, list(fn = f_gen, nvars = d, print.level = verbose)))$par
    }
    Xnew = cbind(Xnew, as.matrix(SUR_point))
    sur$x = cbind(sur$x, SUR_point)
  }
  return(sur)
  #' @return a list containing the points minimising the criterion
}

t.list <- function(l){
  lapply(split(do.call("c", l), names(l[[1]])), unname)
}
