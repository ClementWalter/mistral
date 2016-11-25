# Estimate a SUR criterion
fSUR <- function(x,
                  # param x the points used to estimate the criterion
                 u,
                 # param u the corresponding realisation of u from the Poisson process
                 u_lpa,
                 # param u_lpa the sequence of the repeated events of the superposed
                 # Poisson process. For faster estimation of the integrated criterion.
                 xnew,
                 # param xnew the points added to the DoE
                 Xnew = NULL,
                 # param Xnew a matrix of points added to the Doe. In total, the criterion is
                 # estimated with cbind(Xnew, xnew) new points to the DoE. This is for
                 # sequential optimal search for multi-point SUR.
                 xi_x,
                 # param xi_x the output of xi(x), xi from \code{ok} class
                 xi_xnew = krig$xi(xnew),
                 # param xi_xnew the output of xi(xnew), xi from \code{ok} class
                 xi_Xnew,
                 # param xi_Xnew the output of xi(Xnew), xi from \code{ok} class
                 krig,
                 # param krig a kriging model from \code{ok} class
                 integrated,
                 # param integrated if the integrated criterion is to be estimated or the usual SUR.
                 q,
                 # param q the threshold of interest
                 approx = FALSE,
                 # param approx should an approximation of the criterion be computed instead. This
                 # is approximately 2x faster.
                 pnorm = stats::pnorm,
                 # param pnorm The cdf of a standard Gaussian random variable.
                 N
                 # param N the parameter of the Poisson process
){

  # if(missing(approx.pnorm)) approx.pnorm = pnorm

  # get points added to the model
  if(!is.null(Xnew)) {
    if(missing(xi_Xnew)) xi_Xnew = krig$xi(Xnew)
    xi_Xnew$mean = c(xi_Xnew$mean, xi_xnew$mean)
    xi_Xnew$sd = c(xi_Xnew$sd, xi_xnew$sd)
    xi_Xnew$kn = cbind(xi_Xnew$kn, xi_xnew$kn)
    xi_Xnew$Kinv_kn = cbind(xi_Xnew$Kinv_kn, xi_xnew$Kinv_kn)
    Xnew = cbind(Xnew, xnew)
  } else {
    Xnew <- as.matrix(xnew)
    xi_Xnew <- xi_xnew
  }

  # calculate new kriging variance
  sd.new <- krig$updateSd(x = x,
                          Xnew = Xnew,
                          xi_x = xi_x,
                          xi_Xnew = xi_Xnew)

  # simulate U_2
  delta.sd <- 1 - (sd.new/xi_x$sd)^2
  m_u2 <- delta.sd*(u - xi_x$mean) + xi_x$mean
  sd_u2 <- sqrt(sd.new^2*(1+delta.sd))
  # sur <- list(int=NA, q=NA)

  # do MC estimation
  if(integrated==TRUE){
    m.LPA <- getLPA(m_u2, as.numeric(names(u)), length(unique(names(u))), length(m_u2)-N+1)[,-1];
    sd.LPA <- getLPA(sd_u2, as.numeric(names(u)), length(unique(names(u))), length(m_u2)-N+1)[,-1];
    u_tmp <- (u_lpa - c(m.LPA))/c(sd.LPA)
    # if(approx) u_tmp <- tail(sort(u_tmp), N)
    if(approx!=FALSE){
      sur <- switch(approx,
                    s_u2 = mean(sd_u2),
                    s_norm = mean(u_tmp),
                    {u_tmp = tail(sort(u_tmp), 2*N);
                    sum(pnorm(u_tmp)/N^2)}
      )
    } else {
      sur <- sum(pnorm(u_tmp)/N^2)
    }
  } else {
    sur = mean(pnorm((q- tail(m_u2, N))/tail(sd_u2, N)))
  }

  return(sur)
}
