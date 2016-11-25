estimateSURfast <- function(PPP, xi_PPP_X, krig, approx=FALSE){
  # fir foreach NOTE
  m_s <- NULL
  
  x = PPP$final_X
  u = PPP$final_U
  u_lpa <- NULL
  N <- length(PPP$final_U)
  xi_x <- xi_PPP_X

  x = cbind(PPP$X, x)
  u = c(PPP$U, u)
  names_event = as.numeric(names(u))
  u_lpa <- rep(PPP$U, each = N)

  Nsur <- length(u)

  # get updated sd
  sd_new_2 <- krig$updateSdfast(x = x, xi_x = xi_x)

  # transform vectors to matrices
  sd_old_2 <- t(matrix(1, nrow = Nsur, ncol = Nsur)*xi_x$sd^2)
  u <- t(matrix(1, nrow = Nsur, ncol = Nsur)*u)
  mean <- t(matrix(1, nrow = Nsur, ncol = Nsur)*xi_x$mean)

  # simulate U_2
  delta.sd <- 1 - sd_new_2/sd_old_2
  m_u2 <- delta.sd*(u - mean) + mean
  sd_u2 <- sqrt(sd_new_2*(1+delta.sd))

  if(approx==FALSE){
    N_batch = foreach::getDoParWorkers()
    m_s_tot = cbind(m_u2, sd_u2)
    sur <- foreach::foreach(m_s = iterators::iter(m_s_tot, by = 'row', chunksize = ceiling(Nsur/N_batch)), .combine = c, .verbose=TRUE, .export='getLPA', .packages = "iterators") %dopar% {
      sur <- sapply(1:dim(m_s)[1], function(i) {
        m = head(c(m_s[i,]), Nsur)
        s = tail(c(m_s[i,]), Nsur)
        m.LPA <- getLPA(m, names_event, Nsur-N+1)[,-1];
        sd.LPA <- getLPA(s, names_event, Nsur-N+1)[,-1];
        sum(pnorm(u_lpa, mean =c(m.LPA), sd = c(sd.LPA)))/N^2
      })
    }
  } else {
    sur <- rowMeans(sd_u2)
  }

  ind_min <- which.min(sur)
  return(list(x = x[,ind_min], u = u[ind_min], t = ind_min/N))
}
