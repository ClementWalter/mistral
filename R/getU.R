getU = function(xi, g_prev, p_0 = 0.1, u_high, u_low, rho = p_0*1e-2, optim=FALSE, failure){
  if(p_0=='LPA') res = min(xi$mean)
  else{
    if(u_low==-Inf) u_low = min(xi$mean - 2*xi$sd)

    if(optim==TRUE){
      u = function(u_0) (p_0 - mean(pnorm((u_0 - xi$mean)/xi$sd, lower.tail = FALSE)/g_prev))^2
      res = optim(mean(xi$mean), u, method = "BFGS")$par
    }
    else{
      u_cur = (u_high + u_low)/2
      p_cur = mean(pnorm((u_cur - xi$mean)/xi$sd, lower.tail = FALSE)/g_prev)
      while(abs(p_cur-p_0)>rho){
        if(p_cur>p_0) u_low = u_cur
        else{u_high = u_cur}
        u_cur = (u_high + u_low)/2
        p_cur = mean(pnorm((u_cur - xi$mean)/xi$sd, lower.tail = FALSE)/g_prev)
      }
      res = u_cur
    }
  }
  return(min(res, failure))
}
