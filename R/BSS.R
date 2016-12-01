#' @import ggplot2
#' @import mvtnorm
BSS = function(d, # dimension
               lsf,
               failure = 0,
               p_0 = 0.1, # p_0 subset or "LPA" for LPA
               N = 1e3, # MC population
               n_0 = 5*d, # size first DoE
               N1 = N, # size uniform pop for DoE !maximin
               maximin = FALSE, # initial DoE
               critere = "SUR", # critere pour les paliers
               optim.method = "BFGS",
               kernel = "matern5_2", # noyau de krigeage
               eta = 1e-6, # critere arret enrichissement
               K, # noyau de transition pour resampling
               burnin = 20 # burnin pour le resampling (cout 0)
){
  
  # Fix NOTE issue with R CMD check
  x <- y <- z <- ..level.. <- NULL

  ## Initialisation
  # Generate MC sample
  Y = matrix(rnorm(d*N), ncol = N, dimnames = list(c("x", "y")))
  Ru = max(sqrt(colSums(Y^2)))
  Ncall = 0
  xplot <- yplot <- c(-50:50)/10
  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
  p <- ggplot(data = df_plot, aes(x,y))
  m = 1
  lsf_unif <- function(u){
    if(is.null(dim(u))) u <- t(as.matrix(u))
    x <- t(qnorm(u))
    lsf(x)
  }
  
  if(plot==TRUE & d>2){
    message("Cannot plot in dimension > 2")
    plot <- FALSE
  }

  if(missing(K)){
    K = function(x,g){
      for(burn in 1:burnin){
        x_star = x + matrix(rnorm(x, mean = 0, sd = 1), nrow = d)
        rho = g(x_star)*dmvnorm(t(x_star))/g(x)/dmvnorm(t(x))
        accept = runif(rho)<rho
        x[,accept] = x_star[,accept]
      }
      return(x)
    }
  }

  # Initial DoE (maximim)
  ## maximin
  if(maximin){

  }
  ## clustering dans boule uniforme cf MetaIS
  else{
    tmp = runifSphere(d,N1,radius=Ru)
    learn_db = t(kmeans(tmp, centers=n_0,iter.max=20)$centers)
    row.names(learn_db) = c("x", "y")
    rm(tmp)
    G = lsf(learn_db)
    Ncall = Ncall + n_0
  }

  g = list(1)
  pi = list(function(x) 1)
  u = -Inf
  t = 1

  # Estimate first metamodel
  capture.output(meta <- trainModel(design=learn_db,
                                    response=G,
                                    kernel=kernel,
                                    optim.method = optim.method))
  xi = meta$fun(pnorm(Y))
  g_prev = 1

  while(u[t] < failure){
    cat(" * t =", t, '; u =', tail(u, 1), "\n")

    # Choose threshold
    u = c(u, getU(xi, g_prev, u_high = max(c(failure,xi$mean + 2*xi$sd)), u_low = u[t], failure = failure, p_0=p_0))
    g = c(g, list(pnorm((u[t+1] - xi$mean)/xi$sd, lower.tail = FALSE)))
    g[[t+1]][is.nan(g[[t+1]])] = 1*(xi$mean>=u[t+1])

    # Select and evaluate N new points
    tau = (1 - g[[t+1]])*(g[[t+1]]>=0.5) + g[[t+1]]*(g[[t+1]]<0.5)
    cat(" * mean(tau) =", mean(tau),'\n')
    while(mean(tau)>eta){
      if(critere=='SUR'){
        capture.output(sur <- KrigInv::EGI(T = tail(u, 1),
                                  model = meta$model,
                                  method = "sur",
                                  fun = lsf_unif,
                                  kmcontrol = list(optim.method = optim.method),
                                  iter = 1,
                                  lower = rep(0, d), upper = rep(1, d)))
        learn_db = cbind(learn_db, t(sur$par))
        G = c(G, sur$value); Ncall = Ncall+1
        meta <- list(model = sur$lastmodel, fun = limitFunction(sur$lastmodel))
      }
      else{
        ind = which.max(tau)
        learn_db = cbind(learn_db, Y[,ind])
        G = c(G, lsf(Y[,ind])); Ncall = Ncall+1
        capture.output(meta <- trainModel(design=learn_db,
                                          response=G,
                                          kernel=kernel))
      }
      xi = meta$fun(Y)
      g[[t+1]] = pnorm((u[[t+1]] - xi$mean)/xi$sd, lower.tail = FALSE)
      g[[t+1]][is.nan(g[[t+1]])]= 1*(xi$mean>=u[[t=1]])
      tau = (1 - g[[t+1]])*(g[[t+1]]>=0.5) + g[[t+1]]*(g[[t+1]]<0.5)
      cat("   - mean(tau) =", mean(tau), "\n")
    }

    # Recompute threshold
    u[t+1] = getU(xi, g_prev, u_high = max(c(failure, xi$mean + 2*xi$sd)), u_low = u[t], failure = failure, p_0 = p_0)
    pi[[t+1]] = local({
      mf = meta$fun
      uu = u[t+1]
      function(x){
        xi = mf(x)
        res = pnorm((uu - xi$mean)/xi$sd, lower.tail = FALSE)
        res[is.nan(res)]= 1*(xi$mean>=uu)
        return(res)
      }
    })
    g[[t+1]] = pi[[t+1]](Y)

    # Regenerate MC pop
    ## calculate weights w_i \propto pi_t(\xi > u_t)/pi_{t-1}(\xi > u_{t-1})
    w = g[[t+1]]/g_prev
    m = m*mean(w)
    if(is.nan(m)){
      print(xi$sd[which(is.nan(g[[t+1]]))])
    }

    print(p + geom_contour(aes(z=z), breaks = tail(u,1)) +
            geom_contour(data = data.frame(df_plot[,1:2], z=meta$fun(t(df_plot[,1:2]))$mean), aes(z=z), breaks = tail(u,1)) +
            geom_point(data = data.frame(t(Y), z = lsf(Y)), aes(colour=z)))

    ## generate a MC sample according to weights
    ind = sample(N, size = N, replace = TRUE, prob = w/sum(w))
    Y = Y[,ind]
    ## move the new MC sample with MH
    Y = K(Y, g = pi[[t+1]])
    # xi = meta$fun(pnorm(Y))
    xi = meta$fun(Y)
    g_prev = pi[[t+1]](Y)
    t = t+1
  }
  return(list(m=m, Ncall = Ncall, meta = meta))
}
