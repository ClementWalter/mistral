# To estimate E[\alpha^2]
alpha2 <- function(PPP, model, xi, q){
  X <- PPP$final_X
  U <- PPP$final_U
  X_2 <- X
  X_2[,] <- rnorm(X)
  xi_X1 <- xi(X)
  xi_X2 <- xi(X_2)
  Kinv <- t(xi_X1$Tinv)%*%xi_X1$Tinv
  c_X1_X2 <- diag(DiceKriging::covMat1Mat2(model@covariance, X1 = t(as.matrix(X)), X2 = t(as.matrix(X_2)), nugget.flag = model@covariance@nugget.flag))
  c_xi1_xi2 <- c_X1_X2 - colSums(xi_X2$c.newdata*(Kinv%*%xi_X1$c.newdata))
  m2 <- c_xi1_xi2/xi_X2$sd^2*(U - xi_X1$mean) + xi_X2$mean
  sd2 <- xi_X2$sd^2 - c_xi1_xi2^2/xi_X1$sd^2

  # Monte Carlo estimation
  estim = list(mc = mean(pnorm(q = q, mean = m2, sd = sqrt(sd2), lower.tail = FALSE)))

  # if not precise enought, Poisson process estimation
  sigma <- 0.3
  K = function(x, sigma){
    x = as.matrix(x)
    w = matrix(rnorm(x), nrow = model@d)
    x_star = (x + sigma*w)/(sqrt(1+sigma^2))
  }
  N <- length(U)
  U_2 <- rnorm(N, mean = m2, sd = sqrt(sd2))
  M <- 0
  burnin <- 40
  while(min(U_2)<q){
    xMin <- which.min(U_2)
    Utmp <- U_2[xMin]
    cat(' * current min =', Utmp, '\n')
    cat(' * current sigma K', sigma, '\n')
    ind <- sample(c(1:N)[-xMin], 1)
    X[,xMin] <- X[,ind]
    U[xMin] <- U[ind]
    U_2[xMin] <- U_2[ind]
    acceptance <- 0
    for(i in 1:burnin){
      Xstar <- K(X[,xMin], sigma)
      xiStar <- xi(Xstar)
      Ustar <- rnorm(1, mean = xiStar$mean, sd = xiStar$sd)
      if(Ustar>q){
        X2 <- as.matrix(rnorm(model@d))
        xiX2 <- xi(X2)
        c_X1_X2 <- diag(DiceKriging::covMat1Mat2(model@covariance, X1 = t(as.matrix(Xstar)), X2 = t(X2), nugget.flag = model@covariance@nugget.flag))
        c_xi1_xi2 <- c_X1_X2 - colSums(xiX2$c.newdata*(Kinv%*%xiStar$c.newdata))
        m2 <- c_xi1_xi2/xiX2$sd^2*(Ustar - xiStar$mean) + xiX2$mean
        sd2 <- xiX2$sd^2 - c_xi1_xi2^2/xiStar$sd^2
        U2star <- rnorm(1, mean = m2, sd = sd2)
        if(U2star>Utmp){
          X[,xMin] <- Xstar
          U[xMin] <- Ustar
          U_2[xMin] <- U2star
          acceptance <- acceptance + 1
        }
      }
    }
    acceptance <- acceptance/burnin
    if(acceptance<0.2){sigma <- sigma*0.9}
    if(acceptance>0.6){sigma <- sigma*1.1}
    M <- M+1
  }
  estim$pp <- (1-1/N)^M
}
