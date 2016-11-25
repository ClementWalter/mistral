getXi <- function(object){
  T <- object@T
  Tinv <- backsolve(t(T), diag(dim(T)[1]), upper.tri = FALSE)
  total.sd2 <- object@covariance@sd2
  if (object@covariance@nugget.flag) {
    total.sd2 <- total.sd2 + object@covariance@nugget
  }

  xi <- local({
    covariance <- object@covariance
    nugget.flag = covariance@nugget.flag
    X <- object@X
    beta <- object@trend.coef
    z <- object@z
    M <- object@M
    Tinv <- Tinv
    T.M2 <- sum(M^2)
    total.sd2 <- total.sd2
    function(x){
      newdata <- t(as.matrix(x))
      c.newdata <- DiceKriging::covMat1Mat2(covariance, X1 = X, X2 = newdata, nugget.flag = nugget.flag)
      F.newdata <- as.matrix(rep(1, dim(newdata)[1]))

      Tinv.c.newdata <- Tinv%*%c.newdata
      y.predict <- c(beta + t(Tinv.c.newdata) %*% z)

      s2.predict.1 <- c(rep(1, dim(Tinv.c.newdata)[1])%*%Tinv.c.newdata^2)
      s2.predict.2 <- c(F.newdata - t(Tinv.c.newdata) %*% M)^2/T.M2
      s2.predict <- c(pmax(total.sd2 - s2.predict.1 + s2.predict.2, 0))
      return(list(mean=y.predict, sd = sqrt(s2.predict), lambda = Tinv.c.newdata, Tinv = Tinv, c.newdata = c.newdata))
    }
  })
  return(xi)
}
