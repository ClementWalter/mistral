
predict_nobias_km <- function(object, newdata, type="UK", se.compute=TRUE, cov.compute=FALSE,low.memory=FALSE,...) {
  #newdata : n x d

  X <- object@X
  y <- object@y

  newdata <- as.matrix(newdata)
  dim.newdata <- dim(newdata)
  m <- dim.newdata[1]
  d.newdata <- dim.newdata[2]

  T <- object@T # K = t(T) %*% T
  z <- object@z # z = t(T)^(-1) %*% (y - F %*% beta)
  M <- object@M # M = t(T)^(-1) %*% F
  beta <- object@trend.coef

  # F.newdata <- model.matrix(object@trend.formula, data=data.frame(newdata))
  F.newdata <- as.matrix(rep(1,m)) # in OK, trend ~ 1
  y.predict.trend <- F.newdata%*%beta


  # compute c(x) for x = newdata ; remark that for prediction (or filtering), cov(Yi, Yj)=0  even if Yi and Yj are the outputs related to the equal points xi and xj.
  c.newdata <- DiceKriging::covMat1Mat2(X1=X, X2=newdata, object=object@covariance, nugget.flag=object@covariance@nugget.flag)

  # t(T) %*% T %*% lambda = K %*% lambda = c.newdata
  Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri=FALSE) # Tinv.c.newdata = t(T)^(-1) %*% c.newdata

  # y.predict.complement est la prediction processus stationnaire
  y.predict.complement <- crossprod(Tinv.c.newdata,z) # t(Tinv.c.newdata) %*% z = t(c.newdata) %*% K %*% (y - F %*% beta)

  # y.predict = trend + prediction stationnaire
  y.predict <- y.predict.trend + y.predict.complement
  y.predict <- as.numeric(y.predict)


  output.list <- list(mean = y.predict,
                      c = c.newdata,
                      Tinv.c = Tinv.c.newdata,
                      F.newdata = F.newdata)


  # s2.predict.1 ~ erreur de krigeage simple
  # s2.predict.1 = t(Tinv.c.newdata) %*% Tinv.c.newdata = t(c.newdata) %*% T^(-1) %*% t(T)^(-1) %*% c.newdata = t(c.newdata) %*% K^(-1) %*% c.newdata

  # ecart type du model sigma_0
  if (object@covariance@nugget.flag) {
    total.sd2 <- object@covariance@sd2 + object@covariance@nugget
  } else total.sd2 <- object@covariance@sd2

  # M est la "racine" de la "constante" de normalisation
  # t(M) %*% M = t(F) %*% T^(-1) %*% t(T)^(-1) %*% F = t(F) %*% K^(-1) %*% F = norm.const

  T.M <- chol(crossprod(M))
  # crossprod(Tinv.c.newdata,M) = t(Tinv.c.newdata) %*% M = t(c.newdata) %*% T^(-1) %*% t(T)^(-1) %*% F = t(c.newdata) %*% K^(-1) %*% F
  # F.newdata - crossprod(Tinv.c.newdata,M) = F.newdata - t(c.newdata) %*% K^(-1) %*% F
  # s2.predict.mat = t(T.M)^(-1) %*% ( F.newdata - t(F) %*% K^(-1) %*% c.newdata )
  s2.predict.mat <- backsolve(t(T.M), t(F.newdata - crossprod(Tinv.c.newdata,M)) , upper.tri=FALSE)

  # s2.predict.2 ~ variance dues a l'estimation de la moyenne
  # s2.predict.2 = (t(F.newdata) - t(c.newdata) %*% K^(-1) %*% F ) %*% norm.const^(-1) %*% (F.newdata - t(F) %*% K^(-1) c.newdata)

  C.newdata <- covMatrix(X=newdata, object=object@covariance)[[1]]

  # crossprod(s2.predict.mat) = (t(c.newdata) %*% K^(-1) %*% F - t(F.newdata)) %*% norm.const^(-1) %*% (F.newdata - t(F) %*% K^(-1) c.newdata)

  cond.cov <- C.newdata - crossprod(Tinv.c.newdata) + crossprod(s2.predict.mat)
# cond.cov = cov - variance KS + variance estim tendance
  output.list$cov <- cond.cov

  s2.predict <- diag(cond.cov) + total.sd2

  output.list$sd <- sqrt(s2.predict)

  return(output.list)
}
