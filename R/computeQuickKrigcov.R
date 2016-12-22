computeQuickKrigcov <- function(model,
                                integration.points,
                                X.new,
                                precalc.data,
                                F.newdata,
                                c.newdata){

  integration.points <- t(as.matrix(integration.points))
  colnames(integration.points) <- colnames(model@X)

  c.xnew.integpoints <- DiceKriging::covMat1Mat2(X1=integration.points,X2=X.new, object=model@covariance, nugget.flag=model@covariance@nugget.flag)

  second.member <- t(F.newdata - crossprod(c.newdata,precalc.data$Kinv.F))
  # second.member ~ erreur due a l'estimation de la tendance sur X.new
  # second.member = F.newdata - t(c.newdata) %*% Kinv.F = F.newdata - t(c.newdata) %*% K^(-1) %*% F

  # first.member ~ la moitiee de la variance d'erreur d'estimation de la tendance, que sur les n anciennes donnees pour integration.points
  # first.member = (f.integration.points - t(c.olddata) %*% K^(-1) %*% F) %*% norm.const^(-1)

  # cov.F est la covariance d'erreur due a l'estimation de la tendance
  cov.F <- precalc.data$first.member%*%second.member

  # c.newdata est la cov entre learn.db et les points a predire
  # crossprod(precalc.data$Kinv.c.olddata,c.newdata) = t(c.newdata) %*% K^(-1) %*% c.newdata
  # cov.std = c.xnew.integpoints - t(c.newdata) %*% K^(-1) %*% c.newdata
  kn <- c.xnew.integpoints - crossprod(precalc.data$Kinv.c.olddata,c.newdata) + cov.F
  return(kn)
}

