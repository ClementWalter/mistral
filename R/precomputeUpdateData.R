#' precomputeUpdateData
#'
#' @param model a object from \code{km}
#' @param integration.points the points onto which the updated variance will be computed
#'
precomputeUpdateData <- function(model,integration.points){

  #precalculates the covariances between the integration points and the design points
  integration.points <- as.matrix(integration.points)
  colnames(integration.points) <- colnames(model@X)

  # c.olddata is the covariance between the learn.db and the points to be predicted
  c.olddata <- DiceKriging::covMat1Mat2(X1=model@X, X2=integration.points, object=model@covariance, nugget.flag=model@covariance@nugget.flag)
  # resolution en deux temps du systeme de krigeage avec n donnees
  Tinv.c.olddata <- backsolve(t(model@T), c.olddata, upper.tri=FALSE)						#
  Kinv.c.olddata <- backsolve(model@T, Tinv.c.olddata, upper.tri=TRUE)					# poids de krigeage pour integration.points sur les n db

  # model@T = upper triangular telle que t(T)%*%T = K la covariance de dim n sur les donnes
  # model@M = t(T)^{-1} %*% F
  # Kinv = K^{-1} l'inverse de la cov en dim n
  # F est la matrice des fonctions de la tendance sur les donnes, ie F = rep(1, n) en Krigeage Ordinaire
  Kinv.F <- backsolve(model@T, model@M)
  # Kinv.F = T^(-1) %*% M = T^(-1) %*% t(T)^(-1) %*% F = K^(-1) %*% F

  inv.tF.Kinv.F<- chol2inv(chol(crossprod(model@M)))
  # t(M) %*% M = t(F) %*% K^(-1) %*% F = norm.const
  # inv.tF.Kinv.F = norm.const^(-1)
  # f.integration.points <- model.matrix(model@trend.formula, data=data.frame(integration.points))
  f.integration.points <- rep(1, dim(integration.points)[1])
  first.member <- (f.integration.points-crossprod(c.olddata,Kinv.F))%*%inv.tF.Kinv.F		#
  # crossprod(c.olddata,Kinv.F) = t(c.olddata) %*% Kinv.F = t(c.olddata) %*% K^(-1) %*% F
  # first.member = (f.integration.points - t(c.olddata) %*% K^(-1) %*% F) %*% norm.const^(-1)
  # first.member = vrai F - F predit sur les nouveaux points

  #' @return A list containing the following elements
  return(list(
    #'     \item{Kinv.c.olddata}{kriging weights for the integrations.points over krig@X}
    Kinv.c.olddata=Kinv.c.olddata,
    #'     \item{Kinv.F}{The matrix product of the inverse covariance and F the matrix of
    #'     the trend functions at model@X}
    Kinv.F=Kinv.F,
    #'    \item{first.member}{}
    first.member=first.member
  ))
}
