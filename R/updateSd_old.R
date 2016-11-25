#' UpdateSd.old
updateSd.old <- function(X.new,
                     #' @param X.new the d x N matrix containing the points added to the model for
                     #' the update of the kriging variance.
                     newdata.oldsd,
                     #' @param newdata.oldsd a vector containing the standard deviation
                     #' of the points to be added to the metamodel learning database.
                     model,
                     #' @param model the current kriging model (a \code{km} object).
                     precalc.data,
                     #' @param precalc.data precomputed data from KrigInv::precomputeUpdateData.
                     integration.points
                     #' @param integration.points points where the updated sd is to be calculated.
){
  ###### Code chunk from predict_nobias_km.R
  #newdata : n x d
  object = model
  X <- object@X
  newdata = t(X.new)

  newdata <- as.matrix(newdata)
  dim.newdata <- dim(newdata)
  m <- dim.newdata[1]
  d.newdata <- dim.newdata[2]

  if (!identical(d.newdata, object@d)) stop(paste("newdata must have the same numbers of columns (",d.newdata,") than the experimental design (",object@d,")"))
  if (!identical(colnames(newdata), colnames(X))) colnames(newdata) <- colnames(X)

  T <- object@T
  M <- object@M
  F.newdata <- stats::model.matrix(object@trend.formula, data=data.frame(newdata))
  c.newdata <- covMat1Mat2(X1=X, X2=newdata, object=object@covariance, nugget.flag=object@covariance@nugget.flag)
  Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri=FALSE)

  if (object@covariance@nugget.flag) {
    total.sd2 <- object@covariance@sd2 + object@covariance@nugget
  } else total.sd2 <- object@covariance@sd2

  C.newdata <- covMatrix(X=newdata, object=object@covariance)[[1]]
  cond.cov <- C.newdata - crossprod(Tinv.c.newdata)

  # if (type=="UK") {
  T.M <- chol(t(M)%*%M)
  s2.predict.mat <- backsolve(t(T.M), t(F.newdata - t(Tinv.c.newdata)%*%M), upper.tri=FALSE)
  cond.cov <- cond.cov + crossprod(s2.predict.mat)
  # }

  Sigma.r = cond.cov
  ######

  kn = KrigInv::computeQuickKrigcov(model,
                                    integration.points,
                                    newdata,
                                    precalc.data,
                                    F.newdata,
                                    c.newdata)

  ###### code chunk from predict_update_km_parallel.R
  chol.Sigma.r <- NULL
  chol.Sigma.r <- try(chol(Sigma.r),TRUE)
  if(!is.numeric(chol.Sigma.r)) return(list(error=TRUE))

  lambda_nplus.r <- kn %*% chol2inv(chol.Sigma.r)

  predict_var <- pmax(0,newdata.oldsd*newdata.oldsd - rowSums(lambda_nplus.r * kn))
  predict_sd <- sqrt(predict_var)
  ######

  return(predict_sd)
}
