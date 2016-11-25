#' UpdateSd
#'
#' Update kriging variance when adding new points to the DoE
updateSd <- function(X.new,
                     #' @param X.new the d x N matrix containing the points added to the model for
                     #' the update of the kriging variance.
                     integration.points.oldsd,
                     #' @param integration.points.oldsd a vector containing the standard deviation
                     #' of the points to be added to the metamodel learning database.
                     model,
                     #' @param model the current kriging model (a \code{km} object).
                     precalc.data,
                     #' @param precalc.data precomputed data from KrigInv::precomputeUpdateData.
                     integration.points
                     #' @param integration.points points where the updated sd is to be calculated.
){
  ###### Function adapted from sur_optim_parallel.R
  d <- model@d
  n <- model@n

  X.new <- t(X.new)
  krig  <- predict_nobias_km(object=model,
                             newdata=as.data.frame(X.new),
                             type="UK",
                             se.compute=TRUE,
                             cov.compute=TRUE)
  mk <- krig$mean
  sk <- krig$sd
  F.newdata <- krig$F.newdata
  c.newdata <- krig$c
  # Sigma.r matrice de covariance conditionnelle rxr des nouveaux points
  Sigma.r <- krig$cov

  # kn covariance conditionnelle de integration.points sur les n+r points ajoutes
  kn = computeQuickKrigcov(model,
                                    integration.points,
                                    X.new,
                                    precalc.data,
                                    F.newdata,
                                    c.newdata)

  ###### code chunk from predict_update_km_parallel.R
  chol.Sigma.r <- NULL
  chol.Sigma.r <- try(chol(Sigma.r),TRUE)
  if(!is.numeric(chol.Sigma.r)) return(list(error=TRUE))
  lambda_nplus.r <- kn %*% chol2inv(chol.Sigma.r)
  predict_var <- pmax(0,integration.points.oldsd^2 - rowSums(lambda_nplus.r * kn))
  predict_sd <- sqrt(predict_var)
  ######

  return(predict_sd)
  #' @return a vector containing the kriging sd at points \code{integration.points}
}
