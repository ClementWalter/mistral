## -----------------------------------------------------------------------------
## Fonction UtoX
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : G. DEFAUX
##    CEA
## -----------------------------------------------------------------------------

##
## FONCTION POUR MODIFIER LA MATRICE DE CORRELATION
##   * ENTREE : MATRICE DE CORRELATION DE SPEARMAN
##   * SORTIE : MATRICE MODIFIEE POUR DISTRIBUTION JOINTE
#' @export
ModifCorrMatrix = function(Rs){
  Rho <- Rs
  R0 <- Rho
  d <- dim(Rho)[1]
  for(i in 1:d){
    for(j in 1:i-1){
      R0[i,j] <- 2.0*sin(Rho[i,j]*pi/6.0)
      R0[j,i] <- R0[i,j]
    }
  }
  return(R0)
}

##
## CALCUL DES PARAMETRES INTERNES A PARTIR DES MOMENTS
## 
#' @export
ComputeDistributionParameter = function(margin){
  ##
  ## LOGNORMAL DISTRIBUTION
  ##
  if(margin$type == "Lnorm"){
    if( is.null(margin$MEAN) ){
      margin$MEAN = exp( margin$P1 + 0.5*(margin$P2^2) )
      margin$STD  = margin$MEAN * sqrt( exp(margin$P2^2)-1.0 )
    }
    if( is.null(margin$P1) ){
      margin$P1 = log(   margin$MEAN^2 
                         / sqrt(margin$MEAN^2 + margin$STD^2) )
      margin$P2 = sqrt( log( 1.0 + (margin$STD/margin$MEAN)^2 ) ) 
    }
  }
  ##
  ## UNIFORM DISTRIBUTION
  ##
  if(margin$type == "Unif"){
    if( is.null(margin$MEAN) ){
      margin$MEAN = 0.5*( margin$P1 + margin$P2 )
      margin$STD  = ( margin$P2 - margin$P1 )/ sqrt(12.0)
    }
  }
  ##
  ## GUMBEL DISTRIBUTION
  ##
  if(margin$type == "Gumbel"){
    if( is.null(margin$P1) ){
      margin$P2 = pi / (sqrt(6.0)*margin$STD)
      margin$P1 = margin$MEAN - 0.5772156649/margin$P2
    }
  }
  return(margin)
}


##
## TRANSFORMATION ISO-PROBABILISTE
##
#' @export
UtoX = function(U, input.margin, L0){
  d <- length(input.margin)
  Z <- as.matrix(U)
  if(identical(L0,diag(d))==FALSE){
    Z <- L0%*%as.matrix(U)
  }
  X <- ZtoX(Z, input.margin)
  return(X)
}

#' @title From X to standard space
#' 
#' @description XtoU lets transform datapoint in the original space X to the standard Gaussian
#' space U with isoprobalisitc transformation.
#' 
#' @author Clement WALTER \email{clement.walter@cea.fr}
#' 
#' @seealso \code{UtoX}
#' 
#' @export
XtoU = function(X, 
                #' @param X the matrix d x n of the input points
                input.margin,
                #' @param input.margin A list containing one or more list defining the marginal distribution
                #' functions of all random variables to be used
                L0
                #' @param L0 the lower matrix of the Cholesky decomposition of correlation matrix R0
                #' (result of \code{\link{ModifCorrMatrix}})
){
  d <- length(input.margin)
  
  U <- XtoZ(X, input.margin)
  if(!identical(L0, diag(1, d))) {
    U = chol2inv(chol(L0))%*%U
  }
  return(U)
}
##
## TRANSFORMATION DES MARGINALES
##

ZtoX = function(Z, input.margin){
  d <- length(input.margin)
  # X <- matrix(NA,d,1)
  X <- matrix(NA, d, dim(Z)[2])
  for (i in 1:d){
    if(input.margin[[i]]$type == "Norm"){
      X[i,] <- input.margin[[i]]$MEAN + input.margin[[i]]$STD*Z[i,] 
    }
    if(input.margin[[i]]$type == "Lnorm"){
      X[i,] <- exp( input.margin[[i]]$P1 + input.margin[[i]]$P2*Z[i,] ) 
    }
    if(input.margin[[i]]$type == "Unif"){
      X[i,] <- qunif( p=pnorm(Z[i,]), min=input.margin[[i]]$P1, max=input.margin[[i]]$P2 )
    }
    if(input.margin[[i]]$type == "Gumbel"){
      X[i,] <- input.margin[[i]]$P1 - log(-log(pnorm(Z[i,])))/input.margin[[i]]$P2
    }
    if(input.margin[[i]]$type == "Weibull"){
      X[i,] <- qweibull( p=pnorm(Z[i,]), shape=input.margin[[i]]$P1, scale=input.margin[[i]]$P2 )
    }
    if(input.margin[[i]]$type == "Gamma"){
      X[i,] <- qgamma( p=pnorm(Z[i,]), shape=input.margin[[i]]$P1, scale=input.margin[[i]]$P2 )
    }
    if(input.margin[[i]]$type == "Beta"){
      X[i,] <- qbeta( p=pnorm(Z[i,]), shape1=input.margin[[i]]$P1, shape2=input.margin[[i]]$P2)
    }
    # else {stop(cat("ERROR in ZtoX : type ",input.margin[[i]]$type," is NOT DEFINED",sep=""))}
  }
  return(X)
}

XtoZ = function(X, input.margin){
  d <- length(input.margin)
  # X <- matrix(NA,d,1)
  Z <- matrix(NA, d, dim(X)[2])
  for (i in 1:d){
    if(input.margin[[i]]$type == "Norm"){
      (X[i,] - input.margin[[i]]$MEAN)/input.margin[[i]]$STD -> Z[i,] 
    }
    if(input.margin[[i]]$type == "Lnorm"){
      (log(X[i,]) - input.margin[[i]]$P1)/input.margin[[i]]$P2 -> Z[i,]
    }
    if(input.margin[[i]]$type == "Unif"){
      Z[i,] <- qnorm((X[i,] - input.margin[[i]]$P1)/(input.margin[[i]]$P2 - input.margin[[i]]$P1))
    }
    if(input.margin[[i]]$type == "Gumbel"){
      qnorm(exp(-exp(-(X[i,] - input.margin[[i]]$P1)*input.margin[[i]]$P2))) -> Z[i,]
    }
    if(input.margin[[i]]$type == "Weibull"){
      qnorm(stats::pweibull(q=X[i,], shape=input.margin[[i]]$P1, scale=input.margin[[i]]$P2 )) -> Z[i,]
    }
    if(input.margin[[i]]$type == "Gamma"){
      qnorm(stats::pgamma(X[i,], shape=input.margin[[i]]$P1, scale=input.margin[[i]]$P2 )) -> Z[i,]
    }
    if(input.margin[[i]]$type == "Beta"){
      qnorm(stats::pbeta(X[i,], shape1=input.margin[[i]]$P1, shape2=input.margin[[i]]$P2)) -> Z[i,]
    }
    # else {stop(cat("ERROR in ZtoX : type ",input.margin[[i]]$type," is NOT DEFINED",sep=""))}
  }
  return(Z)
}

