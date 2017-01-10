#' @title Class of Ordinary Kriging
#' 
#' @description An implementation of Ordinary Kriging based upon a km-class object that should be
#' faster than usual predict method.
#' 
#' @author Clement WALTER \email{clementwalter@icloud.com}
#' 
#' @details The Ordinary Kriging is a special case of kriging where the trend is supposed to be
#' and unknown constant. Consequently some linear algebra operations can be reduced by knowning
#' that the vector of parameter \code{beta} is indeed a real.
#' 
#' The ok class defines three functions: \code{xi} the kriging predictor, \code{updateSd} and
#' \code{updateSdfast} two methods for updating the kriging variance when some poitns are
#' virtually added to the model. These two last functions differ in their implementation: the
#' first one allows for the user to specify which are the predicted points and which are the
#' added points. The second one outputs a matrix where the kriging variances of all the points
#' is updated when each one is iteratively added the the Design of Experiments.
#' 
#' The faster between looping \code{updateSd} and using \code{updateSdfast} is indeed problem
#' dependent (depending on parallel computer, size of the data, etc.) and should be
#' benchmark by the user.
#' 
#' @examples
#' # Generate a dataset
#' X <- data.frame(x1 = rnorm(10), x2 = rnorm(10))
#' y <- cos(sqrt(rowSums(X^2)))
#' 
#' # Learn a model
#' krig <- DiceKriging::km(design=X, response=y)
#' 
#' # Create Ordinary Kriging object
#' OK <- ok(krig)
#' 
#' # Microbenchmark
#' # create a dataset
#' X = data.frame(x1 = rnorm(100), x2 = rnorm(100))
#' microbenchmark::microbenchmark(OK$xi(t(X)), predict(krig, X, type="UK"))
#' 
#' # Check identical results
#' X <- rnorm(2)
#' OK$xi(X)[c('mean', 'sd')]
#' predict(krig, data.frame(x1=X[1], x2=X[2]), type="UK")[c('mean', 'sd')]
#' 
#' @export
ok <- function(model,
               #' @param model a kriging model object from \code{DiceKriging::km-class}
               beta = NULL
               #' @param beta the trend of the model
){

  # get data from km model
  d <- model@d # dimension of the problem
  n <- model@n # number of data
  X <- model@X # matrix of the data: size n x d
  y <- model@y # real valued response: length n
  if(is.null(beta)){
    beta <- model@trend.coef # estimation of the trend beta = t(rep(1, n)) %*% Kinv %*% y * sigma_beta^2
  }
  y_centred <- y - beta
  T <- model@T # choleski decompostion of K: t(T) %*% T = K the covariance of the n data
  cov_function <- model@covariance
  sigma_0_2 <- model@covariance@sd2 + ifelse(model@covariance@nugget.flag, model@covariance@nugget, 0)
  sigma_0 <- sqrt(sigma_0_2)

  # get quantities used for OK
  Kinv <- chol2inv(T)
  sigma_beta = sqrt(1/sum(Kinv))

  # define kriging predictor function
  xi <- local({
    X <- X
    Kinv <- Kinv
    y_centred <- y_centred
    beta <- beta
    sigma_0 <- sigma_0
    sigma_beta <- sigma_beta
    cov_function <- cov_function

    function(x){
      # calculate covariance between x and data
      kn <- DiceKriging::covMat1Mat2(object = cov_function, X1 = X, X2 = t(x), nugget.flag = TRUE) # kn is a matrix: n x dim(x)[2]

      # calculate most computationnally intensive product
      Kinv_kn <- Kinv %*% kn # matrix n x dim(x)[2]

      # calculate kriging mean
      mean <- c(t(Kinv_kn) %*% y_centred + beta)

      # calculate kriging sd
      sd_sk_2 <- sigma_0^2 - colSums(kn*Kinv_kn)
      if(sum(sd_sk_2<0)>0){
        # message("Found negative squarred standard deviation(s):", sd_sk_2[sd_sk_2<0],
        #         "rounded to 0 (corresponding squarred coefficient of variation(s) =", sd_sk_2[sd_sk_2<0]/mean[sd_sk_2<0]^2)
        sd_sk_2[sd_sk_2<0] <- 0
      }
      sd_beta_2 <- sigma_beta^2*(colSums(Kinv_kn) - 1)^2
      sd <- sqrt(sd_sk_2 + sd_beta_2)
      return(list(mean=mean,sd=sd, kn=kn, Kinv_kn=Kinv_kn))
    }
  })

  # define update sd when new points are added for given x
  updateSd <- local({
    Kinv <- Kinv
    y_centred <- y_centred
    beta <- beta
    sigma_0 <- sigma_0
    sigma_beta <- sigma_beta
    cov_function <- cov_function

    function(x, Xnew, xi_x, xi_Xnew){
      # Xnew: new data added to the DoE: matrix d x r

      r <- dim(as.matrix(Xnew))[2]

      # get covaraince matrix of the r new data
      if(r==1) {
        K_rr <- sigma_0^2
      } else {
        K_rr <- DiceKriging::covMatrix(object = cov_function, X = t(Xnew))$C # get a r x r matrix
      }

      # get covariance between new data and old data
      K_nr <- as.matrix(xi_Xnew$kn)

      # calculate most computationnally intensive product
      Kinv_Knr <- as.matrix(xi_Xnew$Kinv_kn)

      # get vector for sd_beta_2
      sd_beta <- (colSums(Kinv_Knr) - 1)
      if(r==1){
        sd_beta_2 <- sigma_beta^2*sd_beta^2
      } else {
        sd_beta_2 <- sigma_beta^2*as.matrix(sd_beta)%*%sd_beta
      }

      # get conditional covariance matrix r x r
      Sigma_rr <- K_rr - t(K_nr) %*% Kinv_Knr + sd_beta_2

      # get covariance between x and the r data points
      k_r <- DiceKriging::covMat1Mat2(object = cov_function, X1 = t(Xnew), X2 = t(x), nugget.flag = TRUE) # chi_r is matrix a r x dim(x)[2]

      # get covariance between x and the n old data
      k_n <- as.matrix(xi_x$kn)

      # get conditional covariance
      chi_r = k_r - t(Kinv_Knr) %*% k_n + sigma_beta^2*as.matrix(sd_beta)%*%(colSums(xi_x$Kinv_kn) - 1)

      # inverse covariance matrix
      if( r==1) {
        Sigma_rr_inv = 1/Sigma_rr
      } else {
        Sigma_rr_inv <- chol2inv(chol(Sigma_rr))
      }

      # get SK variance on the r new data
      if(r==1){
        sd_sk_2 <- colSums(chi_r^2)*Sigma_rr_inv
      } else {
        sd_sk_2 <- colSums(chi_r*(Sigma_rr_inv%*%chi_r))
      }
      # get updated sd
      sd <- sqrt(pmax(xi_x$sd^2 - sd_sk_2, 0))

      return(sd)
    }
  })

  updateSdfast <- local({
    sigma_beta <- sigma_beta
    cov_function <- cov_function

    function(x, xi_x){

      Nsur <- dim(x)[2]
      Sigma_r_inv <- matrix(1, nrow = Nsur, ncol = Nsur)/xi_x$sd^2

      K_xx <- DiceKriging::covMatrix(object = cov_function, X = t(x))$C

      chi_r <- K_xx - t(xi_x$kn) %*% xi_x$Kinv_kn + sigma_beta^2 * (colSums(xi_x$Kinv_kn) - 1) %*% t((colSums(xi_x$Kinv_kn) - 1))

      sd_old_2 <- t(matrix(1, nrow = Nsur, ncol = Nsur)*xi_x$sd^2)
      sd_nr_2 <- pmax(sd_old_2 - chi_r^2*Sigma_r_inv, 0)

      return(sd_nr_2)
    }
  })

  #' @return An object of S3 class 'ok' containing
  res <- list(Kinv = Kinv,
              #'               \item{Kinv}{the inverse of the covariance matrix of the data}
              beta = beta,
              #'              \item{beta}{the estimated coefficient of the trend}
              y_centred = y_centred,
              #'              \item{y_centred}{the data centred according to the estimated trend}
              sigma_beta = sigma_beta,
              #'              \item{sigma_beta}{the standard deviation of the estimation of beta}
              xi = xi,
              #'              \item{xi}{the kriging predictor}
              updateSd = updateSd,
              #'              \item{updateSd}{a function to calculate the updated kriging variance when
              #'              \code{Xnew} points are added to the Design of Experiments}
              updateSdfast = updateSdfast
              #'              \item{updateSdfast}{a function to calculate the update kriging variance
              #'              when the SUR criterion is minimised over a population which is also the one
              #'              used to estimate it.}
  )
  class(res) <- 'ok'
  return(res)
}
