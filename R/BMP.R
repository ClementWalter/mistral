#' @title Bayesian Moving Particles
#'
#' @description This function runs the Bayesian Moving Particles algorithm for estimating extreme probability
#' and quantile.
#'
#' @author Clement WALTER \email{clementwalter@icloud.com}
#'
#' @aliases BMP
#'
#' @details
#' The Bayesian Moving Particles algorithm uses the point process framework for rare event to iteratively estimate
#' the conditional expectation of the (random) limit-state function, to quantify the quality of the learning
#' and to propose a new point to be added to the model with a SUR criterion.
#'
#' @note
#' Probleme should be defined in the standard space. Transformations can be made using \code{UtoX} and \code{XtoU}
#' functions.
#' 
#' @seealso
#' \code{\link{SubsetSimulation}}
#' \code{\link{MonteCarlo}}
#' \code{\link{IRW}}
#' \code{\link{MP}}
#'
#' @references
#'   \itemize{
#'    \item A. Guyader, N. Hengartner and E. Matzner-Lober:\cr
#'     \emph{Simulation and estimation of extreme quantiles and extreme
#'     probabilities}\cr
#'     Applied Mathematics \& Optimization, 64(2), 171-196.\cr
#'    \item C. Walter:\cr
#'    \emph{Moving Particles: a parallel optimal Multilevel Splitting
#'    method with application in quantiles estimation and meta-model
#'    based algorithms}\cr
#'    Structural Safety, 55, 10-25.\cr
#'    \item J. Bect, L. Li and E. Vazquez:\cr
#'    \emph{Bayesian subset simulation}\cr
#'    arXiv preprint arXiv:1601.02557
#'  }
#'
#' @examples
#' # Estimate P(g(X)<0)
#' \dontrun{p <- BMP(dimension = 2, lsf = kiureghian, q = 0, N = 100, N.iter = 30, plot = TRUE)}
#' 
#' # More extreme event
#' \dontrun{p <- BMP(dimension = 2, lsf = waarts, q = -4, N = 100, N.iter = 50, plot = TRUE)}
#' 
#' # One can also estimate probability of the form P(g(X)>q)
#' \dontrun{p <- BMP(dimension = 2, lsf = cantilever, q = 1/325, N = 100, N.iter = 30, plot = TRUE)}
#'
#' @import foreach
#' @useDynLib mistral
#' @importFrom Rcpp sourceCpp
#' @export

BMP <- function(dimension,
                #' @param dimension the dimension of the input space.
                lsf,
                #' @param lsf the function defining the RV of interest Y = lsf(X).
                q,
                #' @param q a given quantile to estimate the corresponding probability.
                N = 1000,
                #' @param N the total number of Poisson processes during the refinement step.
                N.final = N,
                #' @param N.final the total number of Poisson processes for the final alpha estimate.
                #                 N.batch = foreach::getDoParWorkers(),
                #                 #' @param N.batch the number of parallel batches for the algorithm. Each batch will then
                #                 #' have \code{N/N.batch} particles.
                N.iter = 30,
                #' @param N.iter the total number of iteration of the algorithm, ie that total number of
                #' calls to the \code{lsf} will be \code{N.DoE + N.iter*r}.
                adaptive = FALSE,
                #' @param adaptive if the algorithm should stop automatically if the stopping criterion
                #' is verified, precisely the mean probability of misclassification of the particles
                #' being over a given threshold.
                N.DoE = 5*dimension,
                #' @param N.DoE the number of points for the initial Design of Experiment
                firstDoE = "uniform",
                #' @param firstDoE default is "uniform" for a random
                #' uniform sampling over a sphere of radius \code{radius}. Also available "maximim" for a maximim LHS.
                radius = qnorm(1e-10, lower.tail = FALSE),
                #' @param radius the size of the radius of the sphere for uniform DoE or the semi length
                #' of the interval on each dimension for maximin LHS
                X,
                #' @param X (optional) a first Design of Experiemnt to be used instead of building
                #' a new DoE
                y,
                #' @param y the value of \code{lsf} on the \code{X}
                covariance = NULL,
                #' @param covariance (optional) to give a covariance kernel for the \code{km} object.
                learn_each_train = Inf,
                #' @param learn_each_train a integer: after this limit the covariance parameters are not
                #' learnt any more and model is just updated with the new datapoints.
                km.param = list(nugget.estim = TRUE, multistart = 1, optim.method="BFGS", coef.trend=q),
                #' @param km.param (optional) list of parameters to be passed to \code{DiceKriging::km}.
                burnin = 20,
                #' @param burnin a burnin parameter for Markov Chain drawing of the metamodel based
                #' Poisson process (this does not change the number of calls to \code{lsf}).
                fast = TRUE,
                #' @param fast in current implementation it appears that the call to the metamodel is
                #' faster when doing batch computation. This parameter lets do the Markov chain the other way
                #' around: instead of first selecting a starting point and then applying \code{burnin} times the
                #' transition kernel, it creates a working population by apply the kernel to all the
                #' particles and then makes some moves with the generated discretised distribution.
                sur = list(integrated=TRUE, r=1, approx.pnorm=FALSE),
                #' @param sur a list containing any parameters to be passed to \code{estimateSUR}. Default is
                #' \code{sur$integrated=TRUE} and \code{sur$r=1} for a one step ahead integrated SUR criterion.
                lower.tail = TRUE,
                #' @param lower.tail as for pxxxx functions, TRUE for estimating P(lsf(X) < q), FALSE
                #' for P(lsf(X) > q).
                save.dir,
                #' @param save.dir (optional) a directory to save the \code{X} and \code{y} at each iteration.
                plot = FALSE,
                #' @param plot to plot the DoE and the updated model.
                plot.lsf = TRUE,
                #' @param plot.lsf to plot the contour of the true \code{lsf}. Note that this requires its
                #' evaluation on a grid and should be used only on toy examples.
                plot.lab = c("x_1", "x_2"),
                #' @param plot.lab the labels of the axis for the plot.
                chi2 = FALSE,
                #' @param chi2 for a chi2 test on the number of events.
                verbose = 1,
                #' @param verbose controls the level of outputs of the algorithm.
                breaks) {
  #' @param breaks optional, for the final histogram if \code{chi2 == TRUE}.

  cat('========================================================================\n')
  cat('               Beginning of BMP algorithm\n')
  cat('========================================================================\n')
  # Fix NOTE issue with R CMD check
  x <- z <- ..level.. <- NULL
  
  if(lower.tail==TRUE){
    lsf_dec = lsf
    lsf = function(x) -1*lsf_dec(x)
    if(!is.null(km.param$coef.trend)) km.param$coef.trend = -1*km.param$coef.trend
    q <- -q
    if(!missing(y)) y = -1*y
  }

  if(plot==TRUE & dimension>2){
    message("Cannot plot in dimension > 2")
    plot <- FALSE
  }
  
  # set default values for SUR if missing
  if(is.null(sur$approx)) sur$approx = FALSE
  if(is.null(sur$approx.pnorm)) sur$approx.pnorm = FALSE

  cat(" * parallel backend registered:", foreach::getDoParRegistered(),"\n")
  cat(" * parallel backend version:", foreach::getDoParName(), "\n")
  cat("   - number of workers:", foreach::getDoParWorkers(), "\n")
  if(!missing(X)) {
    cat(" ============================================================== \n")
    cat(" BEGINNING : FIRST DoE given in inputs:", dim(X)[2],"samples \n")
    cat(" ============================================================== \n\n")
  } else {
    cat(" ============================================================== \n")
    cat(" BEGINNING : FIRST DoE with",firstDoE,"design:",N.DoE,"samples \n")
    cat(" ============================================================== \n\n")
  }
  # First DoE. Change that afterward. For now random uniform sampling
  if(missing(X)){
    X <- switch(firstDoE,
                       "uniform" = t(runifSphere(dimension, N.DoE , radius)),
                       getLHS(n = N.DoE, dimension, radius = radius))
  }
  if(missing(y)){
    y <- lsf(X)
  }
  dimnames(X) <- list(rep(c('x', 'y'), length.out = dimension))
  if(!missing(save.dir)) save(list = c('X', 'y'), file = save.dir)

  # First model
  coef.trend = km.param$coef.trend
  km.param$coef.trend <- NULL
  arg.km = c(list(design = t(X),
                  response = as.matrix(y)),
             km.param)
  if(!is.null(covariance)){
    arg.km = c(arg.km, list(
      covtype = covariance@name,
      coef.cov = covariance@range.val,
      coef.var = covariance@sd2
    ))
  }
  time.krig <- system.time({
    capture.output({model <- do.call(DiceKriging::km, arg.km)})
  })
  cat(" *", length(y), "points added to the model in",time.krig[3],"sec \n")
  cat("   - covtype =", model@covariance@name,"
   - coef.cov =",model@covariance@range.val,"
   - coef.var =",model@covariance@sd2,"
   - coef.trend =",ifelse(is.null(coef.trend), model@trend.coef, coef.trend), "\n")

  # get Ordinary Kriging optimisation
  krig <- ok(model = model, beta = coef.trend)
  xi <- krig$xi
  model.first <- model

  # Plot initialisation
  if(plot==TRUE){
    xplot <- yplot <- c(-80:80)/10
    df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = NA)
    if(plot.lsf==TRUE){df_plot$z = lsf(t(df_plot[,1:2]))}
    p <- ggplot2::ggplot(data = df_plot, aes(x,y)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlim(-8, 8) + ggplot2::ylim(-8, 8) + xlab(plot.lab[1]) + ylab(plot.lab[2]) +
      ggplot2::geom_point(data = data.frame(t(X), z = y), aes(color=z))
    if(plot.lsf==TRUE){p <- p + ggplot2::geom_contour(aes(z=z, color=..level..), breaks = q)}
    print(p)
    z_meta = xi(t(df_plot[,1:2]))
    df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(q - z_meta$mean)/z_meta$sd)
    print(p_meta <- p + ggplot2::geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = q) +
            ggplot2::geom_contour(data = df_plot_meta, aes(z=crit), breaks = 2, linetype = 'dotted'))
  }

  cat(" ================================================== \n")
  cat(" ENRICHMENT STEP:",N.iter*sur$r,"samples to be added to the DoE \n")
  cat(" ================================================== \n\n")

  alpha <- NULL
  cv <- NULL
  J <- NULL
  J_int <- NULL
  crit <- FALSE
  sur_min <- list(u=NULL, t=NULL)
  sur_stat <- data.frame(Nsur=NULL, time=NULL)

  while((N.iter > 0) && crit==FALSE){
    # Simulate N Poisson Processes
    cat(" Remaining iter. :", N.iter, "\n")
    time.PPP <- system.time({PPP = getPPP(dimension = dimension,
                                          N = N,
                                          burnin = burnin,
                                          q = q,
                                          xi = xi,
                                          fast = fast,
                                          keep_process = TRUE,
                                          verbose = verbose)})
    cat(" * Poisson process N =", N, "generated in", time.PPP[3], "sec with",foreach::getDoParWorkers(),"workers \n")

    # Get alpha the conditional expectation of p
    alpha = c(alpha, (1-1/N)^PPP$M)
    cat(" * Current alpha estimate =", tail(alpha, 1), "\n")
    cv <- c(cv, sqrt(-log(tail(alpha, 1))/N))
    cat(" * Current cv estimate =", tail(cv, 1), "\n")

    # Adaptive stopping criterion before enrichment of the model
    xi_PPP_X <- xi(cbind(PPP$X, PPP$final_X))
    J <- c(J, mean(pnorm(q, mean = tail(xi_PPP_X$mean, N), sd = tail(xi_PPP_X$sd, N)))) # J ~ var/alpha
    cat(" * Current h estimate =", tail(J, 1), "\n")

    m.LPA <- getLPA(xi_PPP_X$mean, as.numeric(names(c(PPP$U, PPP$final_U))), N, length(xi_PPP_X$mean)-N+1)[,-1];
    sd.LPA <- getLPA(xi_PPP_X$sd, as.numeric(names(c(PPP$U, PPP$final_U))), N, length(xi_PPP_X$sd)-N+1)[,-1];
    J_int <- c(J_int, sum(pnorm(rep(PPP$U, each=N), mean =c(m.LPA), sd = c(sd.LPA)))/N^2)
    cat(" * Current I estimate =", tail(J_int, 1), "\n")

    if(sur$approx.pnorm==TRUE){
      sur$J = tail(J, 1)
    }

    if(adaptive==TRUE){
      crit <- tail(J, 1) < 0.1*tail(cv, 1)
      if(crit) {
        cat("   - crit:",tail(J, 1), "<",0.1*tail(cv, 1),"\n")
        cat(" Adaptive criterion verified, stop of enrichment step\n")
        break;
      } else {
        cat("   - crit:",tail(J, 1), ">",0.1*tail(cv, 1),"\n")
        cat(" Adaptive criterion not verified\n")
      }
    }

    # SUR criterion
    sur.opt <- c(list(
      # d = dimension,
      q = q,
      krig = krig,
      PPP = PPP,
      xi_PPP_X = xi_PPP_X,
      verbose = verbose), sur)
    time.SUR <- system.time({
      SUR <- do.call(estimateSUR, sur.opt)
    })
    cat(" * SUR criterion:", length(xi_PPP_X$mean), "points tested in", time.SUR[3], "sec \n")
    sur_stat <- rbind(sur_stat, data.frame(Nsur=length(xi_PPP_X$mean), t=time.SUR[3]))

    # update model
    # NewX <- do.call(cbind,SUR)
    NewX <- as.matrix(SUR$x)
    dimnames(NewX) = NULL
    sur_min$u <- c(sur_min$u, SUR$u)
    sur_min$t <- c(sur_min$t, SUR$t)
    rownames(NewX) = rep(c('x', 'y'), l = dimension)
    cat(" * Call the lsf on the proposed point(s)\n")
    capture.output({time.lsf <- system.time({NewY <- lsf(NewX)})})
    cat(" * Lsf evaluated in", time.lsf[3], "sec \n")
    time.krig <- system.time({
      if(!is.null(covariance) | length(model@y)>learn_each_train) {
        capture.output({model <- DiceKriging::update(model,
                                                     newX = t(NewX),
                                                     newy = NewY,
                                                     cov.reestim = FALSE,
                                                     trend.reestim = is.null(coef.trend))})
      } else {
        arg.km$design = rbind(model@X, t(NewX))
        arg.km$response = c(model@y, NewY)
        arg.km$coef.trend = NULL
        arg.km$covtype = NULL
        arg.km$coef.cov = NULL
        arg.km$coef.var = NULL
        capture.output({model <- do.call(DiceKriging::km, arg.km)})
      }

      # get Ordinary Kriging optimisation
      krig <- ok(model = model, beta = coef.trend)
      xi <- krig$xi
    })
    cat(" *", length(NewY), "points added to the model in",time.krig[3],"sec \n")
    cat("   - covtype =", model@covariance@name,"
   - coef.cov =",model@covariance@range.val,"
   - coef.var =",model@covariance@sd2,"
   - coef.trend =",ifelse(is.null(coef.trend), model@trend.coef, coef.trend), "\n")

    if(!missing(save.dir)) {
      X <- model@X
      y <- model@y
      save(list = c('X', 'y'), file = save.dir)
    }

    # increment
    N.iter <- N.iter - 1

    # update plot
    if(plot==TRUE){
      p <- p + ggplot2::geom_point(data = data.frame(t(NewX), z = NewY), col = 'red')
      z_meta = xi(t(df_plot[,1:2]))
      df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(q - z_meta$mean)/z_meta$sd)
      print(p_meta <- p + ggplot2::geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = q) +
              ggplot2::geom_contour(data = df_plot_meta, aes(z=crit), breaks = 2, linetype = 'dotted'))
    }
  }

  # update kriging model with estimated trend
  if(!is.null(covariance)){
    capture.output({model <- DiceKriging::km(design = model@X, response = model@y,
                                             covtype = covariance@name,
                                             coef.cov = covariance@range.val,
                                             coef.var = covariance@sd2)})
  } else {
    arg.km$coef.trend = NULL
    arg.km$design = model@X
    arg.km$response = model@y
    capture.output({model <- do.call(DiceKriging::km, arg.km)})

  }

  # get Ordinary Kriging optimisation
  krig <- ok(model = model)
  xi <- krig$xi

  if(plot==TRUE){
    z_meta = xi(t(df_plot[,1:2]))
    df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean, crit = abs(q - z_meta$mean)/z_meta$sd)
    print(p_meta <- p + ggplot2::geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = q) +
            ggplot2::geom_contour(data = df_plot_meta, aes(z=crit), breaks = 2, linetype = 'dotted'))
  }

  # Simulate N Poisson Processes
  N <- N.final
  time.PPP <- system.time({PPP = getPPP(dimension = dimension,
                                        N = N,
                                        burnin = burnin,
                                        q = q,
                                        xi = xi,
                                        fast = fast,
                                        keep_process = TRUE,
                                        verbose = verbose)})

  cat(" * Final Poisson process N =",N,"generated in", time.PPP[3], "sec with",foreach::getDoParWorkers(),"workers \n")

  cat("========================================================================\n")
  cat("                     End of BMP algorithm \n")
  cat("========================================================================\n\n")

  # Get alpha the conditional expectation of p
  alpha = c(alpha, (1-1/N)^PPP$M)
  alpha.final = tail(alpha, 1)
  cat(" * Current alpha estimate =", tail(alpha, 1), "\n")
  cv <- c(cv, sqrt(-log(tail(alpha, 1))/N))
  cat(" * Current cv estimate =", tail(cv, 1), "\n")

  # Adaptive stopping criterion before enrichment of the model
  xi_PPP_X <- xi(cbind(PPP$X, PPP$final_X))
  J <- c(J, mean(pnorm(q, mean = tail(xi_PPP_X$mean, N), sd = tail(xi_PPP_X$sd, N)))) # J ~ var/alpha
  cat(" * Current h estimate =", tail(J, 1), "\n")

  m.LPA <- getLPA(xi_PPP_X$mean, as.numeric(names(c(PPP$U, PPP$final_U))), N, length(xi_PPP_X$mean)-N+1)[,-1];
  sd.LPA <- getLPA(xi_PPP_X$sd, as.numeric(names(c(PPP$U, PPP$final_U))), N, length(xi_PPP_X$sd)-N+1)[,-1];
  J_int <- c(J_int, sum(pnorm(rep(PPP$U, each=N), mean =c(m.LPA), sd = c(sd.LPA)))/N^2)
  cat(" * Current I estimate =", tail(J_int, 1), "\n")

  cv.seq <- cv
  cv = tail(cv.seq, 1)

  alpha_int = c( alpha.final*exp(-2*cv) , alpha.final*exp(2*cv) )
  if(lower.tail==TRUE){
    q <- -q
  }

  L_max <- q
  q_int <- c(q, q)
  L <- PPP$U
  L <- (-1)^lower.tail*L
  ecdf_MP <- local({
    lower.tail <- lower.tail
    L_max <- L_max
    L <- L;
    N <- N
    function(q) {
      if(lower.tail==TRUE){
        Nevent <- sapply(q, function(q) sum(L>q))
        Nevent[q<L_max] <- NA
        if(sum(q<L_max)>0) message(paste('Empirical cdf valid only for q >=', L_max, "NAs inserted"))
      }
      else{
        Nevent <- sapply(q, function(q) sum(L<q))
        Nevent[q>L_max] <- NA
        if(sum(q>L_max)>0) message(paste('Empirical cdf valid only for q <=', L_max, "NAs inserted"))
      }
      (1-1/N)^Nevent
    }
  })

  if(chi2 == TRUE) {
    moves <- rle(sort(as.numeric(names(PPP$U))))$lengths
    if(missing(breaks)){
      res = hist(moves, freq = FALSE)
    }
    else{
      res = hist(moves, breaks = breaks, freq = FALSE)
    }
    cont = res$counts; l = length(cont)
    chi2_p = 1:l;
    br = res$breaks
    br[1] = -Inf; br[length(br)] = Inf;
    for(i in 1:l){
      chi2_p[i] = ppois(br[i+1],mean(moves)) - ppois(br[i], mean(moves))
    }
    capture.output(chi2.test <- chisq.test(x = cont, p = chi2_p))
    p.val = 1 - pchisq(chi2.test$statistic, df = (l-2))
  }

  if(length(sur_stat)>0) sur_stat$N <- N

  cat("   - alpha =",alpha.final,"\n")
  cat("   - cv =", cv, "\n")
  cat("   - q =",q,"\n")
  cat("   - 95% confidence intervalle :",alpha_int[1],"< alpha <",alpha_int[2],"\n")
  cat("   - Total number of calls =",length(model@y),"\n")
  if(chi2 == TRUE) {cat("   - Chi-2 test =", chi2.test$statistic,"; p-value =", p.val,"\n")}

  #' @return An object of class \code{list} containing the outputs described below:
  res = list(alpha = alpha.final,
             #' \item{alpha}{the estimated conditional expectation of the probability.}
             alpha.seq = alpha,
             #' \item{alpha.seq}{the sequence of estimated alpha during the refinement step.}
             cv = cv,
             #' \item{cv2}{an estimate of the squarred coefficient of variation of alpha.}
             cv.seq = cv.seq,
             #' \item{cv.seq}{the sequence of the estimated coefficients of variations.}
             h = J,
             #' \item{h}{the sequence of the estimated upper bound of the conditional variance
             #' divided by estimated alpha.}
             I = J_int,
             #' \item{I}{the sequence of the estimated integrated h.}
             sur_min = sur_min,
             #' \item{sur_min}{a list containing the the sequence of corresponding thresholds and
             #' -log probability of the sample minimising the SUR criterion.}
             sur_stat = sur_stat,
             #' \item{sur_stat}{a list containing at each iterations number of points tried for the SUR
             #' criterion as well as the computational spent.}
             q = q,
             #' \item{q}{the reference quantile for the probability estimate.}
             ecdf = ecdf_MP,
             #' \item{ecdf}{the empirical cdf, i.e. the estimation of the function q -> E(alpha(q)).}
             L_max = L_max,
             #' \item{L_max}{the farthest state reached by the random process. Validity range
             #' for the \code{ecdf} is then (-Inf, L_max] or [L_max, Inf).}
             PPP = PPP,
             #' \item{PPP}{the last Poisson process generated with \code{N.final} particles.}
             meta_fun = function(x) {
               g = xi(x)
               g$mean = (-1)^(lower.tail)*g$mean
               return(g)
             },
             #' \item{meta_fun}{the metamodel approximation of the \code{lsf}. A call output is a
             #' list containing the value and the standard deviation.}
             model = model,
             #' \item{model}{the final metamodel. An S4 object from \pkg{DiceKriging}. Note
             #' that the algorithm enforces the problem to be the estimation of P[lsf(X)>q]
             #' and so using \sQuote{predict} with this object will return inverse values if
             #' \code{lower.tail==TRUE}; in this scope prefer using directly \code{meta_fun} which
             #' handles this possible issue.}
             model.first = model.first,
             #' \item{model.first}{the first metamodel with the intial DoE.}
             alpha_int = alpha_int)
  #' \item{alpha_int}{a 95\% confidence intervalle on the estimate of alpha.}
  if(chi2 == TRUE) {res = c(res, list(
    moves = moves,
    #' \item{moves}{a vector containing the number of moves for each one of the \code{N.batch} particles.}
    chi2 = chi2.test))}
  #' \item{chi2}{the output of the chisq.test function.}

  return(res)
}
