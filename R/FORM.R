#' @title First-order reliability method
#' 
#' @description The First-Order Reliability Method computes an estimation of the failure
#' probability by approximating the limit-state function at the Most Probable Failure Point
#' with a hyperplane.
#' 
#' @details The FORM method has to be used in the standard Gaussian input space. It is designed
#' to estimate probability of the form \eqn{P[g(\mathbf{X}) < 0]}{P[g(X) < 0]} with g the limit-state function.
#' This function has to be modified accordingly to fit into this framework
#' 
#' Furthermore, it should be able to handle matrix input of column vectors. See the mistral vignette
#' for more info about \code{lsf} definition
#' 
#' @author Vincent MOUTOUSSAMY and Clement WALTER \email{clementwalter@icloud.com}
#' 
#' @references
#' \itemize{
#' \item
#' O. Ditlevsen and H.O. Madsen. Structural reliability methods, Wiley, 1996\cr
#' \item
#' M. Lemaire, A. Chateauneuf and J. Mitteau. Structural reliability, Wiley Online Library, 2009.\cr
#' }
#' 
#' @examples 
#' \dontrun{
#' # u.dep is a starting point for the research of the Most Probable Failing Point
#' # N.calls is a total number of calls
#' form <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000,
#'                      u.dep = c(0,0))
#' form$p
#' 
#' # use IS=TRUE to use an Importance Sampling scheme with a Gaussian standard
#' # proposal distribution centred at the MPFP
#' form.IS <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000,
#'                         u.dep = c(0,0),
#'                         IS = TRUE)
#' form.IS$p
#' }
#' 
#' @return A list containing the following objects
#' \item{p}{Failure probability}
#' \item{indice.reliab}{Reliability index}
#' \item{Ncall}{Number of calls to f}
#' \item{Design.Point}{Coordinates of the design point}
#' \item{fact.imp}{Importance factors}
#' \item{variance}{Standard error of the probability estimator (if IS = TRUE)}
#' \item{Interval.conf}{Confidence interval of the estimator at 0.95 (if IS = TRUE)}
#' \item{DOE}{List which contains the design of experiments}
#'
#' @import ggplot2 
#'
#' @export
FORM <- function(dimension,
                 #' @param dimension the dimension of the input space.
                 lsf,
                 #' @param lsf the limit-state function.
                 u.dep = rep(0, dimension),
                 #' @param u.dep the starting point for the MPFP search.
                 N.calls = 1e2,
                 #' @param N.calls the total number of calls for the whole algorithm.
                 eps = 1e-7,
                 #' @param eps stopping criterion: distance of two points between two iterations.
                 Method = "HLRF",
                 #' @param Method choice of the method to search the design point: "AR" for Abdo-Rackwitz
                 #' and "HLRF" for Hasofer-Lindt-Rackwitz-Fiessler.
                 IS = FALSE,
                 #' @param IS "TRUE" for using importance Sampling method with an standard Gaussian importance
                 #' density centred at the MPFP.
                 IS.ratio = 0.5,
                 #' @param IS.ratio ratio of N.calls for the search of the design point by FORM. Default = 0.5.
                 #' 1-\code{IS.ratio} = the remaining ratio to be used for importance sampling.
                 plot = FALSE,
                 #' @param plot to plot the generated samples.
                 plot.lsf = FALSE,
                 #' @param plot.lsf a boolean indicating if the \code{lsf} should be added to the
                 #' plot. This requires the evaluation of the \code{lsf} over a grid and
                 #' consequently should be used only for illustation purposes.
                 plot.lab = c('x_1', 'x_2')
                 #' @param plot.lab the x and y labels for the plot.
){
  
  # Few error messages : 
  if (mode(lsf) != "function") stop("lsf is not a function")
  if(plot==TRUE & dimension>2){
    message("Cannot plot dimension>2 => plot <- FALSE")
    plot <- FALSE
  }
  
  # Fix check NOTE
  x <- y <- z <- NULL
  
  #------ Gradient------#
  GRAD <- function(func, X, f.X){
    
    X <- as.matrix(X)
    dX <- X[,rep(1, dimension)] + 1e-6*diag(dimension)
    f.dX <- lsf(dX)
    dfunc <- (f.dX - f.X)/1e-6
    
    return(list(g = dfunc,
                dX = dX,
                f.dX = f.dX))
  }
  
  #----------------------------------------------------------------------#
  
  List.res <- list()				 # List of the output
  u.dep <- as.matrix(u.dep)
  dimnames(u.dep) = list(rep(c('x', 'y'), length.out = dimension))
  u.new    <- 0
  cp       <- 0   			       #Numbers of calls to f
  fact.imp <- 0					 #importance factor
  N.FORM   <- ifelse(IS, floor(N.calls*IS.ratio), N.calls) 		 #Number of calls reserved to FORM
  
  u.ans <- lsf(u.dep); cp <- cp + 1
  DOE <- list(x = u.dep, y = u.ans)                       # Design of experiments
  
  if(plot==TRUE){
    xplot <- yplot <- c(-80:80)/10
    df_plot = data.frame(expand.grid(x=xplot, y=yplot))
    df.points <- data.frame(t(DOE$x), z = DOE$y)
    p <- ggplot(data = df.points, aes(x,y)) +
      theme(legend.position = "none") +
      xlim(-8, 8) + ylim(-8, 8) + xlab(plot.lab[1]) + ylab(plot.lab[2]) +
      geom_point(aes(color=z))
    if(plot.lsf==TRUE){
      zplot = lsf(t(df_plot))
      p <- p + geom_contour(data = data.frame(df_plot, z = zplot), aes(z=z, color=z), breaks = 0)
    }
    print(p)
  }
  
  res     <- GRAD(lsf, u.dep, u.ans); cp <- cp + dimension
  
  g <- res[[1]]
  
  DOE$x <- cbind(DOE$x, res$dX)
  DOE$y <- c(DOE$y, res$f.dX)
  
  err   <- 1					       #distance of two tested points
  
  
  
  #----------------------------------Algorithme HLRF :----------------- #
  if(Method == "HLRF"){
    
    a <- g/sqrt(sum(g^2))
    
    while( (err > eps) & (cp <= (N.FORM - (1 + 2*dimension)) ) ){ # modif BIS
      
      u.new <- (sum(u.dep*a)*a - u.ans*a/sqrt(sum(g^2)))
      err   <- sqrt(sum((u.dep-u.new)^2))
      u.dep <- u.new
      u.ans <- lsf(u.dep); cp <- cp + 1
      DOE$x <- cbind(DOE$ x, u.dep)
      DOE$y <- c(DOE$y, u.ans)
      
      res   <- GRAD(lsf, u.dep, u.ans); cp <- cp + dimension
      
      g <- res[[1]]
      
      DOE$x <- cbind(DOE$x, res$dX)
      DOE$y <- c(DOE$y, res$f.dX)
      
      a <- g*ifelse( sum(g^2) == 0, 0, 1/sqrt(sum(g^2)))
    }
    
  }
  #-----------------------------End of the algorithm HLRF :--------------------------#
  
  #-----------------------------Algorithm Abdo-Rackwitz :--------------------------#
  if(Method == "AR"){
    
    lambda <- 2*(u.ans - (t(g)%*%u.dep))/(g%*%g)
    
    while( (err > eps) & (cp <= (N.FORM - (1 + 2*dimension)) ) ){ # modif BIS
      
      u.new <- -lambda*g/2
      err   <- sqrt(t(u.dep-u.new)%*%(u.dep-u.new))
      u.dep <- u.new
      
      u.ans <- lsf(u.dep); cp <- cp + 1      
      DOE$x <- cbind(DOE$ x, u.dep)
      DOE$y <- c(DOE$y, u.ans)   
      
      res   <- GRAD(lsf, u.dep, u.ans); cp <- cp + dimension
      
      g <- res[[1]]
      
      DOE$x <- cbind(DOE$x, res$dX)
      DOE$y <- c(DOE$y, res$f.dX)
      
      lambda <- 2*(u.ans - sum(g*u.dep))/sum(g^2)
      if(cp >= N.FORM){break}
    }
    
  }
  #------------------------End of the algorithm of Abdo-Rackwitz----------------------#
  
  if(plot==TRUE){
    p$data <- data.frame(t(DOE$x), z = DOE$y, row.names = NULL, check.names = FALSE)
    print(p)
  }
  
  a <- g/sqrt(sum(g^2))
  B <- -a%*%u.dep			
  P <- c(pnorm(-B))
  
  fact.imp <- a^2				
  List.res <- list(P, B, cp, u.dep, fact.imp, DOE);
  
  names(List.res) <- c("p", "indice.reliab", "Ncall", "Design.Point", "fact.imp", "DOE")
  
  #---------------------------Importance Sampling :---------------------------------#
  if(IS == TRUE){
    alpha  <- 0.05				
    N.IS   <- N.calls - cp			                 #Number of call reserved to the importance sampling
    norm.u <- u.dep%*%u.dep
    P.temp <- 0
    P 	 <- 0
    VAR 	 <- 0
    s	 <- 0
    Z 	 <- numeric(N.IS)
    
    v.temp <- matrix(rnorm(N.IS*dimension,0,1),ncol=dimension)
    v.temp <- v.temp  +  matrix(rep(u.dep,N.IS),ncol=dimension,byrow=TRUE)     #One simulate points around u.dep
    
    zz <- v.temp -  matrix(rep(u.dep,N.IS), ncol = dimension, byrow = TRUE)
    s  <- exp( u.dep%*%u.dep/2 - as.numeric(v.temp%*%u.dep) )
    
    test.v.temp <- lsf(t(v.temp))
    DOE$x <- cbind(DOE$x, t(v.temp))
    DOE$y <- c(DOE$y, test.v.temp)
    
    indic.Def <- 1*(test.v.temp <=0) 
    P.temp <- indic.Def*s
    
    List.res$Ncall  <- cp + N.IS					            #Numbers of call to the failure fonction
    List.res$p   <- P <- c(sum(P.temp)/N.IS)    			            #failure probability
    List.res$indice.reliab   <- -qnorm(P)  	        			            #indices of reliability
    List.res$variance <- VAR <- (sum(P.temp^2)/N.IS - P^2)/(N.IS-1)              #standard error of the estimator
    
    IC.inf     <- P - qnorm(1-alpha/2)*sqrt(VAR) 		#low bound of confidence interval
    IC.sup     <- P + qnorm(1-alpha/2)*sqrt(VAR)		#high bound of confidence interval
    
    List.res$Interval.conf <- c(IC.inf,IC.sup)
    List.res$fact.imp   <- a^2
    List.res$DOE <- DOE  
  }
  #----------------------------End Importance Sampling----------------------------#
  return(List.res)
}
