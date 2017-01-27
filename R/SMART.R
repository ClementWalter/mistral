## -----------------------------------------------------------------------------
## Fonction SMART
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

#' @title Support-vector Margin Algoritm for Reliability esTimation
#' 
#' @description Calculate a failure probability with SMART method. This
#' should not be used by itself but only through S2MART.
#'
#' @author Clement WALTER \email{clementwalter@icloud.com}
#' 
#' @aliases SMART
#'
#' @details 
#' \code{SMART} is a reliability method proposed by J.-M. Bourinet et al. It makes
#' uses of a SVM-based metamodel to approximate the limit state function and calculates
#' the failure probability with a crude Monte-Carlo method using the metamodel-based
#' limit state function. As SVM is a classification method, it makes use of limit state
#' function values to create two classes : greater and lower than the failure threshold.
#' Then the border is taken as a surogate of the limit state function.
#' 
#' Concerning the refinement strategy, it distinguishes 3 stages, known as Localisation,
#' Stalibilsation and Convergence stages. The first one is proposed to reduce the margin
#' as much as possible, the second one focuses on switching points while the last one works
#' on the final Monte-Carlo population and is designed to insure a strong margin;
#' see F. Deheeger PhD thesis for more information.
#' 
#' @note
#' Problem is supposed to be defined in the standard space. If not, use \code{\link{UtoX}}
#' to do so.
#' 
#' Furthermore, each time a set of vector is defined as a matrix,
#' \sQuote{nrow} = \code{dimension} and \sQuote{ncol} = number of vector.
#' 
#' @references
#'   \itemize{
#' \item
#' J.-M. Bourinet, F. Deheeger, M. Lemaire:\cr
#' \emph{Assessing small failure probabilities by combined Subset Simulation and Support Vector Machines}\cr
#' Structural Safety (2011)
#' 
#' \item
#' F. Deheeger:\cr
#' \emph{Couplage mecano-fiabiliste : 2SMART - methodologie d'apprentissage stochastique en fiabilite}\cr
#'   PhD. Thesis, Universite Blaise Pascal - Clermont II, 2008
#' }
#' 
#' @seealso 
#' \code{\link{SubsetSimulation}}
#' \code{\link{MonteCarlo}}
#' \code{\link[e1071]{svm}} (in package \pkg{e1071})
#' \code{\link{S2MART}}
#' 
#' @export
SMART = function(dimension,
                 #' @param dimension the dimension of the input space
                 lsf,
                 #' @param lsf the limit-state function
                 N1 = 10000,
                 #' @param N1 Number of samples for the (L)ocalisation step
                 N2 = 50000,
                 #' @param N2 Number of samples for the (S)tabilisation step
                 N3 = 200000,
                 #' @param N3 Number of samples for the (C)onvergence step
                 Nu = 50,
                 #' @param Nu Size of the first Design of Experiments
                 lambda1 = 7,
                 #' @param lambda1 Relaxing parameter for MH algorithm at step L
                 lambda2 = 3.5,
                 #' @param lambda2 Relaxing parameter for MH algorithm at step S
                 lambda3 = 1,
                 #' @param lambda3 Relaxing parameter for MH algorithm at step C
                 tune_cost = c(1,10,100,1000),
                 #' @param tune_cost Input for tuning cost paramter of the SVM
                 tune_gamma = c(0.5,0.2,0.1,0.05,0.02,0.01),
                 #' @param tune_gamma Input for tuning gamma parameter of the SVM
                 clusterInMargin = TRUE,
                 #' @param clusterInMargin Enforce selected clusterised points to be in margin
                 alpha_margin = 1,
                 #' @param alpha_margin a real value defining the margin. While
                 #' 1 is the \sQuote{real} margin for a SVM, one can decide here to
                 #' stretch it a bit.
                 #Localization, Convergence and Stabilization bounds
                 k1 = round(6*(dimension/2)^(0.2)),
                 #' @param k1 Rank of the first iteration of step S
                 k2 = round(12*(dimension/2)^(0.2)),
                 #' @param k2 Rank of the first iteration of step C
                 k3 = k2 + 16,
                 #' @param k3 Rank of the last iteration of step C
                 #Arguments for SMART used in Subset simulation
                 X  = NULL,
                 #' @param X Coordinates of alredy known points
                 y = NULL,
                 #' @param y Value of the LSF on these points
                 failure   = 0,
                 #' @param failure Failure threshold
                 limit_fun_MH = NULL,
                 #' @param limit_fun_MH Define an area of exclusion with a limit function
                 sampling_strategy = "MH",
                 #' @param sampling_strategy Either MH for Metropolis-Hastings of AR for accept-reject
                 seeds = NULL,
                 #' @param seeds If some points are already known to be in the subdomain defined
                 #' by \code{limit_fun_MH}
                 seeds_eval = NULL,
                 #' @param seeds_eval Value of the metamodel on these points
                 burnin = 20,
                 #' @param burnin Burnin parameter for MH
                 thinning = 4,
                 #' @param thinning Thinning parameter for MH
                 plot = FALSE,
                 #' @param plot Set to TRUE for a full plot, ie. refresh at each iteration
                 limited_plot = FALSE,
                 #' @param limited_plot Set to TRUE for a final plot with final DOE, metamodel and LSF
                 add = FALSE,
                 #' @param add If plots are to be added to the current device
                 output_dir = NULL,
                 #' @param output_dir If plots are to be saved in jpeg in a given directory
                 z_MH = NULL,
                 #' @param z_MH For plots, if the limit_fun_MH has already been evaluated on the grid
                 z_lsf = NULL,
                 #' @param z_lsf For plots, if LSF has already been evaluated on the grid
                 verbose = 0) {
  #' @param verbose Either 0 for almost no output, 1 for medium size output and 2
  #' for all outputs
  
  
  cat("===================================\n")
  cat("   Beginning of SMART algorithm\n")
  cat("===================================\n")
  
  # Fix NOTE for R CMD check
  x <- z <- ..level.. <- NULL
  
  ## STEP 0 : INITIALISATION
  
  Ncall = 0;
  z_meta = NA;
  
  #Define the proportion of margin, switching and closest points at the beginning and at the end of each stage
  #current proportion at each iteration is given by linear regression between these two extremities
  proportion = list(L_stage=list(Nmargin=c(100,100),Nswitch=c(0,0),Nclose=c(0,0)),
                    S_stage=list(Nmargin=c(80,40),Nswitch=c(20,50),Nclose=c(0,10)),
                    C_stage=list(Nmargin=c(0,0),Nswitch=c(90,90),Nclose=c(10,10)))
  
  #Define the list variables used in the core loop
  #U stands for the points coordinates, one column per points, one row per dimension
  U = list(N1=matrix(nrow=dimension,ncol=N1),
           N2=matrix(nrow=dimension,ncol=N2),
           N3=matrix(nrow=dimension,ncol=N3),
           Nu=matrix(nrow=dimension,ncol=Nu),
           Nmargin=NA,
           Nswitch=NA,
           Nclose=NA)
  
  #G stands for the value of the limit state function on these points
  G = list(g=NA,#value on X points, ie all the value already calculated at a given iteration
           Nmargin=NA,
           Nswitch=NA,
           Nclose=NA,
           Nu=NA)
  #G_meta stands for the value of the surrogate on these points
  G_meta = list(N1=NA,
                N2=NA,
                N3=NA)
  
  # cat(" ============================================= \n")
  cat(" STEP 1 : EVALUATION OF A FIRST METAMODEL \n")
  # cat(" ======================================== \n")
  
  if(is.null(seeds) || sampling_strategy=="AR"){
    if(verbose>0){cat("* Get Ru maximum radius of N3 standard Gaussian samples \n")}
    tmp = matrix(rnorm(dimension*N3), nrow = dimension)
    radius = sqrt(colSums(tmp^2))
    Ru = max(radius)
    if(verbose>1){cat("   - Ru maximum radius of",N3,"standard samples =",Ru,"\n")}
  }
  
  if(is.null(limit_fun_MH)) {
    if(verbose>0){cat(" * Generate Nu =",Nu,"uniformly distributed samples in a sphere of radius Ru =",Ru,"\n")}
    U$Nu = t(runifSphere(dimension=dimension,N=Nu,radius=Ru))
  }
  else{
    if(sampling_strategy=="MH"){
      if(verbose>0){cat(" * Generate N1 =",N1,"points from seeds with l1r-mM algorithm\n")}
      gen_pop = generateWithlrmM(seeds=seeds$N1,
                                 seeds_eval=seeds_eval$N1,
                                 N=N1,
                                 lambda=lambda1,
                                 limit_f=limit_fun_MH,
                                 burnin=0,
                                 thinning=0)
      U$N1 = gen_pop$points
      G_meta$N1 = gen_pop$eval
      # U$N1 <- generateK(X = seeds$N1, N = N1, lsf = function(x) limit_fun_MH(x)$mean < 0)
      # G_meta$N1 <- limit_fun_MH(U$N1)$mean
      if(verbose>0){cat(" * Get Nu =",Nu,"points by clustering of the N1 points\n")}
      U$Nu = t(kmeans(t(U$N1), centers=Nu,iter.max=30)$centers)
    }
    else{
      if(verbose>0){cat(" * Generate Nu =",Nu,"uniformly distributed samples with an accept-reject strategy\n")}
      rand = function(dimension,N) {t(runifSphere(dimension,N,radius=Ru))}
      U$Nu = generateWithAR(dimension=dimension,N=Nu,limit_f=limit_fun_MH,rand=rand)
    }
  }
  rownames(U$Nu) <- rep(c('x', 'y'), length.out = dimension)
  
  if(verbose>0){cat(" * Assessment of g on these points\n")}
  G$Nu = lsf(U$Nu);Ncall = Ncall + Nu
  
  #plotting part
  if(plot==TRUE){
    
    xplot <- yplot <- c(-80:80)/10
    df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
    
    if(add==TRUE){
      p <- last_plot()
    } else {
      p <- ggplot(data = df_plot, aes(x,y)) +
        geom_contour(aes(z=z, color=..level..), breaks = failure) +
        theme(legend.position = "none") +
        xlim(-8, 8) + ylim(-8, 8)
    }
    if(verbose>0){cat(" * 2D PLOT : FIRST DoE \n")}
    if(is.null(limit_fun_MH)) {
      circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
        r = diameter / 2
        tt <- seq(0,2*pi,length.out = npoints)
        xx <- center[1] + r * cos(tt)
        yy <- center[2] + r * sin(tt)
        return(data.frame(x = xx, y = yy))
      }
      
      dat <- circleFun(diameter = Ru*2,npoints = 100)
      #geom_path will do open circles, geom_polygon will do filled circles
      p <- p + geom_path(data = dat, aes(alpha = 0.3))
      print(p)
    }
    else{
      z_MH = limit_fun_MH(t(df_plot[,1:2]))
      df_plot_MH <- data.frame(df_plot[,1:2], z = z_MH$mean)
      p <- p + geom_contour(data = df_plot_MH, aes(z=z, color=..level..), breaks = 0)
    }
    p <- p + geom_point(data = data.frame(t(U$Nu), z = lsf(U$Nu)), aes(color = z))
    print(p)
  }
  
  #add points to the learning database
  if(verbose>0){cat(" * Add points to the learning database\n")}
  if(is.null(X)){
    X = cbind(seq(0,0,l=dimension),U$Nu); dimnames(X) <- NULL;
    dimnames(X) <- NULL
    g0 = lsf(as.matrix(seq(0,0,l=dimension)));Ncall = Ncall + 1;#this is to insure consistency between model sign and lsf sign
    G$g = c(g0,G$Nu)
  }
  else{
    X = cbind(X,U$Nu); dimnames(X) <- NULL;
    dimnames(X) <- NULL
    G$g = c(y,G$Nu)
  }
  
  #Train the model
  if(verbose>0){cat(" * Train the model\n")}
  if(verbose<1){capture.output(meta <- trainModel(design=X,response=(G$g-failure),type="SVM",cost=tune_cost,gamma=tune_gamma))}
  else{meta = trainModel(design=X,response=(G$g-failure),type="SVM",cost=tune_cost,gamma=tune_gamma)}
  
  #Update meta_fun & meta_model
  if(verbose>0){cat(" * UPDATE quantities based on surrogate model \n")}
  meta_model = meta$model
  meta_fun = meta$fun
  
  #plotting part
  if(plot==TRUE){
    if(verbose>0){cat(" * 2D PLOT : FIRST METAMODEL\n\n")}
    z_meta = meta_fun(t(df_plot[,1:2]))
    df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean)
    print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = c(0, -1, 1)))
  }
  
  # cat(" ============================= \n")
  cat(" STEP 2 : REFINEMENT PROCEDURE \n")
  # cat(" ============================= \n\n")
  
  for (k in c(1:k3)) {
    stage = 1*(k<k1) + 2*(k>=k1)*(k<k2) + 3*(k>=k2)
    stage_str = switch(stage, "1" = {"Localisation"}, "2" = {"Stabilisation"}, "3" = {"Convergence"})
    
    cat(" * ITERATION",k,"of", k3,":",stage_str,"stage \n")
    # cat("    ======================================================== \n\n")
    
    #Get the parameters corresponding to the current stage
    if(verbose>0){cat(" * Get the parameters corresponding to the current stage\n")}
    k_start = 1*(stage==1) + k1*(stage==2) + k2*(stage==3)
    k_end = (k1-1)*(stage==1) + (k2-1)*(stage==2) + k3*(stage==3)
    N = N1*(stage==1) + N2*(stage==2) + N3*(stage==3)
    if(verbose>1){cat("   - Number of points to sample =",N,"\n")}
    lambda = lambda1*(stage==1) + lambda2*(stage==2) + 1*(stage==3)
    
    #Get the parameters corresponding to the current iteration
    if(verbose>0){cat(" * Get the parameters corresponding to the current iteration\n")}
    Nsup = round((3+(k-1)/(k1-1))*sqrt(dimension))*(stage==1) +
      round((4+(k-k1)/(k2-k1))*sqrt(dimension))*(stage==2) +
      round(5*sqrt(dimension))*(stage==3)
    aMargin = proportion[[stage]]$Nmargin[1]; bMargin = proportion[[stage]]$Nmargin[2]
    aSwitch = proportion[[stage]]$Nswitch[1]; bSwitch = proportion[[stage]]$Nswitch[2]
    aClose = proportion[[stage]]$Nclose[1]; bClose = proportion[[stage]]$Nclose[2]
    Nmargin = round(Nsup/100*(aMargin+(k-k_start)/(k_end-k_start)*(bMargin-aMargin)))
    Nswitch = round(Nsup/100*(aSwitch+(k-k_start)/(k_end-k_start)*(bSwitch-aSwitch)))
    Nclose = round(Nsup/100*(aClose+(k-k_start)/(k_end-k_start)*(bClose-aClose)))
    if(verbose>1){cat("   - Nsup =",Nsup,"\n")
      cat("   - Nmargin =",Nmargin,"\n")
      cat("   - Nswitch =",Nswitch,"\n")
      cat("   - Nclose =",Nclose,"\n")}
    
    if(is.null(limit_fun_MH)){
      if(stage==3) {
        #Generate standard gaussian samples
        if(verbose>0){cat(" * Generate standard gaussian samples\n")}
        U[[stage]] = matrix(rnorm(dimension*N, mean=0, sd=1),dimension,N)
      }
      else {
        #Generate samples with a uniform distribution on a sphere of radius Ru
        if(verbose>0){cat(" * Generate samples with a uniform distribution on a sphere of radius Ru =",Ru,"\n")}
        U[[stage]] = t(runifSphere(dimension=dimension,N=N,radius=Ru))
      }
    }
    else{
      if(sampling_strategy=="MH"){
        if(verbose>0){cat(" * Generate N =",N,"points with lr-mM algorithm\n")}
        U[[stage]] = generateWithlrmM(seeds=seeds[[stage]],
                                      seeds_eval=seeds_eval[[stage]],
                                      N=N,
                                      lambda=lambda,
                                      limit_f=limit_fun_MH,
                                      burnin=0,thinning=0)$points
        # U[[stage]] = generateK(X = seeds[[stage]], N = N, lsf = function(x) limit_fun_MH(x)$mean < 0)
      }
      else{
        if(verbose>0){cat(" * Generate N =",N,"points with an accept-reject strategy\n")}
        if(stage==3) {
          rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}
          U[[stage]] = generateWithAR(dimension=dimension,N=N,limit_f=limit_fun_MH,rand=rand)
        }
        else {
          rand = function(dimension,N) {t(runifSphere(dimension=dimension,N=N,radius=Ru))}
          U[[stage]] = generateWithAR(dimension=dimension,N=N,limit_f=limit_fun_MH,rand=rand)
        }
      }
    }
    
    #Evaluate g on U[[stage]] using the meta model
    if(verbose>0){cat(" * Evaluate the metamodel on these points \n")}
    meta_pred = meta_fun(U[[stage]])
    G_meta[[stage]] = meta_pred$mean
    
    #Selection of Nsup = Nmargin + Nswitch + Nclose new points where g is to be assessed
    if(Nmargin>0) {
      #Selection of Nmargin new points where g is to be assessed
      if(verbose>0){cat(" * Selection of Nmargin =",Nmargin,"new points where g is to be assessed\n")}
      isMargin = inMargin(meta_pred,type="SVM",alpha=alpha_margin)
      if(verbose>1){cat("   - Number of margin points =",sum(isMargin*1),"\n")}
      if(!sum(isMargin)==0){
        U$Nmargin = tryCatch(
          clusterize(data=U[[stage]][,isMargin],Ncluster=Nmargin,inMargin=clusterInMargin),
          error = function(cond) {
            message(cond)
            message("\nall margin points are kept")
            return(U[[stage]][,isMargin])
          })
        U$Nmargin = as.matrix(U$Nmargin)
        rownames(U$Nmargin) <- rep(c('x', 'y'), length.out = dimension)
        
        #Assessment of g
        if(verbose>0){cat(" * Assessment of g on these points \n")}
        G$Nmargin <- lsf(U$Nmargin);Ncall = Ncall + Nmargin
        #Add points U$Nmargin to the learning database
        if(verbose>0){cat(" * Add points to the learning database\n")}
        X <- cbind(X,U$Nmargin); dimnames(X) <- NULL;
        G$g = c(G$g,G$Nmargin)
        
        #plotting part
        if(plot==TRUE){
          if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
          print(p_meta <- p_meta + geom_point(data = data.frame(t(U$Nmargin), z = G$Nmargin), aes(color = z)))
          p <- p + geom_point(data = data.frame(t(U$Nmargin), z = G$Nmargin), aes(color = z))
        }
        
        #Train the model
        if(verbose>0){cat(" * Train the model\n")}
        meta = trainModel(meta_model,design=X,response=(G$g-failure),updesign=U$Nmargin,upresponse=(G$Nmargin-failure),type="SVM")
        #Update meta_fun & meta_model
        meta_model.prev = meta_model
        meta_model = meta$model
        meta_fun.prev = meta_fun
        meta_fun = meta$fun
      }
    }
    
    if(Nswitch>0) {
      #Selection of Nswitch new points where g is to be assessed
      if(verbose>0){cat(" * Selection of Nswitch =",Nswitch,"new points where g is to be assessed\n")}
      G_meta_prev = meta_fun.prev(U[[stage]])$mean
      isSwitched = as.logical(sign(G_meta_prev*G_meta[[stage]])<0)
      if(verbose>1){cat("   - Number of switching points =",sum(isSwitched*1),"\n")}
      if(!sum(isSwitched)==0){
        U$Nswitch <- tryCatch(t(kmeans(t(U[[stage]][,isSwitched]), centers=Nswitch, iter.max=30)$centers),
                              error = function(cond) {
                                message(cond);
                                message("\nall switching points are kept");
                                return(U[[stage]][,isSwitched]);
                              })
        U$Nswitch = as.matrix(U$Nswitch)
        rownames(U$Nswitch) <- rep(c('x', 'y'), length.out = dimension)
        #Assessment of g
        if(verbose>0){cat(" * Assessment of g\n")}
        G$Nswitch = lsf(U$Nswitch);Ncall = Ncall + Nswitch
        #Add points U$Nmargin + U$Nswitch to the learning database
        if(verbose>0){cat(" * Add points U$Nswitch  to the learning database\n")}
        X = cbind(X,U$Nswitch); dimnames(X) <- NULL;
        G$g = c(G$g,G$Nswitch)
        
        if(plot==TRUE){
          if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
          print(p_meta <- p_meta + geom_point(data = data.frame(t(U$Nswitch), z = G$Nswitch), aes(color = z)))
          p <- p + geom_point(data = data.frame(t(U$Nswitch), z = G$Nswitch), aes(color = z))
        }
        
        #Train the model
        if(verbose>0){cat(" * Train the model\n")}
        meta = trainModel(meta_model,design=X,response=(G$g-failure),updesign=U$Nswitch,upresponse=(G$Nswitch-failure),type="SVM")
        #Update meta_fun & meta_model
        meta_model.prev = meta_model
        meta_model = meta$model
        meta_fun.prev = meta_fun
        meta_fun = meta$fun
      }
    }
    
    if(Nclose>0) {
      #Selection of Nclose new points where g is to be assessed
      if(verbose>0){cat(" * Selection of Nclose =",Nclose,"new points where g is to be assessed\n")}
      isClose = sort(abs(G_meta[[stage]]),index=TRUE)$ix[1:Nclose]
      U$Nclose = as.matrix(U[[stage]][,isClose])
      rownames(U$Nclose) <- rep(c('x', 'y'), length.out = dimension)
      
      #Assessment of g
      if(verbose>0){cat(" * Assessment of g\n")}
      G$Nclose = lsf(U$Nclose);Ncall = Ncall + Nclose
      #Add points U$Nclose to the learning database
      if(verbose>0){cat(" * Add points U$Nclose to the learning database\n")}
      X = cbind(X,U$Nclose); dimnames(X) <- NULL;
      G$g = c(G$g,G$Nclose)
      
      if(plot==TRUE){
        if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
        print(p_meta <- p_meta + geom_point(data = data.frame(t(U$Nclose), z = G$Nclose), aes(color = z)))
        p <- p + geom_point(data = data.frame(t(U$Nclose), z = G$Nclose), aes(color = z))
      }
      
      #Train the model
      if(verbose>0){cat(" * Train the model\n")}
      meta = trainModel(meta_model,design=X,response=(G$g-failure),updesign=U$Nclose,upresponse=(G$Nclose-failure),type="SVM")
      #Update meta_fun & meta_model
      meta_model.prev = meta_model
      meta_model = meta$model
      meta_fun.prev = meta_fun
      meta_fun = meta$fun
    }
    
    #plotting part
    if(plot==TRUE){
      if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
      z_meta = meta_fun(t(df_plot[,1:2]))
      df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean)
      print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = c(0, -1, 1)))
    }
  }
  
  if( (plot || limited_plot) & add==FALSE & !is.null(output_dir)) {dev.off()}
  
  # cat("\n ======================================= \n")
  cat(" STEP 3 : FAILURE PROBABILITY ESTIMATION \n")
  # cat(" ======================================= \n\n")
  
  #estimation of probability P using metamodel classifier values
  # if(verbose>0){cat(" * estimation of probability P using metamodel classifier values\n")}
  if(!is.null(limit_fun_MH)) {
    # capture.output(
    if(verbose>0){cat(" * Generate N =",N,"points with lr-mM algorithm\n")}
    gen_pop <- generateWithlrmM(seeds=seeds$N3,
                                seeds_eval=seeds_eval$N3,
                                N=N,
                                lambda=lambda,
                                limit_f=limit_fun_MH,
                                burnin=burnin,thinning=thinning,
                                VA_function=function(x){1*(meta_fun(x)$mean<0)})$points
    # gen_pop <- generateK(X = seeds$N3, N = N, lsf = function(x) limit_fun_MH(x)$mean < 0)
    # )
    meta_eval <- limit_fun_MH(gen_pop)$mean
    fail_points = (meta_eval<0)
    points=gen_pop[,fail_points]
    # meta_eval= meta_fun(points)$mean
    # MC_est = gen_pop$est
    MC_est <- mean(fail_points)
  }
  else{
    G_meta$N3 = meta_fun(U$N3)$mean
    fail_points = (G_meta$N3<0)
    points=U$N3[,fail_points]
    meta_eval=G_meta$N3[fail_points]
    MC_est = sum(1*fail_points)/N3
  }
  MC_var = MC_est*(1-MC_est)/N3
  MC_delta <- sqrt(MC_var)/MC_est
  
  G_meta$N1 = meta_fun(U$N1)$mean
  indN1 = (G_meta$N1<0)
  G_meta$N2 = meta_fun(U$N2)$mean
  indN2 = (G_meta$N2<0)
  points = list(N1 = U$N1[,indN1],
                N2 = U$N2[,indN2],
                N3 = points)
  meta_eval = list(N1=G_meta$N1[indN1],N2=G_meta$N2[indN2],N3=meta_eval)
  
  cat("===================================\n")
  cat("      End of SMART algorithm\n")
  cat("===================================\n")
  
  cat(" * proba =",MC_est,"\n")
  # cat(" * variance =",MC_var,"\n")
  cat(" * cov =",MC_delta,"\n")
  # cat(" * Markov chain gamma =",MC_gamma,"\n")
  
  #' @return
  #'   An object of class \code{list} containing the failure probability and some more outputs as described below:
  
  res = list(proba = MC_est,
             #' \item{proba}{The estimated failure probability.}
             cov = MC_delta,
             #' \item{cov}{The coefficient of variation of the Monte-Carlo probability estimate.}
             Ncall = Ncall,
             #' \item{Ncall}{The total number of calls to the \code{limit_state_function}.}
             X = X,
             #' \item{X}{The final learning database, ie. all points where \code{lsf} has been calculated.}
             y = G$g,
             #' \item{y}{The value of the \code{limit_state_function} on the learning database.}
             meta_fun = meta_fun,
             #' \item{meta_fun}{The metamodel approximation of the \code{limit_state_function}.
             #' A call output is a list containing the value and the standard deviation.}
             meta_model = meta_model,
             #' \item{meta_model}{The final metamodel.}
             points = points,
             #' \item{points}{Points in the failure domain according to the metamodel.}
             meta_eval = meta_eval)
  #'   \item{meta_eval}{Evaluation of the metamodel on these points.}
  
  if(plot+limited_plot) {
    res = c(res,list(z_meta=z_meta$mean))
    #'   \item{z_meta}{If \code{plot}==TRUE, the evaluation of the metamodel on the plot grid.}
  }
  
  return(res)
}
