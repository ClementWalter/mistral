## -----------------------------------------------------------------------------
## Fonction generateWithMH
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

generateWithlrmM = function(seeds,
                            seeds_eval,
                            N,
                            lambda = 1,
                            limit_f,
                            burnin=20,
                            thinning=4,
                            modified = TRUE,
                            VA_function = NULL,
                            VA_esp = NULL,
                            VA_var = NULL) {

	#VA_function is an optional function to add if samples are to be used for a MC estimator as sum(VA_function(x)) to get correlation between samples
	#VA_esp & VA_var are the esperance and the variance of the random variable VA_function(x)

	tic = proc.time()[3]
	
	# Fix foreach NOTE
	j <- NULL

	# select only non duplicated seeds
	seeds = as.matrix(seeds)
	dup <- duplicated(t(seeds))
	seeds <- as.matrix(seeds[,!dup])
	d = dim(seeds)
	n_seeds = d[2];
	d = d[1]
	seeds_eval <- ifelse(missing(seeds_eval), rep(-1, n_seeds), seeds_eval[!dup])
	
	#set variables
	chain_length = ceiling(N/n_seeds)
	U = matrix(NA,d,chain_length*n_seeds)
	G = 1:chain_length*n_seeds
	Ncall = 0;

	if(missing(limit_f)) {
		limit_fun = function(x) {-1}
	}
	else {
		limit_fun = function(x) {tryCatch(limit_f(x)$mean, error = function(cond) {return(limit_f(x))})}
	}

	#core loop
	# for (j in 1:n_seeds){
	 MH <- foreach(j=1:n_seeds, .combine = cbind) %dopar% {
		# start = (j-1)*chain_length+1
		# end = j*chain_length
		seed = seeds[,j]
		eval_seed = seeds_eval[j]
		MH = MetropolisHastings(x0 = seed,
		                        eval_x0 = eval_seed,
		                        chain_length = (chain_length-1),
		                        modified = modified,
		                        lambda = lambda,
		                        limit_fun = limit_fun,
		                        burnin = burnin,
		                        thinning = thinning)
		
		length(MH$Ncall) <- length(MH$eval)
    rbind(MH$points, MH$eval, MH$Ncall)
		# U[,start:end] = MH$points
		# G[start:end] = MH$eval
		# Ncall = Ncall + MH$Ncall
	}
	U <- MH[1:d,]
	G <- MH[d+1,]
	Ncall <- sum(MH[d+2,], na.rm = TRUE)
	 
	if(dim(U)[2]>N) {
	  sel <- sample(1:dim(U)[2], size = N, replace = FALSE)
	  U = U[,sel]; G = G[sel]
	}
	dimnames(U) <- NULL
	toc = proc.time()[3]-tic
	cat("   -",(burnin + (thinning+1)*chain_length)*n_seeds,"points generated in",toc,"sec. with ",n_seeds,"seeds,",N,"points kept : burnin =",burnin,"thinning =",thinning,"\n")
	cat("   -",sum(duplicated(t(U))), "duplicated samples \n")

	res = list(points=U,
	           eval=G,
	           chain_length=chain_length,
	           Ncall=Ncall)
	
	if(!is.null(VA_function)) {
		# cat("====== Beginning of Monte-Carlo estimation ======\n")
		cat(" * Calculate Monte-Carlo estimator\n")
		VA_values = VA_function(U)
		MC_est = mean(VA_values)
		cat("   - MC_est =",MC_est,"\n")

		if(MC_est==0){
		  print(dim(U))
		}
		
		# cat(" * Calculate the covariance between samples \n")
		# stat = MCMCcovariance(n_seeds=n_seeds,chain_length=chain_length,VA_values=VA_values,VA_esp=MC_est)
		# cat('   - gamma =', stat$gamma)
		# cat('   - cov =', stat$cov)
		# MC_gamma = stat$gamma
		# MC_var = stat$var
		# MC_delta = stat$cov

		res = c(res, list(VA_values = VA_values,
                       est = MC_est
                       # var = MC_var,
                       # delta = MC_delta,
                       # gamma = MC_gamma
		                  ))
	}
	return(res)
}