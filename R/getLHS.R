getLHS <- function(n, dimension, Q = 1e4, radius = qnorm(1e-5, lower.tail = FALSE)){
  lhs <- foreach(icount(Q)) %do% {
    lhsDesign(n, dimension)$design
  }
  crit <- sapply(lhs, function(l) min(stats::dist(l)))
  ind <- which.max(crit)
  lhs <- t(lhs[[ind]])
  rownames(lhs) <- rep(c("x", "y"), l = dimension)
  # 2*radius*lhs - radius
  qnorm(lhs)
}

# rewrite DiceDesign::lhsDesign because of random generation issue
lhsDesign <- function(n, dimension, randomized=TRUE, seed=NULL){
  # Randomized LHS: U[0,1]-sampling of n x dim values
  if (randomized) ran = matrix(runif(n*dimension),nrow=n,ncol=dimension)
  # Centered LHS
  else ran = matrix(0.5,nrow=n,ncol=dimension)

  x = matrix(0,nrow=n,ncol=dimension)  # initializing matrix x

  for (i in 1:dimension) {
    idx = sample(1:n)        # vector of permutations of [1 to n]
    P = (idx-ran[,i]) / n    # vecteur of probabilities
    x[,i] <- P  }

  # Outputs:
  return(list(n=n,dimension=dimension,design=x,randomized=randomized,seed=seed))
}
