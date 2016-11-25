## -----------------------------------------------------------------------------
## Lancement du fichier : executer la commande
## R --no-save < *.R
## -----------------------------------------------------------------------------
##    Copyright (C) 2016
##    Gilles DEFAUX
##    CEA / DIF / DCSA / BSE
##    gilles.defaux@cea.fr
## -----------------------------------------------------------------------------

## -----------------------------------
## PARAMETRES ET OPTIONS
## -----------------------------------
set.seed(123456)
DT     = 0.1
Tfin   = 40

NbTsample = as.integer(Tfin/DT)
times = seq(0, Tfin, by=DT)

params = c(sigma=3., rho=26., beta=1., alpha=1.)
state0 = c(X=5.5, Y=5.5, Z=25.5)
rps    = params['rho'] + params['sigma']

## -----------------------------------
## MODEL DETERMINISTE
## -----------------------------------

Lorentz <- local({
  DT     = DT
  Tfin   = Tfin
  NbTsample = NbTsample
  times = times
  params = params
  state0 = state0
  rps    = rps
  function(t, state, params) {
    with(as.list(c(state,params)), {
      dX <- sigma*(Y-X)
      dY <- X*(rho-Z) - Y
      dZ <- X*Y - beta*Z
      list(c(dX,dY,dZ))
    })
  }
})

LorentzWithBM <-local({
  DT     = DT
  Tfin   = Tfin
  NbTsample = NbTsample
  times = times
  params = params
  state0 = state0
  rps    = rps
  function(t, state, params, ZZ) {
    with(as.list(c(state,params)), {
      k  <- as.integer(t/DT)
      ZZ[1] = 0.
      Uk <- alpha*sqrt(DT)*sum(ZZ[1:k])
      dX <- sigma*(Y-X) + Uk
      dY <- X*(rho-Z) - Y
      dZ <- X*Y -beta*Z
      list(c(dX,dY,dZ))
    })
  }
})

Gt <- local({
  DT     = DT
  Tfin   = Tfin
  NbTsample = NbTsample
  times = times
  params = params
  state0 = state0
  rps    = rps
  function(out) {
    temp = out[,'X']^2/(rps^2*params['beta']/params['sigma']) +
      out[,'Y']^2/(rps^2*params['beta']) +
      (out[,'Z']-rps)^2/(rps^2)
    return(temp)
  }
})
##
## SANS EXCITATION
##

out  <- deSolve::ode(y = state0, times=times, func = Lorentz, parms=params)

scatterplot3d::scatterplot3d(x=out[,'X'], y=out[,'Y'], z=out[,'Z'], type='l',
              highlight.3d=TRUE, col.grid="lightblue", box=FALSE,
              xlab='X', ylab='Y', zlab='Z')

plot(times,Gt(out),type='l',ylim=c(0.,1.2))

rm(out)

##
## AVEC EXCITATION PAR MOUVEMENT BROWNIEN
##

out2 <- deSolve::ode(y = state0, times=times, func = LorentzWithBM, parms=params, ZZ=rnorm(NbTsample))

plot(times,Gt(out2),type='l',ylim=c(0.,1.2))

rm(out2)


## -----------------------------------
## FONCTION DE PERFORMANCE
## -----------------------------------

PlotTraj = 0
Debug    = 0

## VERSION VECTORIELLE
myCode <- local({
  DT     = DT
  Tfin   = Tfin
  NbTsample = NbTsample
  times = times
  params = params
  state0 = state0
  rps    = rps
#   PlotTraj = PlotTraj
#   Debug = Debug
  LorentzWithBM = LorentzWithBM
  Lorentz = Lorentz

  function(ZZ) {
    out2 <- deSolve::ode(y = state0, times=times, func = LorentzWithBM, parms=params, ZZ=ZZ)
    # if(PlotTraj == 1) { lines(times, Gt(out2), type='l', lty=2, ylim=c(0.,1.2)) }
    resu <- 1.0 - max(Gt(out2))
    # if(Debug == 1) { print(resu) }
    return(resu)
  }
})

## VERSION MATRICIELLE
myLSF <- local({
  myCode = myCode
  function(UU) {
    UU <- as.matrix(UU)
    apply(UU,2,myCode)
  }
})

## VERSION PARRALLELE
myParLSF = function(UU) {
  UU <- as.matrix(UU)
  Resu <- foreach( u=iter(UU, by='col'), .combine = c) %dopar% {myCode(u)}
  if(Debug==1) { print(resu) }
  return(Resu)
}

# 
# ## TEST
# Do.Test.PAR = FALSE
# Do.Test     = FALSE
# 
# set.seed(123456)
# NbTest   = 10
# 
# if(Do.Test.PAR == TRUE){
#   PlotTraj = 0
#   U0 = matrix(rnorm(NbTest*NbTsample), nrow = NbTsample, ncol = NbTest)
#   test = myParLSF(U0)
#   print(test)
# }
# 
# if(Do.Test == TRUE){
#   PlotTraj = 1
#   U0 = matrix(rnorm(NbTest*NbTsample), nrow = NbTsample, ncol = NbTest)
#   test = myLSF(U0)
#   print(test)
# }


## -----------------------------------
## CALCUL
## -----------------------------------

PlotTraj = 0
Debug    = 0

set.seed(123456)

require(doParallel)
require(foreach)
registerDoParallel(cores = parallel::detectCores())

res <- list()
for(i in 1:1){
  res[[i]] = mistral::MP(dimension = NbTsample,
                         lsf = myLSF,
                         q = 0.0,
                         lower.tail = TRUE,
                         N = 1800,
                         plot = FALSE,
                         verbose = 2)
}
# Close backend
# stopCluster(cl)


## -----------------------------------
## RESULTATS
## * REFERENCE :
## -----------------------------------

# cat('PROBABILITE DE DEFAILLANCE = ',res$p)
# cat('COEFFICIENT DE VARIATION   = ',res$cov)

save(res, file='Resu-Lorentz-MP')


