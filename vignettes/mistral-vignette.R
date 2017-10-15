## ----init, cache=FALSE, include=FALSE------------------------------------
library(knitr)
knitr::opts_chunk$set(fig.path='figure/',
               cache.path='cache/',
               cache=TRUE)

knitr::opts_knit$set(animation.fun = hook_scianimator)

## ----input-margin--------------------------------------------------------
distX1 <- list(type='Norm',  MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X1')
distX2 <- list(type='Norm',  MEAN=0.0, STD=1.0, P1=NULL, P2=NULL, NAME='X2')
input.margin <- list(distX1,distX2)

## ----modifCorrMatrix-----------------------------------------------------
input.Rho    <- matrix( c(1.0, 0.5,
                          0.5, 1.0),nrow=2)
input.R0     <- mistral::ModifCorrMatrix(input.Rho)
L0           <- t(chol(input.R0))

## ----UtoX----------------------------------------------------------------
U <- rnorm(2)
U <- cbind(U, U)
X <- mistral::UtoX(U, input.margin, L0)
X

## ----lsf-UtoX------------------------------------------------------------
lsf_U = function(U) {   
    X <- mistral::UtoX(U, input.margin, L0)
    lsf(X)
}

## ----lsf-vect------------------------------------------------------------
lsf <- function(x){
  x[1] + x[2]
}

## ----lsf-vect-example----------------------------------------------------
x <- c(1,2)
lsf(x)

## ----lsf-vect-example-matrix---------------------------------------------
x <- as.matrix(x)
lsf(x)

## ----lsf-vect-fail-matrix------------------------------------------------
x <- cbind(x,x)
lsf(x)

## ----lsf-matrix----------------------------------------------------------
lsf <- function(x) {
  x[1,] + x[2,]
}
lsf(x)

## ----lsf-matrix-x--------------------------------------------------------
x <- 1:2
as.matrix(x)

## ----lsf-matrix-error, error=TRUE----------------------------------------
x <- cbind(x,x)
lsf(x[,1])

## ----lsf-asmatrix--------------------------------------------------------
lsf <- function(x) {
  x <- as.matrix(x)
  x[1,] + x[2,]
}
lsf(x[,1])

## ----lsf-apply-----------------------------------------------------------
lsf <- function(x){
  x <- as.matrix(x)
  apply(x, 2, myCode)
}

## ----lsf-foreach---------------------------------------------------------
require(foreach)
lsf <- function(x){
  x <- as.matrix(x)
  foreach(x = iterators::iter(x, by = 'col'), .combine = 'c') %dopar% {
    myCode(x)
  }
}

## ----foreach-nobackend, warning=TRUE-------------------------------------
myCode <- function(x) x[1] + x[2]
x <- 1:2
lsf(x)

## ----foreach-registerDoSeq-----------------------------------------------
foreach::registerDoSEQ()

## ----foreach-registerDoSNOW----------------------------------------------
# return the number of cores of the computer
n <- parallel::detectCores()
# default behaviour if n not specified explained in the help page
cl <- parallel::makeCluster(1)
doSNOW::registerDoSNOW(cl)

# Control that everything is set properly
foreach::getDoParName()
foreach::getDoParWorkers()

## ----foreach-stopCluster-------------------------------------------------
parallel::stopCluster(cl)

## ----basic-MC------------------------------------------------------------
X <- matrix(rnorm(2e5), nrow = 2) # generate 1e5 standard Gaussian samples
Y <- mistral::kiureghian(X) # evaluate to model to get 1e5 iid samples
q <- 0 # define the threshold
(p <- mean(Y<q)) # estimate P[g(X)<0]

## ----mistral-MC----------------------------------------------------------
mc <- mistral::MonteCarlo(dimension = 2, lsf = mistral::kiureghian, N_max = 1e5, q = q,
                          # these first parameters are exactly the one used above
                           N_batch = 1e4) # define the batch size

## ----mistral-MC-nolimit--------------------------------------------------
mc <- mistral::MonteCarlo(dimension = 2, lsf = mistral::kiureghian, N_max = Inf, q = q,
                           N_batch = 1e4, # define the batch size
                           verbose = 1) # control the level of log messages

## ----mistral-MC-plot-----------------------------------------------------
mc <- mistral::MonteCarlo(dimension = 2, lsf = mistral::kiureghian, N_max = 1e4, q = q,
                           N_batch = 1e4, # define the batch size
                           plot = TRUE)

## ----MC-ecdf-plot--------------------------------------------------------
require(ggplot2)
y = seq(-5, 10, l = 200)
ggplot(data.frame(y=y, ecdf = mc$ecdf(y)), aes(y,ecdf)) + geom_line()

## ----SubsetSimulation, fig.keep='all', fig.show='animate'----------------
ss <- mistral::SubsetSimulation(dimension = 2, lsf = mistral::kiureghian, q = 0,
                                N = 1e4,
                                plot = TRUE)

## ----MP------------------------------------------------------------------
foreach::registerDoSEQ()
mp <- mistral::MP(dimension = 2,  lsf = mistral::kiureghian, q = 0, N = 1e2)

## ----SS-MP-benchmark-----------------------------------------------------
ss$cv^2*ss$Ncall / (mp$cv^2*sum(mp$Ncall))

## ----MP-ecdf-------------------------------------------------------------
y = seq(0, 10, l = 200)
ggplot(data.frame(y=y, ecdf = mp$ecdf(y)), aes(y,ecdf)) + geom_line()

## ----MP-quantile---------------------------------------------------------
mp <- mistral::MP(dimension = 2,  lsf = mistral::kiureghian, p = mp$p, N = 1e2)

## ----MP-quantile-confidence----------------------------------------------
mp <- mistral::MP(dimension = 2,  lsf = mistral::kiureghian, p = mp$p, N = 1e2,
                  compute_confidence = TRUE)

## ----SVM-tuto------------------------------------------------------------
require(e1071)
X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
Y <- rowSums(X^2)
(svm.model <- svm(X, (Y>1), type = "C-classification"))
X.test <- data.frame(x1=rnorm(1), x2=rnorm(1))
predict(svm.model, X.test)
sum(X.test^2)

## ----kriging-tuto--------------------------------------------------------
require(DiceKriging)
X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
Y <- rowSums(X^2)
km.model <- km(design = X, response = Y)
x.new <- data.frame(x1=rnorm(1), x2=rnorm(1))
print(sum(x.new^2))
predict(km.model, x.new, type = "UK")[c('mean', 'sd')]

## ----FORM----------------------------------------------------------------
form <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000, u.dep = c(0,0))
form$p

## ----FORM-IS-------------------------------------------------------------
form.IS <- mistral::FORM(dimension = 2, mistral::kiureghian, N.calls = 1000,
                         u.dep = c(0,0),
                         IS = TRUE)
form.IS$p

## ----MetaIS, fig.keep="all", fig.show="animate"--------------------------
metais <- mistral::MetaIS(dimension = 2, lsf = mistral::waarts, N = 3e5, K_alphaLOO = 5,
                          plot = TRUE)

## ----AKMCS, fig.keep="all", fig.show="animate"---------------------------
akmcs <- mistral::AKMCS(dimension = 2, lsf = mistral::waarts, N = 3e5, plot = TRUE, Nmax = 10)

## ----BMP-----------------------------------------------------------------
bmp <- mistral::BMP(dimension = 2, lsf = mistral::waarts, q = 0, N = 100,
                    N.iter = 0, X = akmcs$X, y = akmcs$y)

## ----BMP-learn, fig.keep='all', fig.show="animate"-----------------------
bmp <- mistral::BMP(dimension = 2, lsf = mistral::waarts, q = -4, N = 100,
                    N.iter = 2, plot = TRUE)

