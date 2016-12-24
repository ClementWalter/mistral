## ---- cache=FALSE, include=FALSE-----------------------------------------
library(knitr)
knitr::opts_chunk$set(fig.path='figure/',
               cache.path='cache/',
               cache=TRUE)

knitr::opts_knit$set(animation.fun = hook_scianimator)

## ------------------------------------------------------------------------
# return the number of cores of the computer
n <- parallel::detectCores()
# default behaviour if n not specified explained in the help page
cl <- parallel::makeCluster(1)
doSNOW::registerDoSNOW(cl)

# Control that everything is set properly
foreach::getDoParName()
foreach::getDoParWorkers()

## ------------------------------------------------------------------------
doMC::registerDoMC(1)
# Control that everything is set properly
foreach::getDoParName()
foreach::getDoParWorkers()

