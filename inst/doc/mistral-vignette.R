## ----init, cache=FALSE, include=FALSE------------------------------------
library(knitr)
knitr::opts_chunk$set(fig.path='figure/',
               cache.path='cache/',
               cache=TRUE)

knitr::opts_knit$set(animation.fun = hook_scianimator)

## ----MP------------------------------------------------------------------
foreach::registerDoSEQ()
mp <- mistral::MP(dimension = 2,  lsf = mistral::kiureghian, q = 0, N = 1e2)

