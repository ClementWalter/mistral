## ----init, cache=FALSE, include=FALSE------------------------------------
library(knitr)
knitr::opts_chunk$set(fig.path='figure/',
               cache.path='cache/',
               cache=TRUE)

knitr::opts_knit$set(animation.fun = hook_scianimator)

## ----S2MART, fig.keep="all", fig.show="animate"--------------------------
smart <- mistral::S2MART(dimension = 2, lsf = mistral::waarts, failure = -2,
                        k1 = 3, k2 = 5, k3 = 6,
                        plot = TRUE)

## ----BMP-----------------------------------------------------------------
bmp <- mistral::BMP(dimension = 2, lsf = mistral::waarts, q = 0, N = 100,
                    N.iter = 0, X = akmcs$X, y = akmcs$y)

## ----BMP-learn, fig.keep='all', fig.show="animate"-----------------------
bmp <- mistral::BMP(dimension = 2, lsf = mistral::waarts, q = -4, N = 100,
                    N.iter = 2, plot = TRUE)

