# salmonIPM

This is the development repo for **salmonIPM**, an R package for fitting integrated population models to salmon data.

## Installation

1. Install and configure **rstan** by following the instructions in the [RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) vignette. In particular, **salmonIPM** requires **rstan** version 2.26 or higher, which (as of 2023-05-10) is only available from [https://mc-stan.org/r-packages/](https://mc-stan.org/r-packages/). If you have previously installed earlier versions of **rstan** and **StanHeaders** from CRAN, you will need to uninstall and update them as shown in the vignette.

2. Install the current version of **salmonIPM** from GitHub using **devtools**. Because the repo is private for the time being, it is necessary to [generate a personal access token](https://github.com/settings/tokens) (PAT) and pass it to `install_github()` as discussed [here](https://stackoverflow.com/questions/21171142/how-to-install-r-package-from-private-repo-using-devtools-install-github).

```r
if(!require("devtools")) {
  install.packages("devtools")
  library("devtools")
}
devtools::install_github("ebuhle/salmonIPM", auth_token = "my_PAT")
```

We recommend using multiple cores if available when installing **salmonIPM** to reduce compilation time. You can do this by setting the R environment variable `MAKEFLAGS = -jX`, where `X` is the number of cores. This can be done interactively using `Sys.setenv(MAKEFLAGS = "-jX")` or it can be specified in `.Renviron`.
