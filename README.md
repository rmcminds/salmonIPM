**note that this fork contains experimental features by rmcminds**

# salmonIPM

This is the development repo for a fork of **salmonIPM**, an R package that fits integrated population models to data from anadromous Pacific salmonid populations using a hierarchical Bayesian framework implemented in [Stan](https://mc-stan.org/). Various models are available, representing alternative life histories and data structures as well as independent or hierarchically pooled populations. Users can specify stage-specific covariate effects and hyper-priors using formula syntax.

## Installation

1. Install and configure **rstan** (version 2.26 or higher) by following the instructions in the [RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) vignette.

2. Install the current version of **salmonIPM** from GitHub using **devtools**. 

```r
if(!require("devtools")) install.packages("devtools")
devtools::install_github("rmcminds/salmonIPM")
```

We recommend using multiple cores if available when installing **salmonIPM** to reduce compilation time. You can do this by setting the R environment variable `MAKEFLAGS` to `-jX`, where `X` is the number of cores. This can be done interactively using `Sys.setenv(MAKEFLAGS = "-jX")` or it can be specified in `.Renviron`.

## Citing **salmonIPM**

Buhle, E. R. and M. D. Scheuerell. 2024. salmonIPM: Integrated population models for Pacific salmonids (version x.x.x). https://doi.org/10.5281/zenodo.14511463

Buhle, E. R., M. D. Scheuerell, T. D. Cooney, M. J. Ford, R. W. Zabel, and J. T. Thorson. 2018. Using integrated population models to evaluate fishery and environmental impacts on Pacific salmon viability. U.S. Department of Commerce, NOAA Technical Memorandum, NMFS‐NWFSC‐14. http://doi.org/10.7289/V5/TM-NWFSC-140

Ohlberger, J., E. R. Buhle, T. W. Buehrens, N. W. Kendall, T. Harbison, A. M. Claiborne, J. P. Losee, J. Whitney, and M. D. Scheuerell. 2025. Declining marine survival of steelhead trout linked to climate and ecosystem change. Fish and Fisheries x:x-x. https://doi.org/10.1111/faf.12878

Buonanduci, M. S., E. R. Buhle, M. J. Case, E. R. Howe, J. C. Robertson, N. VanBuskirk, and A. K. Ettinger. 2025. Forest restoration can bolster salmon population persistence under climate change. Biological Conservation 305: 111099. https://doi.org/10.1016/j.biocon.2025.111099

[![DOI](https://zenodo.org/badge/84359284.svg)](https://doi.org/10.5281/zenodo.14511463)
