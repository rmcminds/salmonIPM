# salmonIPM

This is the development repo for **salmonIPM**, an R package for fitting integrated population models to salmon data.

## Installation

You can install the current version from GitHub using `devtools`. Because the repo is private for the time being, it is necessary to generate a personal access token (PAT) and pass it to `install_github()` as described [here](https://stackoverflow.com/questions/21171142/how-to-install-r-package-from-private-repo-using-devtools-install-github).

```{r }
if(!require("devtools")) {
  install.packages("devtools")
  library("devtools")
}
devtools::install_github("ebuhle/salmonIPM", auth_token = "my_PAT")
```

We strongly recommend compiling **salmonIPM** with the following `MAKEVARS`:

```
CXX14FLAGS += -Wno-unused -Wno-ignored-attributes -Wno-sign-compare -Wno-deprecated-declarations
CXX14FLAGS += -Wno-attributes -Wno-parentheses
CXX14FLAGS += -mfpmath=sse -mstackrealign  # https://github.com/kaskr/adcomp/issues/321#issuecomment-780945941
CXX14FLAGS += -mtune=native -O3 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2
PKG_CPPFLAGS += -DUSE_STANC3 
```

For more information on how to  `MAKEVARS`, see the "Configuring C++ Toolchain" section of the [RStan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) vignette for your operating system.
