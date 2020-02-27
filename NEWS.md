# rivr 1.2-2

This patch release modifies internal behavior to handle an upcoming change
to R 4.0.0, which makes `data.frame()` and `read.table()` use
a `stringsAsFactors = FALSE` default. See
https://developer.r-project.org/Blog/public/2020/02/16/stringsasfactors/index.html
for more information about this change.


# rivr 1.2-1

This patch release recompiles the C++ code in response to 
updates to package Rcpp in order to generate new binaries on CRAN.


# rivr 1.2

This minor version update adds new documentation.

* New vignette which allows reproduction of R Journal article. 
* New shiny app for interactive exploration of gradually-varied flow 
  profiles.


# rivr 1.0

This major version update defines an S3 class for `rivr` and methods 
for printing, summarizing and plotting results of gradually-varied and 
unsteady flow simulations. A technical vignette has also been added 
that provides nitty-gritty details on the derivations and numerical
schemes used in the package.


# rivr 0.9.2

This minor version update improves error messages and clarifies
some of the documentation.


# rivr 0.9.1

This package is on CRAN now! The binaries should be available 
soon, but if you get an error try 
`install.packages("rivr", type = "source")`.
