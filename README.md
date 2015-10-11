rivr: An R package for Open Channel Hydraulics
==============================================

[![Build Status](https://travis-ci.org/mkoohafkan/rivr.svg)](https://travis-ci.org/mkoohafkan/rivr)

This package is designed as an educational tool for students and instructors 
of undergraduate and graduate courses in open channel hydraulics. Functions are 
provided for computing flow and channel geometry, normal and critical depth, 
gradually-varied water-surface profiles (e.g. backwater curves) and unsteady 
flow (e.g. flood wave routing).

UPDATE 10/11/2015
---------------

Version 1.0 now on CRAN! This major version update defines an S3 class for 
`rivr` and methods for printing, summarizing and plotting results of 
gradually-varied and unsteady flow simulations. A technical vignette has also
been added that provides nitty-gritty details on the derivations and numerical
schemes used in the package.

UPDATE 4/2/2015
---------------

Version 0.9.2 now on CRAN! This version improves error messages and clarifies 
some of the documentation.



UPDATE 2/22/2015
----------------

This package is on CRAN now! The binaries should be available soon, but if you
get an error try `install.packages("rivr", type = "source")`.
