Submission of rivr 1.2-2. 

This patch release modifies internal behavior to handle an upcoming change
to R 4.0.0, which makes `data.frame()` and `read.table()` use
a `stringsAsFactors = FALSE` default.

## Test environments

* Local Windows 10 install, R 3.6.2
* Ubuntu 14.04 (on travis-ci), R-oldrel, R-release, R-devel

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTES