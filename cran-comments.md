Submission of rivr 1.2-3. 

This patch release modifies the vignettes
to produce images in PNG rather than SVG
format. SVG support is optional in R and
thus not suitable for use in when generating package vignettes.

## Test environments

R-CMD-check using GitHub Actions and the following environments:

- {os: macOS-latest,   r: 'devel'}
- {os: macOS-latest,   r: 'release'}
- {os: macOS-latest,   r: 'oldrel'}
- {os: windows-latest, r: 'devel'}
- {os: windows-latest, r: 'release'}
- {os: windows-latest, r: 'oldrel'}
- {os: ubuntu-16.04,   r: 'devel', rspm: "https://packagemanager.rstudio.com/cran/__linux__/xenial/latest", http-user-agent: "R/4.0.0 (ubuntu-16.04) R (4.0.0 x86_64-pc-linux-gnu x86_64 linux-gnu) on GitHub Actions" }
- {os: ubuntu-16.04,   r: 'release', rspm: "https://packagemanager.rstudio.com/cran/__linux__/xenial/latest"}
- {os: ubuntu-16.04,   r: 'oldrel',  rspm: "https://packagemanager.rstudio.com/cran/__linux__/xenial/latest"}


## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTES