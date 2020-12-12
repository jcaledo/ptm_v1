## Resubmission

This is a  resubmission. Now, in the current version:

* I've elaborated a bit on the description field as requested.
* Code lines in examples (in parse.dssp.Rd) have been uncommented.
* The line in man/foldx.assembly.Rd has been fixed as requested.
* We are only using suppressWarning() when necessary.
* In ddG.R and ptmplot.R we make use of on.exit().
* In the functions that offer the option of saving the results in a file, the user has to explicitly pass an argument to allow it.

## Test environments

* local OS X install, R 4.0.2

* win-builder (devel and release)

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies

There are currently no downstream dependencies for this package.
