## Resubmission

This is a resubmission. Examples that make use of Internet resources may, occasionally, exhibit elapsed times greater than 10 seconds depending on the server state, which can lead to the appearance of NOTEs. To avoid that, now I have made use of ‘\donrunt{<example>}’ in the documentation of those function that make extensive use of Internet resources. The code of those functions that use Internet resources have also been modified to fail gracefully with an informative message if the resource is not available, instead of giving a check warning or error.



## Test environments

* local OS X install, R 4.0.2

* local Linux install, R 4.0.4

* win-builder (devel and release)

## R CMD check results

There were no ERRORs, WARNINGs or NOTEs. 


## Downstream dependencies

There are currently no downstream dependencies for this package.
