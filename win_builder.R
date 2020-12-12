library(devtools)
devtools::check_win_devel(pkg = ".",
                          args = NULL,
                          manual = TRUE,
                          email = NULL,
                          quiet = FALSE)

devtools::check_win_release(pkg = ".",
                          args = NULL,
                          manual = TRUE,
                          email = NULL,
                          quiet = FALSE)

devtools::run_examples(pkg = ".",
                       start = NULL,
                       run_dontrun = FALSE,
                       run_donttest = FALSE,
                       document = TRUE)

devtools::spell_check(pkg = ".")
devtools::check()
rhub::validate_email().
devtools::check_rhub(pkg = ".", email = "caledo@uma.es")

## --------- Eventually --------- ##
devtools::release(pkg = ".")
devetools::submit_cran(pkg = ".")
## ------------------------------ ##
# Run devtools::submit_cran() to re-submit the package without working through all the release() questions a second time.
