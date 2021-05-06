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

devtools::check_rhub()

rhub::platforms()
rhub::check_on_windows(path = ".")

devtools::spell_check(pkg = ".")

rhub::validate_email()

## -----------------------------------------------------##
devtools::check()
devtools::check_rhub(pkg = ".", email = "caledo@uma.es")
## -----------------------------------------------------##

## --------- Eventually --------- ##
devtools::release(pkg = ".")
devtools::submit_cran(pkg = ".")
## ------------------------------ ##
# Run devtools::submit_cran() to re-submit the package without working through all the release() questions a second time.
