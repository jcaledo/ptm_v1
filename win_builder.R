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


