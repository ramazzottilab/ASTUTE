Sys.setenv("R_TESTS" = "")

library("testthat")
library("ASTUTE")

test_check("ASTUTE")
