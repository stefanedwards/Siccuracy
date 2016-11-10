library(testthat)
library(Siccuracy)


if (TRUE) {
  test_check("Siccuracy")
} else {
  test_path <- file.path(getwd(), 'tests', 'testthat')
  testthat:::run_tests('Siccuracy', test_path, filter='plink', reporter='check')
}