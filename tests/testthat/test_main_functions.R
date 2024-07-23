context("ASTUTE")

set.seed(12345)
data(datasetExample)
resExample <- ASTUTE( alterations = datasetExample$alterations, 
                   expression = datasetExample$expression, 
                   regularization = TRUE, 
                   nboot = NA, 
                   num_processes = NA, 
                   verbose = FALSE )

data("resExample")
test_that("ASTUTE returns the correct output", {
    expect_equal(names(resExample),c("input_data","inference","parameters","goodness_fit","fold_changes","pvalues","qvalues"))
})
