#File: tests/testhtat/test-inst-app.R
library(shinytest2)
test_that("tinyapp2_test",{
  app <- AppDriver$new(name= "tinyapp2")
  app$set_inputs(gene = "DSP")
  app$expect_values()
} )
