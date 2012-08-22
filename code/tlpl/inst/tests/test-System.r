context("System.r")

test_that("test.system throws error",
{
    expect_error(test.system(2),"Test system doesn't exist.")
    expect_error(test.system(3),"Test system doesn't exist.")
})

