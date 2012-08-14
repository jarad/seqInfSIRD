
context("Resampling")

test_that("is.increasing works properly", {
    expect_true(is.increasing(1:3))
    expect_true(is.increasing(1:3,"C"))
    expect_false(is.increasing(c(2,1,3)))
    expect_false(is.increasing(c(2,1,3),"C"))
    
    for (i in 1:10) 
    {
        v <- rnorm(3)
        expect_identical(is.increasing(v), is.increasing(v,"C"))
    }
})




