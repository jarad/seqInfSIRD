
context("Resampling")

test_that("is.increasing works properly", {
    expect_true(is.increasing(1:3))
    expect_true(is.increasing(1:3,"C"))
    expect_false(is.increasing(c(2,1,3)))
    expect_false(is.increasing(c(2,1,3),"C"))
    
    for (i in 1:10) 
    {
        v <- rnorm(rpois(1,1)+3)
        expect_identical(is.increasing(v), is.increasing(v,"C"), 
                         info=paste(v))
    }
})


test_that("cusum works properly", {
    expect_equal(cusum(1:3),c(1,3,6))
    expect_equal(cusum(1:3,"C"),c(1,3,6))

    for (i in 1:10) 
    {
        v = rnorm(rpois(1,10)+1)
        expect_equal(cusum(v),cusum(v,"C"))
    }
})




