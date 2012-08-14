
context("Resampling utility functions")

test_that("is.increasing works properly", {
    inc = 1:3
    expect_true(is.increasing(inc    ))
    expect_true(is.increasing(inc,"R"))
    expect_true(is.increasing(inc,"C"))

    not = c(2,1,3)
    expect_false(is.increasing(not    ))
    expect_false(is.increasing(not,"R"))
    expect_false(is.increasing(not,"C"))
    
    for (i in 1:10) 
    {
        v <- rnorm(rpois(1,1)+3)
        expect_identical(is.increasing(v,"R"), 
                         is.increasing(v,"C"), 
                         info=paste(v))
    }
})


test_that("cusum works properly", {
    v = 1:3
    cs = c(1,3,6)
    expect_equal(cusum(v    ),cs)
    expect_equal(cusum(v,"R"),cs)
    expect_equal(cusum(v,"C"),cs)

    for (i in 1:10) 
    {
        v = rnorm(rpois(1,10)+1)
        expect_equal(cusum(v,"R"),
                     cusum(v,"C"),
                     info=paste(v))
    }
})


test_that("rep2id works properly", {
    rep = c(3,2,1)
    id = c(1,1,1,2,2,3)
    expect_equal(rep2id(rep,   ),id-1) # Default is C
    expect_equal(rep2id(rep,"R"),id)
    expect_equal(rep2id(rep,"C"),id-1)

    for (i in 1:10) 
    {
        rep = rpois(100,1)
        expect_equal(rep2id(rep,"R")-1,
                     rep2id(rep,"C"))
    }
})


