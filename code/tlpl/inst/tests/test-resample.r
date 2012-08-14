
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
    id = c(0,0,0,1,1,2)
    expect_equal(rep2id(rep    ),id)
    expect_equal(rep2id(rep,"R")-1,id)
    expect_equal(rep2id(rep,"C"),id)

    for (i in 1:10) 
    {
        rep = rpois(100,1)
        expect_equal(rep2id(rep,"R")-1,
                     rep2id(rep,"C"))
    }
})



test_that("inverse.cdf.weights works properly", {
    w = rep(1/4,4)
    u = c(.1,.3,.6,.8)
    id = 0:3
    expect_equal(inverse.cdf.weights(w,u    )  ,id) 
    expect_equal(inverse.cdf.weights(w,u,"R")-1,id)
    expect_equal(inverse.cdf.weights(w,u,"C")  ,id)

    u = c(.3,.1,.8,.6)
    expect_equal(inverse.cdf.weights(w,u    )  ,id)
    expect_equal(inverse.cdf.weights(w,u,"R")-1,id)
    expect_equal(inverse.cdf.weights(w,u,"C")  ,id)


    for (i in 1:10) 
    {
        w = runif(rpois(1,10)+1)
        w = w/sum(w)
        u = runif(rpois(1,10)+1)
        expect_equal(inverse.cdf.weights(w,u,"R")-1,
                     inverse.cdf.weights(w,u,"C"))
    }

    u = numeric(0) 
    expect_error(inverse.cdf.weights(w,u,"R"))
    expect_error(inverse.cdf.weights(w,u,"C"))
    u = runif(4); u[2] = -u[2]
    expect_error(inverse.cdf.weights(w,u,"R"))
    expect_error(inverse.cdf.weights(w,u,"C"))
    u = runif(4)
    w = numeric(0)
    expect_error(inverse.cdf.weights(w,u,"R"))
    expect_error(inverse.cdf.weights(w,u,"C"))
    w = rep(1/4,4); w[3] = -w[3]
    expect_error(inverse.cdf.weights(w,u,"R"))
    expect_error(inverse.cdf.weights(w,u,"C"))
})

