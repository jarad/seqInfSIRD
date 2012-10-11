context("tlpl")

test_that("tlpl throws errors", {
    expect_error(tlpl())
    expect_error(tlpl(1:6))
    expect_error(tlpl(1:6,1:6))
})


test_that("tlpl completes", {
    sys = test.system(1)
    dat = list()
    dat$y = matrix(1:6, 6, 1)
    dat$tau = 1
    expect_true({tlpl(dat,sys); TRUE})
})


test_that("tlpl: constitutive production", {
    sys = test.system(1)
    dat = list()
    n = 6
    dat$y = matrix(1, n, 1)
    dat$tau = 1
    res = tlpl(dat,sys,nonuniformity="ess",threshold=0)
    np = dim(res$X)[3];
    expect_equal(res$hyper$prob$a, array(1:(n+1), dim=c(n+1,1,np)))

    z =  apply(res$X, 3, diff)-matrix(dat$y, n, np)                   # unobserved transitions
    b = array(apply(rbind(rep(1,np), z), 2, cumsum), dim=c(n+1,1,np)) # updated hyper parameter
    expect_equal(res$hyper$prob$b, b)

    expect_equal(res$hyper$rate$a, res$X) 
    expect_equal(res$hyper$rate$b, array(1:(n+1), dim=c(n+1,1,np)))
})




