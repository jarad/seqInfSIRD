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
    dat$y = matrix(1:n, n, 1)
    dat$tau = 1
    res = tlpl(dat,sys)
    print(res)
    np = dim(res$X)[3]; 
    expect_equal(res$hyper$rate$b, array(c(0,dat$y)+1, dim=c(n+1,1,np)))
    for (i in 1:np) expect_equal(res$hyper$rate$b[,,i], c(0,dat$y)+1)
})




