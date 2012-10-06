context("tlpl")

test_that("tlpl throws errors", {
    expect_error(tlpl())
    expect_error(tlpl(1:6))
    expect_error(tlpl(1:6,1:6))
})

test_that("tlpl completes", {
    sys = test.system(1)
    dat = list()
    dat$y = tau.leap(sys)
    dat$tau = 1
    expect_equal(tryCatch(tlpl(dat,sys)),0)
})

