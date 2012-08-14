n.reps = 9

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
    
    for (i in 1:n.reps) 
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

    for (i in 1:n.reps) 
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

    for (i in 1:n.reps) 
    {
        rep = rpois(100,1)
        expect_equal(rep2id(rep,"R")-1,
                     rep2id(rep,"C"))
    }
})

test_that("inverse.cdf.weights throws proper errors", {
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


    for (i in 1:n.reps) 
    {
        w = runif(rpois(1,10)+1)
        w = w/sum(w)
        u = runif(rpois(1,10)+1)
        expect_equal(inverse.cdf.weights(w,u,"R")-1,
                     inverse.cdf.weights(w,u,"C"))
    }
})




test_that("renormalize throws errors", {
    w = runif(4); w[2]=-w[2]
    expect_error(renormalize(w))
})

test_that("renormalize works properly", {
    w = rep(1/4,4)
    expect_equal(renormalize(w,   ),w)   
    expect_equal(renormalize(w,engine="R"),w)   
    expect_equal(renormalize(w,engine="C"),w)   

    w = runif(4)
    lw = log(w) 
    w = w/sum(w)
    expect_equal(renormalize(lw,T,   ),w)   
    expect_equal(renormalize(lw,T,"R"),w)   
    expect_equal(renormalize(lw,T,"C"),w)   

    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1)
        expect_equal(renormalize(w,engine="R"),
                     renormalize(w,engine="C"))

        lw = log(w)
        expect_equal(renormalize(w,T,"R"),
                     renormalize(w,T,"C"))
    }
})





###############################################################
# Effective sample size functions
###############################################################

context("Sample size functions")


test_that("ess throws proper errors", {
    w = numeric(0)
    expect_error(ess(w    ))
    expect_error(ess(w,"R"))
    expect_error(ess(w,"C"))

    w = runif(4); w[2] = -w[2]
    expect_error(ess(w    ))
    expect_error(ess(w,"R"))
    expect_error(ess(w,"C"))
})

test_that("ess works properly", {
    n = 4
    w = rep(1/n, n)
    expect_equal(ess(w,   ), n) 
    expect_equal(ess(w,"R"), n) 
    expect_equal(ess(w,"C"), n) 

    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1); w=w/sum(w)
        expect_equal(ess(w,engine="R"),
                     ess(w,engine="C"))
    }

})


test_that("entropy throws proper errors", {
    w = numeric(0)
    expect_error(entropy(w    ))
    expect_error(entropy(w,"R"))
    expect_error(entropy(w,"C"))

    w = runif(4); w[2] = -w[2]
    expect_error(entropy(w    ))
    expect_error(entropy(w,"R"))
    expect_error(entropy(w,"C"))
})

test_that("entropy works properly", {
    n = 4
    w = rep(1/n, n)
    ent = -log2(1/n)
    expect_equal(entropy(w    ), ent) 
    expect_equal(entropy(w,"R"), ent) 
    expect_equal(entropy(w,"C"), ent) 

    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1); w=w/sum(w)
        expect_equal(entropy(w,engine="R"),
                     entropy(w,engine="C"))
    }

})


test_that("cov2 throws proper errors", {
    w = numeric(0)
    expect_error(cov2(w    ))
    expect_error(cov2(w,"R"))
    expect_error(cov2(w,"C"))

    w = runif(4); w[2] = -w[2]
    expect_error(cov2(w    ))
    expect_error(cov2(w,"R"))
    expect_error(cov2(w,"C"))
})

test_that("cov2 works properly", {
    n = 4
    w = rep(1/n, n)
    cov2 = 0
    expect_equal(cov2(w    ), cov2) 
    expect_equal(cov2(w,"R"), cov2) 
    expect_equal(cov2(w,"C"), cov2) 

    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1); w=w/sum(w)
        expect_equal(cov2(w,engine="R"),
                     cov2(w,engine="C"))
    }

})


###############################################################
# Resampling functions
###############################################################


context("Resampling functions")

test_that("multinomial resampling works properly", {
    set.seed(1)
    w = runif(4); w=w/sum(w)
    n = 5
    id = c(1,1,3,3,3)
    set.seed(2); expect_equal(multinomial.resample(w,n,   ),   id)
    set.seed(2); expect_equal(multinomial.resample(w,n,"R")-1, id)
    set.seed(2); expect_equal(multinomial.resample(w,n,"C"),   id)

    w = runif(rpois(1,10)+2); w=w/sum(w)
    n = rpois(1,10)+1
    seed = proc.time()
    set.seed(seed)
    
})


