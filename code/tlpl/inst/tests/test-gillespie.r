
n.reps = 9

context("gillespie")


test_that("random.system passes check.system",
{
    for (i in 1:n.reps) 
    {
        expect_true({check.system(random.system());TRUE}, info=paste("i=",i))
    }
})

sys = test.system(2)
hp = c(1,0,1,2,0,0,2,0)

test_that("hazard.part passes test cases", {
    expect_equal(hazard.part(sys,engine="R"),hp)
    expect_equal(hazard.part(sys,engine="C"),hp)
})


h = hp*sys$theta
test_that("hazard passes test cases", {
    expect_equal(hazard(sys,engine="R"),list(h=h,hp=hp))
    expect_equal(hazard(sys,engine="C"),list(h=h,hp=hp))
})

test_that("hazard R and C match", {
    for (i in 1:n.reps) {
        sys = random.system()
        expect_equal(hazard(sys,engine="R"), 
                     hazard(sys,engine="C"))
    }
})

########### update.species ##########
nr = 0:7
xn = c(11,10,9)
test_that("update.species passes test cases", {
    expect_equal(update.species(sys,nr,engine="R"),xn)
    expect_equal(update.species(sys,nr,engine="C"),xn)
})

test_that("update.species R and C match", {
    for (i in 1:n.reps) {
        sys = random.system()
        nr = rpois(sys$r,1)
        expect_equal(update.species(sys,nr,engine="R"), 
                     update.species(sys,nr,engine="C"))
    }
})


########### tau.leap.one.step ###########
nr = c(0,0,1,2,0,0,4,0)
X = c(7,3,3)

test_that("tau.leap.one.step passes test cases", {
    expect_equal({set.seed(1); tau.leap.one.step(sys,engine="R")}, list(X=X,nr=nr))
    expect_equal({set.seed(1); tau.leap.one.step(sys,engine="C")}, list(X=X,nr=nr))
})

test_that("tau.leap.one.step R and C match", {
    for (i in 1:n.reps) {
        seed = sample(1e6,1)
        expect_equal({set.seed(seed); tau.leap.one.step(sys, engine="R")},
                     {set.seed(seed); tau.leap.one.step(sys, engine="C")})
    }
})

########### tau.leap ###########

#test_that("tau.leap passes test cases", {
#    expect_equal({set.seed(1); tau.leap.one.step(sys,engine="R")}, list(X=X,nr=nr))
#    expect_equal({set.seed(1); tau.leap.one.step(sys,engine="C")}, list(X=X,nr=nr))
#})

context("tau.leap")    
set.seed(proc.time()*1e6)
test_that("tau.leap R and C match", { 
    for (i in 1:n.reps) {
        seed = sample(1e6,1)
        #sys = random.system()
        n = rpois(1,1)+1
        tau = rgamma(1,100,100)
        expect_equal({set.seed(seed); tau.leap(sys, n, tau, engine="R")},
                     {set.seed(seed); tau.leap(sys, n, tau, engine="C")}, 
                     info=paste("seed= ",seed, ", n=", n, ", tau= ", tau))
    }
})


