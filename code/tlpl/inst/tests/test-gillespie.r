
n.reps = 9

context("Gillespie")


test_that("random.system passes check.system",
{
    for (i in 1:n.reps) 
    {
        expect_true({check.system(random.system());TRUE}, info=paste("i=",i))
    }
})


sys = list()
sys$X = c(0,1,2)
sys$Pre = rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
sys$s = length(sys$X)
sys$r = nrow(sys$Pre)
sys$theta =  rgamma(sys$r, 100, 100)
sys$Post = matrix(0,nrow=sys$r,ncol=sys$s)
sys$stoich = t(sys$Post-sys$Pre)
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
        expect_equal(hazard(sys,engine="R"), hazard(sys,engine="C"))
    }
})

