
n.reps = 9

context("Gillespie")


test_that("random.system passes check.system")
{
    for (i in 1:n.reps) 
    {
        check.system(random.system(), info=paste("i=",i))
    }
}


sys = list()
sys$X = c(0,1,2)
sys$Pre = rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
sys$s = length(sys$X)
sys$r = nrow(sys$Pre)
sys$theta = numeric(sys$r)
sys$Post = matrix(0,nrow=sys$r,ncol=sys$s)
sys$stoich = t(sys$Post-sys$Pre)
res = c(1,0,1,2,0,0,2,0)

test_that("hazard.part passes test cases", {
    out = hazard.part(sys,engine="R")
    for (i in 1:length(res)) expect_equal(out[i],res[i],info=paste("R: i=",i))
    out = hazard.part(sys,engine="C")
    for (i in 1:length(res)) expect_equal(out[i],res[i],info=paste("C: i=",i))
})



