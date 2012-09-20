context("pl-discrete")


sys = test.system()
part = test.particle()
y = c(1,0,1,1,0,0,1,0)
tau = 1


test_that("calc.pred.like passes test case",
{
    answer = -5.7807435157923
    expect_equal(calc.pred.like(y,tau,sys,part), answer)
    expect_equal(calc.pred.like(y,tau,sys,part,engine="C"), answer, info=paste(calc.pred.like(y,tau,sys,part,engine="C")))
})

