harare = c(2, 2, 6, 8, 8, 8, 9, 10, 10, 24, 29, 52, 52, 71, 89, 100, 103, 105, 109, 110, 111, 112, 112, 113, 117, 119, 119, 126, 129, 131, 132, 134, 134, 135, 135, 146, 148, 149, 150, 151, 151, 152, 152, 152, 153, 153, 153, 153, 153, 154, 154, 154, 154, 154, 155, 156, 156, 156, 156, 156)

plot(harare)


############################### Important dates #####################

# http://www.who.int/hac/crises/zwe/sitreps/24may_2june2010/en/index.html
week.meeting     = 27 # meeting when mass vaccination campaign was decided upon
week.vaccination = 31 # start of mass vaccination campaign




############################### Priors ############################

# http://www.who.int/entity/hac/crises/zwe/sitreps/zimbabwe_epi_w42_17oct2009.pdf
coverage1 = 0.72
coverage2 = 0.92 # vaccination coverage during National Immunization days (June 2009)

# http://en.wikipedia.org/wiki/Harare
harare.pop = 1.6e6
metrop.pop = 2.8e6

# expected susceptibles
exp.sus = harare.pop*(1-coverage2) # metrop.pop*(1-coverage1) 




# expected S->I sampling probability
p.stoi = 517/8173 # confirmed to suspected ratio
# cheating since this uses data until the end of the outbreak

m.p.stoi   = 10
p.stoi.a = m.p.stoi*p.stoi
p.stoi.b = m.p.stoi

curve(dbeta(x,p.stoi.a,p.stoi.b))
legend("topright", legend = c(paste("  2.5%=", round(qbeta(.025,p.stoi.a,p.stoi.b),4)), 
                              paste("97.5%=",  round(qbeta(.975,p.stoi.a,p.stoi.b),4))))

# expected I->R sampling probability
# equal to zero
p.itor.a = 1e-6
p.itor.b = 1e6


# expected infected
observed = 2
exp.inf = observed/p.stoi


# expected S->I rate


# expected I->R rate
# www.cdc.gov/vaccines/pubs/pinkbook/downloads/meas.pdf
# 'Rash onset 2-4 days after prodrome'
# 'Measles virus is shed from the nasopharynx beginning with the prodrome until 3–4 days after rash onset. '

r.itor = 1 # virus shedding lasts for approximately one week 
m.r.itor = 50
r.itor.a = m.r.itor * r.itor
r.itor.b = m.r.itor

curve(dgamma(x,r.itor.a, r.itor.b), 0,2, main="Prior for I->R rate", xlab="I->R rate", ylab="density")
legend("topright", legend = c(paste("  2.5%=", round(qgamma(.025,r.itor.a,r.itor.b),2)), 
                              paste("97.5%=",  round(qgamma(.975,r.itor.a,r.itor.b),2))))



# expected S->I rate
# http://en.wikipedia.org/wiki/Basic_reproduction_number
# set so that basic reproductive number (R0) is 12-18
r.stoi   = 15
m.r.stoi = 20
r.stoi.a = m.r.stoi * r.stoi
r.stoi.b = m.r.stoi

curve(dgamma(x,r.stoi.a, r.stoi.b), 10,20, main="Prior for S->I rate", xlab="S->I rate", ylab="density")
legend("topright", legend = c(paste("  2.5%=", round(qgamma(.025,r.stoi.a,r.stoi.b),2)), 
                              paste("97.5%=",  round(qgamma(.975,r.stoi.a,r.stoi.b),2))))

# Implied R0
n = 1e6
r0 = rgamma(n, r.stoi.a, r.stoi.b)/rgamma(n, r.itor.a, r.itor.b)
hist(r0, 1001, freq=F, main="Implied prior on R0", xlab="R0", ylab="density")
legend("topright", legend=c(paste("  2.5%=", round(quantile(r0,.025),1)), 
                            paste("97.5%=" , round(quantile(r0,.975),1))))


