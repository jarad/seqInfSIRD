SIRpath.dat contains a sample path of the discrete-time
stochastic SIRD model. See Ludkovski and Niemi (Proceedings of
2011 Winter Simulation Conference) for full details. The
parameters are:

S->I : 0.75
I->R : 0.5
S->R : 0
I->D : 0.02

This scenario lasts for 76 periods. At each period we record:

X_t = (S,    I,  R,  D) states (first 4 columns)
Y_t = (-dS, dI, dR, dD) sampled states (last 4 columns)

The sampling proportions are fixed as c(0.25, 0.25, 0.5, 0.5). 
