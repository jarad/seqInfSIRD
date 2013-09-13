
# Overall settings
set.seed(20121218)

# Simulation settings
n.sims = 100
n = 60 # Number of time points

# Filter settings
n.particles = 500
resampling.function = "multinomial"

# Quantile settings
# Equal-tail credible intervals
ci = c(.8,.95)
lt = (1-ci)/2
ut = 1-lt

probs = sort(c(lt,.5,ut))

