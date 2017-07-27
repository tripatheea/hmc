import numpy as np
from pyhmc import hmc
import corner


# define your probability distribution
def logprob(x, ivar):
    logp = -0.5 * np.sum(ivar * x**2)
    grad = -ivar * x
    return logp, grad


# run the sampler
ivar = 1. / np.random.rand(5)
samples = hmc(logprob, x0=np.random.randn(5), args=(ivar,), n_samples=1e4)


# Optionally, plot the results
figure = corner.corner(samples)
figure.savefig('triangle.pdf')
