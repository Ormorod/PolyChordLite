import os
from numpy import pi, log
import numpy as np
import pypolychord
from pypolychord.priors import UniformPrior
try:
    from mpi4py import MPI
except ImportError:
    pass


#| Define a four-dimensional spherical gaussian likelihood,
#| width sigma=0.1, centered on the 0 with one derived parameter.
#| The derived parameter is the squared radius

nDims = 4
nDerived = 1
sigma = 0.1

def likelihood(theta):
    """ Simple Gaussian Likelihood"""

    nDims = len(theta)
    r2 = sum(theta**2)
    logL = -log(2*pi*sigma*sigma)*nDims/2.0
    logL += -r2/2/sigma/sigma

    return logL, [r2]

#| Define a box uniform prior from -1 to 1

def prior(hypercube):
    """ Uniform prior from [-1,1]^D. """
    return UniformPrior(-1, 1)(hypercube)

#| Optional dumper function giving run-time read access to
#| the live points, dead points, weights and evidences

def dumper(live, dead, logweights, logZ, logZerr):
    print("Last dead point:", dead[-1])

#| Optional cluster function allow user-defined clustering

def cluster(points):
    """
    <do some clustering algorithm to assign clusters>
    - clusters should be an array of cluster labels for each point
    - each cluster should have at least one point
    - thus max(clusters)+1 should be the number of clusters
    - i.e. clusters are 0-indexed
    - work with the above numpy integer array
    """
    npoints = points.shape[0]
    clusters = np.full(npoints, -1, dtype=int)

    return clusters

#| Make an anesthetic plot 

paramnames = [('p%i' % i, r'\theta_%i' % i) for i in range(nDims)]
paramnames += [('r*', 'r')]

#| Run PolyChord

output = pypolychord.run(
    likelihood,
    nDims,
    nDerived=nDerived,
    prior=prior,
    dumper=dumper,
    file_root='gaussian',
    nlive=200,
    do_clustering=True,
    read_resume=False,
    paramnames=paramnames,
)

#| Make an anesthetic plot 

fig, axes = output.plot_2d(['p0','p1','p2','p3','r'])
fig.savefig('posterior.pdf')
