from pypolychord.output import PolyChordOutput
import os
import _pypolychord
import numpy as np


def default_prior(cube):
    return cube.copy()


def default_dumper(live, dead, logweights, logZ, logZerr):
    pass


def run_polychord(loglikelihood, nDims, nDerived,
                  prior=default_prior, dumper=default_dumper, **kwargs):
    """
    Runs PolyChord.

    For full details of nested sampling and PolyChord, please refer to:

    * PolyChord paper: http://arxiv.org/abs/1506.00171
    * Nested Sampling paper: http://projecteuclid.org/euclid.ba/1340370944

    To run in mpi mode, just run your script with mpirun as usual.
    Make sure that pypolychord is compiled with MPI:
    $ make veryclean
    $ make pypolychord MPI=1

    Users are also required to cite the PolyChord papers:
    arXiv:1502.01856
    arXiv:1506.00171
    in their publications.


    Parameters
    ----------

    loglikelihood: function
        This function computes the log-likelihood of the model and derived
        parameters (phi) from the physical coordinates (theta).

        Parameters
        ----------
        theta: numpy.array
            physical coordinate. Length nDims.

        Returns
        -------
        (logL, phi): (float, array-like)
            return is a 2-tuple of the log-likelihood (logL) and the derived
            parameters (phi). phi length nDerived.

        Returns
        -------
        logL: float
            log-likelihood
        TODO: what about []?

    nDims: int
        Dimensionality of the model, i.e. the number of physical parameters.

    nDerived: int
        The number of derived parameters (can be 0).

    settings: settings.Settings
        Object containing polychord settings

    Optional Arguments
    ------------------

    prior: function
        This function converts from the unit hypercube to the physical
        parameters.
        (Default: identity function => prior(cube) == cube )

        Parameters
        ----------
        cube: numpy.array
            coordinates in the unit hypercube. Length nDims.

        Returns
        -------
        theta: array-like
            physical coordinates. Length nDims.

    dumper: function
        This function gives run-time access to the posterior and live points.

        Parameters
        ----------
        live: numpy.array
            The live points and their loglikelihood birth and death contours
            Shape (nlive, nDims+nDerived+2)
        dead: numpy.array
            The dead points and their loglikelihood birth and death contours
            Shape (ndead, nDims+nDerived+2)
        logweights: numpy.array
            The posterior weights of the dead points
            Shape (ndead)
        logZ: float
            The current log-evidence estimate
        logZerr: float
            The current log-evidence error estimate
        
    Keyword arguments
    -----------------
    nlive: int
        (Default: nDims*25)
        The number of live points.
        Increasing nlive increases the accuracy of posteriors and evidences,
        and proportionally increases runtime ~ O(nlive).

    num_repeats : int
        (Default: nDims*5)
        The number of slice slice-sampling steps to generate a new point.
        Increasing num_repeats increases the reliability of the algorithm.
        Typically
        * for reliable evidences need num_repeats ~ O(5*nDims).
        * for reliable posteriors need num_repeats ~ O(nDims)

    nprior : int
        (Default: nlive)
        The number of prior samples to draw before starting compression.

    nfail : int
        (Default: nlive)
        The number of failed spawns before stopping nested sampling.

    do_clustering : boolean
        (Default: True)
        Whether or not to use clustering at run time.

    feedback : {0,1,2,3}
        (Default: 1)
        How much command line feedback to give

    precision_criterion : float
        (Default: 0.001)
        Termination criterion. Nested sampling terminates when the evidence
        contained in the live points is precision_criterion fraction of the
        total evidence.

    logzero : float
        (Default: -1e30)
        The loglikelihood value at which PolyChord considers points
        'unphysical', and excludes them at the prior level.

    max_ndead : int
        (Default: -1)
        Alternative termination criterion. Stop after max_ndead iterations.
        Set negative to ignore (default).

    boost_posterior : float
        (Default: 0.0)
        Increase the number of posterior samples produced.  This can be set
        arbitrarily high, but you won't be able to boost by more than
        num_repeats
        Warning: in high dimensions PolyChord produces _a lot_ of posterior
        samples. You probably don't need to change this

    posteriors : boolean
        (Default: True)
        Produce (weighted) posterior samples. Stored in <root>.txt.

    equals : boolean
        (Default: True)
        Produce (equally weighted) posterior samples. Stored in
        <root>_equal_weights.txt

    cluster_posteriors : boolean
        (Default: True)
        Produce posterior files for each cluster?
        Does nothing if do_clustering=False.

    write_resume : boolean
        (Default: True)
        Create a resume file.

    read_resume : boolean
        (Default: True)
        Read from resume file.

    write_stats : boolean
        (Default: True)
        Write an evidence statistics file.

    write_live : boolean
        (Default: True)
        Write a live points file.

    write_dead : boolean
        (Default: True)
        Write a dead points file.

    write_prior : boolean
        (Default: True)
        Write a prior points file.

    maximise : boolean
        (Default: False)
        Perform maximisation at the end of the run to find the maximum
        likelihood point and value

    compression_factor : double
        (Default: exp(-1))
        How often to update the files and do clustering

    synchronous : boolean
        (Default: True)
        Parallelise with synchronous workers, rather than asynchronous ones.
        This can be set to False if the likelihood speed is known to be
        approximately constant across the parameter space. Synchronous
        parallelisation is less effective than asynchronous by a factor ~O(1)
        for large parallelisation.

    base_dir : string
        (Default: 'chains')
        Where to store output files.

    file_root : string
        (Default: 'test')
        Root name of the files produced.

    grade_frac : List[float]
        (Default: [1])
        The amount of time to spend in each speed. 
        If any of grade_frac are <= 1, then polychord will time each sub-speed,
        and then choose num_repeats for the number of slowest repeats, and
        spend the proportion of time indicated by grade_frac. Otherwise this
        indicates the number of repeats to spend in each speed.

    grade_dims : List[int]
        (Default: nDims)
        The number of parameters within each speed.

    nlives : dict {double:int}
        (Default: {})
        Variable number of live points option. This dictionary is a mapping
        between loglike contours and nlive.
        You should still set nlive to be a sensible number, as this indicates
        how often to update the clustering, and to define the default value.
    seed : positive int
        (Default: system time in milliseconds)
        Choose the seed to seed the random number generator.
        Note **Positive seeds only**
        a negative seed indicates that you should use the system time in
        milliseconds
    cube_samples: array-like
        (Default: None)
        samples from the unit hypercube to start nested sampling from. This is
        useful if the prior is hard to rejection sample directly from the unit
        hypercube (e.g. if there is a large amount of excluded parameter
        space), but can be done more efficiently manually. A concrete example
        could be include sorted uniform parameters by requiring 0 < x1 < x2 <
        x3 < 1, if one did not want to use SortedUniformPrior.
        Only available in Python interface.
        shape (:, nDims)

    Returns
    -------
    None. (in Python)

    All output is currently produced in the form of text files in <base_dir>
    directory. If you would like to contribute to pypolychord and improve this,
    please get in touch:

    Will Handley: wh260@cam.ac.uk

    In general the contents of <base_dir> is a set of getdist compatible files.

    <root> = <base_dir>/<file_root>

    <root>.txt                                              (posteriors = True)
        Weighted posteriors. Compatible with getdist. Each line contains:

          weight, -2*loglikelihood, <physical parameters>, <derived parameters>

        Note that here the weights are not normalised, but instead are chosen
        so that the maximum weight is 1.0.

    <root>_equal_weights.txt                                    (equals = True)
        Weighted posteriors. Compatible with getdist. Each line contains:
            1.0, -2*loglikelihood, <physical parameters>, <derived parameters>

    <root>_dead.txt                                         (write_dead = True)
        Dead points. Each line contains:
            loglikelihood, <physical parameters>, <derived parameters>

    <root>_phys_live.txt                                    (write_live = True)
        Live points. Each line contains:
            <physical parameters>, <derived parameters>, loglikelihood
        Mostly useful for run-time monitoring.

    <root>.resume
        Resume file. Human readable.

    <root>.stats
        Final output evidence statistics




    """

    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    except ImportError:
        rank = 0

    kwargs.setdefault("base_dir", "chains")
    kwargs.setdefault("cluster_dir", "test")
    kwargs.setdefault("read_resume", True)
    try:
        if rank == 0:
            os.makedirs(kwargs["base_dir"])
    except OSError:
        pass
        
    try:
        if rank == 0:
            os.makedirs(os.path.join(kwargs["base_dir"], kwargs["cluster_dir"]))
    except OSError:
        pass

    if "cube_samples" in kwargs:
        make_resume_file(loglikelihood, prior, **kwargs)
        read_resume = kwargs["read_resume"]
        kwargs["read_resume"] = True

    def wrap_loglikelihood(theta, phi):
        logL, phi[:] = loglikelihood(theta)
        return logL

    def wrap_prior(cube, theta):
        theta[:] = prior(cube)

    kwargs["grade_dims"] = [int(d) for d in kwargs.pop("grade_dims", [nDims])]
    kwargs["nlives"] = {float(logL):int(nlive) for logL, nlive in kwargs.pop("nlives", {}).items()}

    # Run polychord from module library
    _pypolychord.run(wrap_loglikelihood,
        wrap_prior,
        dumper,
        nDims,
        nDerived,
        kwargs.pop('nlive', nDims*25),
        kwargs.pop('num_repeats', nDims*5),
        kwargs.pop('nprior', -1),
        kwargs.pop('nfail', -1), 
        kwargs.pop('do_clustering', True),
        kwargs.pop('feedback', 1),
        kwargs.pop('precision_criterion', 0.001),
        kwargs.pop('logzero', -1e30),
        kwargs.pop('max_ndead', -1),
        kwargs.pop('boost_posterior', 0.0),
        kwargs.pop('posteriors', True),
        kwargs.pop('equals', True),
        kwargs.pop('cluster_posteriors', True),
        kwargs.pop('write_resume', True),
        kwargs.pop('write_paramnames', False),
        kwargs.pop('read_resume', True),
        kwargs.pop('write_stats', True),
        kwargs.pop('write_live', True),
        kwargs.pop('write_dead', True),
        kwargs.pop('write_prior', True),
        kwargs.pop('maximise', False),
        kwargs.pop('compression_factor', np.exp(-1)),
        kwargs.pop('synchronous', True),
        kwargs.pop('base_dir', 'chains'),
        kwargs.pop('file_root', 'test'),
        kwargs.pop('seed', -1),
        kwargs.pop("grade_frac", [1.0]*len(kwargs["grade_dims"])),
        kwargs["grade_dims"],
        kwargs["nlives"], 
        kwargs.pop("seed", -1),
    )

    if "cube_samples" in kwargs:
        kwargs["read_resume"] = read_resume
    try:
        import anesthetic
        return anesthetic.NestedSamples(os.path.join(kwargs["base_dir"], kwargs["file_root"]))
    except ImportError:
        return PolyChordOutput(os.path.join(kwargs["base_dir"], kwargs["file_root"]))



def make_resume_file(loglikelihood, prior, **kwargs):
    import fortranformat as ff
    resume_filename = os.path.join(kwargs["base_dir"],
                                   kwargs["file_root"])+".resume"

    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    except ImportError:
        rank = 0
        size = 1

    lives = []
    logL_birth = kwargs["logzero"]
    for i in np.array_split(np.arange(len(kwargs["cube_samples"])), size)[rank]:
        cube = kwargs["cube_samples"][i]
        theta = prior(cube)
        logL, derived = loglikelihood(theta)
        nDims = len(theta)
        nDerived = len(derived)
        lives.append(np.concatenate([cube,theta,derived,[logL_birth, logL]]))

    try:
        sendbuf = np.array(lives).flatten()
        sendcounts = np.array(comm.gather(len(sendbuf)))
        if rank == 0:
            recvbuf = np.empty(sum(sendcounts), dtype=int)
        else:
            recvbuf = None
        comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)

        lives = np.reshape(sendbuf, (len(kwargs["cube_samples"]), len(lives[0])))
    except NameError:
        lives = np.array(lives)

    if rank == 0:
        with open(resume_filename,"w") as f:
            def write(var):
                var = np.atleast_1d(var)
                if isinstance(var[0], np.integer):
                    fmt = '(%iI12)' % var.size
                elif isinstance(var[0], np.double):
                    fmt = '(%iE24.15E3)' % var.size
                else:
                    fmt = '(A)'
                writer = ff.FortranRecordWriter(fmt)
                f.write(writer.write(var) + '\n')

            write('=== Number of dimensions ===')
            write(nDims)
            write('=== Number of derived parameters ===')
            write(nDerived)
            write('=== Number of dead points/iterations ===')
            write(0)
            write('=== Number of clusters ===')
            write(1)
            write('=== Number of dead clusters ===')
            write(0)
            write('=== Number of global weighted posterior points ===')
            write(0)
            write('=== Number of global equally weighted posterior points ===')
            write(0)
            write('=== Number of grades ===')
            write(len(kwargs["grade_dims"]))
            write('=== positions of grades ===')
            write(kwargs["grade_dims"])
            write('=== Number of repeats ===')
            write(kwargs["num_repeats"])
            write('=== Number of likelihood calls ===')
            write(len(lives))
            write('=== Number of live points in each cluster ===')
            write(len(lives))
            write('=== Number of phantom points in each cluster ===')
            write(0)
            write('=== Number of weighted posterior points in each cluster ===')
            write(0)
            write('=== Number of equally weighted posterior points in each cluster ===')
            write(0)
            write('=== Minimum loglikelihood positions ===')
            write(np.argmin(lives[:,-1]))
            write('=== Number of weighted posterior points in each dead cluster ===')
            write('=== Number of equally weighted posterior points in each dead cluster ===')
            write('=== global evidence -- log(<Z>) ===')
            write(kwargs["logzero"])
            write('=== global evidence^2 -- log(<Z^2>) ===')
            write(kwargs["logzero"])
            write('=== posterior thin factor ===')
            write(kwargs["boost_posterior"])
            write('=== local loglikelihood bounds ===')
            write(lives[:,-1].min())
            write('=== local volume -- log(<X_p>) ===')
            write(0.0)
            write('=== last update volume ===')
            write(0.0)
            write('=== global evidence volume cross correlation -- log(<ZX_p>) ===')
            write(kwargs["logzero"])
            write('=== local evidence -- log(<Z_p>) ===')
            write(kwargs["logzero"])
            write('=== local evidence^2 -- log(<Z_p^2>) ===')
            write(kwargs["logzero"])
            write('=== local evidence volume cross correlation -- log(<Z_pX_p>) ===')
            write(kwargs["logzero"])
            write('=== local volume cross correlation -- log(<X_pX_q>) ===')
            write(0.0)
            write('=== maximum log weights -- log(w_p) ===')
            write(kwargs["logzero"])
            write('=== local dead evidence -- log(<Z_p>) ===')
            write('=== local dead evidence^2 -- log(<Z_p^2>) ===')
            write('=== maximum dead log weights -- log(w_p) ===')
            write('=== covariance matrices ===')
            write('---------------------------------------')
            for x in np.identity(nDims):
                write(x)
            write('=== cholesky decompositions ===')
            write('---------------------------------------')
            for x in np.identity(nDims):
                write(x)
            write('=== live points ===')
            write('---------------------------------------')
            for x in lives:
                write(x)
            write('=== dead points ===')
            write('=== logweights of dead points ===')
            write('=== phantom points ===')
            write('---------------------------------------')
            write('=== weighted posterior points ===')
            write('---------------------------------------')
            write('=== dead weighted posterior points ===')
            write('=== global weighted posterior points ===')
            write('=== equally weighted posterior points ===')
            write('---------------------------------------')
            write('=== dead equally weighted posterior points ===')
            write('=== global equally weighted posterior points ===')
