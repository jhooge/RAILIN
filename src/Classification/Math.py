'''
Created on Sep 9, 2010

@author: jhooge
'''
import numpy as np
import numpy.linalg as lin

EXP_MIN = -1e308
EXP_MAX = +709
LOG_MIN = 1e-300
LOG_MAX = 1e+300

def log(x, x_min=LOG_MIN, x_max=LOG_MAX):
    """
    Save version of  log, clips argument such that overflow does not occur.

    @param x: input
    @type x: numpy array or float or int

    @param x_min: lower value for clipping
    @type x_min: float

    @param x_max: upper value for clipping
    @type x_max: float
    """
    from numpy import log, clip

    x_min = max(x_min, LOG_MIN)
    x_max = min(x_max, LOG_MAX)

    return log(clip(x, x_min, x_max))

def exp(x, x_min=EXP_MIN, x_max=EXP_MAX):
    """
    Save version of exp, clips argument such that overflow does not occur.

    @param x: input
    @type x: numpy array or float or int

    @param x_min: lower value for clipping
    @type x_min: float

    @param x_max: upper value for clipping
    @type x_max: float
    """
    from numpy import exp, clip

    x_min = max(x_min, EXP_MIN)
    x_max = min(x_max, EXP_MAX)

    return exp(clip(x, x_min, x_max))

def log_sum_exp(x, axis=0):
    """
    Returns the logarithm of the sum of exponentials.

    @type x: Numpy array
    """
    xmax = x.max(axis) + EXP_MAX

    return log(exp(x - xmax).sum(axis)) + xmax

def outer_2way(x, y, func):
    """
    Outer product of two two-way tensors (matrices).
    That is, if shape(x) = (m,n), shape(y) = (q,n)
    the resulting tensor z will have shape (m,n,q)
    with elements: z[i,j,k] = func(x[i,j],y[k,j]).
    """
    from numpy import zeros
    
    if x.shape[1] != y.shape[1]:
        raise ValueError('Shapes do not match!')

    z = zeros((x.shape[0], x.shape[1], y.shape[0]))
    for i in range(z.shape[0]):
        z[i, :, :] = func(y, x[i]).T

    return z

# Error classes
class DenError(Exception):
    """Base class for exceptions in this module."""
    
    def __init__(self, message):
        """
        @param message: explanation of the error
        @type message: string
        """
        self.message = message
        Exception.__init__(self)
    
    def __str__(self):
        return self.message

# The following function do all the fancy stuff to check that parameters
# are Ok, and call the right implementation if args are OK.
def gauss_den(x, mu, va, log=False):
    """Compute multivariate Gaussian density at points x for 
    mean mu and variance va.

    @note: Vector are row vectors, except va which can be a matrix
    (row vector variance for diagonal variance).

    @param x: points where to estimate the pdf.  each row of the array is one
    point of d dimension
    @type x: numpy.ndarray
    @param mu: mean of the pdf. Should have same dimension d than points in x.
    @type mu: numpy.ndarray
    @param va: variance of the pdf. If va has d elements, va is interpreted as the
    diagonal elements of the actual covariance matrix. Otherwise,
    should be a dxd matrix (and positive definite).
    @type va: numpy.ndarray
    @param log: if True, returns the log-pdf instead of the pdf.
    @type log: bool
    
    @return: Returns a rank 1 array of the pdf at points x.
    @rtype: numpy.ndarray
    """
    
    lmu = np.atleast_2d(mu)
    lva = np.atleast_2d(va)
    lx = np.atleast_2d(x)
    
    #=======================#
    # Checking parameters   #
    #=======================#
    if len(np.shape(lmu)) != 2:
        raise DenError("mu is not rank 2")
        
    if len(np.shape(lva)) != 2:
        raise DenError("va is not rank 2")
        
    if len(np.shape(lx)) != 2:
        raise DenError("x is not rank 2")
        
    d = np.shape(lx)[1]
    (dm0, dm1) = np.shape(lmu)
    (dv0, dv1) = np.shape(lva)
    
#    print '= '+str(np.shape(lx)[1])
#    print '(dm0, dm1) = '+str(np.shape(lva))
#    print '(dv0, dv1) = '+str(np.shape(lva))
    
    # Check x and mu same dimension
    if dm0 != 1:
        msg = "mean must be a row vector!"
        raise DenError(msg)
    if dm1 != d:
        msg = "x and mu not same dim"
        raise DenError(msg)
    # Check va and mu same size
    if dv1 != d:
        msg = "mu and va not same dim"
        raise DenError(msg)
    if dv0 != 1 and dv0 != d:
        msg = "va not square"
        raise DenError(msg)

    #=============#
    # Computation #
    #=============#
    if d == 1:
        # scalar case
        return _scalar_gauss_den(lx[:, 0], lmu[0, 0], lva[0, 0], log)
    elif dv0 == 1:
        # Diagonal matrix case
        return _diag_gauss_den(lx, lmu, lva, log)
    elif dv1 == dv0:
        # full case
        return  _full_gauss_den(lx, lmu, lva, log)
    else:
        raise DenError("variance mode not recognized, this is a bug")

# Those 3 functions do almost all the actual computation
def _scalar_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    Call gauss_den instead"""
    d = mu.size
    inva = 1 / va
    fac = (2 * np.pi) ** (-d / 2.0) * np.sqrt(inva)
    inva *= -0.5
    y = ((x - mu) ** 2) * inva
    if not log:
        y = fac * np.exp(y)
    else:
        y += np.log(fac)

    return y
    
def _diag_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    Call gauss_den instead """
    # Diagonal matrix case
    d = mu.size
    #n   = x.shape[0]
    if not log:
        inva = 1 / va[0]
        fac = (2 * np.pi) ** (-d / 2.0) * np.prod(np.sqrt(inva))
        inva *= -0.5
        x = x - mu
        x **= 2
        y = fac * np.exp(np.dot(x, inva))
    else:
        # XXX optimize log case as non log case above
        y = _scalar_gauss_den(x[:, 0], mu[0, 0], va[0, 0], log)
        for i in range(1, d):
            y += _scalar_gauss_den(x[:, i], mu[0, i], va[0, i], log)
    return y

def _full_gauss_den(x, mu, va, log):
    """ This function is the actual implementation
    of gaussian pdf in full matrix case. 
    
    It assumes all args are conformant, so it should 
    not be used directly Call gauss_den instead
    
    Does not check if va is definite positive (on inversible 
    for that matter), so the inverse computation and/or determinant
    would throw an exception."""
    d = mu.size
    inva = lin.inv(va)
    fac = 1 / np.sqrt((2 * np.pi) ** d * np.fabs(lin.det(va)))

    # we are using a trick with sum to "emulate" 
    # the matrix multiplication inva * x without any explicit loop
    #y   = -0.5 * np.sum(np.dot((x-mu), inva) * (x-mu), 1)
    y = -0.5 * np.dot(np.dot((x - mu), inva) * (x - mu),
                       np.ones((mu.size, 1), x.dtype))[:, 0]

    if not log:
        y = fac * np.exp(y)
    else:
        y = y + np.log(fac)

    return y
    
def gaussian_pdf(x, mu, sigma, log):
    """ This function is the actual implementation
    of gaussian pdf in scalar case. It assumes all args
    are conformant, so it should not be used directly
    
    Call gauss_den instead"""

    y = - 0.5 * (((x - mu)/sigma)**2 + np.log(2*np.pi*sigma**2))
    if not log:
        return np.exp(y)
    else:
        return y

def laplace_pdf(x, mu, sigma, log):
    """ This function is the implementation
    of Laplace pdf in scalar case.
    """

    b = sigma / 2**0.5

    y = - np.fabs(x - mu)/b - np.log(2/b)
    if not log:
        return np.exp(y)
    else:
        return y

def normalize_rows(A):

    R = A.sum(1)
    A[:,:] = (A.T / (R + 1e-300)).T
    
def maximize(x):

    if not len(x):
        return 0.
    else:
        return max(x)
    
def sample_from_histogram(p, n_samples=1):
    """
    returns the indice of bin according to the histogram p

    @param p: histogram
    @type p: numpy.array
    @param n_samples: number of samples to generate
    @type n_samples: integer
    """
   
    from numpy import add, less, argsort, take, arange
    from numpy.random import random

    indices = argsort(p)
    indices = take(indices, arange(len(p) - 1, -1, -1))

    c = add.accumulate(take(p, indices)) / add.reduce(p)

    return indices[add.reduce(less.outer(c, random(n_samples)), 0)]

if __name__ == "main":
    pass