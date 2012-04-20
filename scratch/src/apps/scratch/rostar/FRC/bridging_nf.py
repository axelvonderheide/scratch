'''
Created on Aug 5, 2010

@author: rostislav
'''

from numpy import exp, linspace, trapz, frompyfunc, arange
from scipy.special import gammaln
from matplotlib import pyplot as plt


def binomln(n, k):
    ''' binomical coefficients with gammaln function
        allows for larger numbers to be computed ''' 
    return exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


def pmf_func(n, k, b):
    '''the exact prob. mass func. evaluated by the formula for
        beta-binomial distribution.'''
    # in the first tenth the distribution is approximated by a
    # uniform distribution bcs of numerical errors in this interval
    if k < n * b / 10:
        return pmf_func(n, n * b / 10, b)
    # the correct solution for larger k
    else:
        p = linspace(b / 1000., b, 1001)
        return binomln(n, k) / b * trapz(p ** k * (1 - p) ** (n - k), x = p)

def pmf_binom(n, k, p):
    '''for comparsion the binomial distribution with fixed probability''' 
    return binomln(n, k) * p ** k * (1 - p) ** (n - k)

def pmf_approx(n, k, b):
    '''Approximation with a uniform distribution - holds for large n,k'''
    if k < n * b:
        return 1. / ((n + 1) * b)
    else:
        return 0

# b is the upper limit of a uniform distribution of fiber lengt projection
# divided by the composite length Lc (e.g. 1m)
b = 0.03
# n is the number of fibers in the composite of length Lc
n = 1000
# vectorize the functions and plot
pmf = frompyfunc(pmf_func, 3, 1)
k = arange(0, 1.5 * n * b)
plt.plot(k, pmf(n, k, b), linewidth = 1, label = 'exact pmf')
frompy_approx = frompyfunc(pmf_approx, 3, 1)
plt.plot(k, frompy_approx(n, k, b), linewidth = 1, label = 'approx uniform pmf')
plt.plot(k, pmf_binom(n, k, b / 2.), label = 'binomial pmf with fixed p')
plt.legend(loc = 'best')
plt.show()
