from numpy.random import randn, rand
from numpy import ones, sum as nsum, linspace, sqrt, exp, log
from matplotlib import pyplot as plt
from scipy.stats import norm, uniform
from math import pi
from scipy.special import erf

mu_K = 50.
a_K = 40.
b_K = 60.
sig_K = 10.



def transf(X):
    Y = 1./X
    return Y

def pdf_K(x):
    return norm(mu_K, sig_K).pdf(x)

def pdf_K_inv(x):
    return pdf_K(transf(x))/abs(1./transf(x)**2)

def pdf_K_invsum(x,N):
    return pdf_K_inv(x/N)/N

def pdf_K_series(x,N):
    return pdf_K_invsum(transf(x),N)/abs(1./transf(x)**2)

N = 5
K = linspace(20,80,200)
K_inv = transf(K)
K_inv_sum = K_inv*N
K_series = transf(K_inv_sum)
    
arr = 1./nsum(1./norm(mu_K, sig_K).ppf(rand(500000)))*100000

#density of the stiffness
plt.plot(K, pdf_K(K), lw = 2, ls = '--', color = 'red')

# density of the inverse stiffness
#plt.plot(K_inv, pdf_K_inv(K_inv))

# density of the sum of inverse stiffnesses
#plt.plot(K_inv_sum, pdf_K_invsum(K_inv_sum,N), lw = 2, ls = '--')

# density of the resulting stiffness of the serially coupled system
plt.plot(K_series, pdf_K_series(K_series,N), lw = 2, color = 'black')

# MC inverse stiffness
#plt.hist(arr, 200, normed = True, range = (0.0,100))

#print 0.25*(log(1) - log(0.2))
plt.show()