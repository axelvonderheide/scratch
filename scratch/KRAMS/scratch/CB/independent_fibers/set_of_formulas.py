import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import weibull_min
from scipy.special import gammainc, gamma
from scipy.integrate import cumtrapz
from math import pi

def H(x):
    return x >= 0.0

Ef = 200e3
r = 0.003
tau = .5
m = .5
s = 50.
L0 = 1.
sV0 = 4.e-8

def CDFw(e, L):
    return 1 - np.exp(-L/L0 * (e/s)**m)

# failure probability between 0 and a at e
def CDFa(e):
    T = 2. * tau / r
    a = e / T
    return 1 - np.exp(-a * 2 * Ef * (e / s) ** m / (m + 1) / L0)

# failure probability between 0 and a at e
def CDFa2(e):
    T = 2. * tau / r
    s0 = ((T * (m+1) * sV0**m)/(2. * Ef * 3.14159 * r ** 2))**(1./(m+1))
    return 1 - np.exp(-(e / s0) ** (m+1))

# failure probability density between 0 and a at e
def PDFa(e):
    T = 2. * tau / r
    a = e * Ef / T
    return 2 * (e / s) ** m / (T/Ef * L0) * np.exp(-a * 2 * (e / s) ** m / (m + 1) / L0)

# failure probability density between 0 and a at e
def PDFa2(e):
    T = 2. * tau / r
    s0 = ((T * (m+1) * sV0**m)/(2. * Ef * 3.14159 * r ** 2))**(1./(m+1))
    return (m+1)/s0 * (e/s0)**m * (1 - CDFa2(e))

# percent point function
def PPFa(p):
    T = 2. * tau / r
    return (-0.5 * np.log(1.-p) * T / Ef * L0 * (m+1) * s**m) ** (1./(m+1))

# failure probability between 0 and L at e
def CDF_L(e, L):
    T = 2. * tau / r / Ef
    a = e / T
    Pf = 1. - np.exp(-2 * a * (e / s) ** m * (1 - (1 - L / a) ** (m + 1)) / (m + 1) / L0)
    return Pf / CDFa(e)

def mean_xi():
    T = 2. * tau / r
    n = (m+1)
    c = 2. * Ef / T / L0 / n / s**m
    mu_xi = c**(-1./n)/n * gamma(1./n)
    return mu_xi 

def mean_xi2():
    T = 2. * tau / r
    s0 = ((T * (m+1) * sV0**m)/(2. * Ef * 3.14159 * r ** 2))**(1./(m+1))
    mu_xi = s0 * gamma(1+1./(m+1))
    return mu_xi 

def h(e, x):
    a = e / 2. / tau * r * Ef
    h_xe = (1 - CDFa(e)) * H(a - x) * 2 * m * (e * (1. - x/a)) ** (m-1) / L0 / s ** m 
    return h_xe

def g_z(e, x):
    T = 2. * tau / r
    a = e * Ef / T
    pdf = m/a * (1-x/a)**(m-1) * H(a - x)
    return pdf

def cdfL(e, x):
    a = e / 2. / tau * r * Ef
    return 1 - (1 - x/a)**m

def muLate(e):
    T = 2. * tau / r / Ef
    a = e / T
    return a/(m+1)

def evans(T):
    suma = (L0/r*(s* Ef)**m*tau*(m+1))**(1./(m+1.))
    alpha = (T/suma)**(m+1)
    return L0/2. * ((s* Ef)/suma)**m * gamma((m+2)/(m+1)) * gammainc((m+2)/(m+1), alpha) / CDFa(T/Ef)

def muH(e):
    T = 2. * tau / r
    n = (m+1)
    c = 2. * Ef / T / L0 / n / s**m
    I = c**(-1./n)/n * gamma(1./n) * gammainc(1./n, c*e**n) - e*(1-CDFa(e))
    return I * Ef / T / n / CDFa(e)

def muH2(e):
    T = 2. * tau / r
    n = m+1
    s0 = ((T * (m+1) * sV0**m)/(2. * Ef * pi * r ** 2))**(1./(m+1))
    I = s0 * gamma(1. + 1./n) * gammainc(1. + 1./n, (e/s0)**n)
    return I * Ef / T / n / CDFa(e)

def simulation(e_max):
    T = 2. * tau / r / Ef
    a = e_max / T
    nz = 500
    z_arr = np.linspace(a/2./nz, a - a/2./nz, nz)
    eps = e_max - T * z_arr
    dz = a/nz
    fu_arr = weibull_min(m, scale=s * (dz/L0)**(-1./m)).ppf(np.minimum(np.random.rand(nz), np.random.rand(nz)))
    max_diff = np.max(eps - fu_arr)
    if max_diff > 0.0:
        L = z_arr[np.argmax(eps - fu_arr)]
        eu = e_max - max_diff
        return L, eu
    else:
        return None, None

def MC(n, e_arr):
    Ls = []
    es = []
    for i in range(n):
        l, e = simulation(e_arr[-1])
        if l == None:
            pass
        else:
            Ls.append(l)
            es.append(e)
    es = np.array(es)
    Ls = np.array(Ls)
    muL = []
    fact = []
    for i, ei in enumerate(e_arr):
        mask = es < ei
        mask2 = es > e_arr[i-1]
        muL.append(np.mean(Ls[mask]))
        a = ei * r * Ef / 2. / tau
        pullouts = Ls[mask * mask2]
        mean = np.mean(pullouts)
        if len(pullouts) != 0:
            plt.hist(pullouts, bins=30, normed=True)
            plt.plot(np.linspace(0,a/1.5,50), g_z(ei, np.linspace(0,a/1.5,50)))
            plt.show()
        fact.append(np.mean(a) / mean)
    mask3 = np.isnan(fact) == False
    print 'mean L/a ratio =', np.mean(np.array(fact)[mask3])
    return muL

e_arr = np.linspace(0.0, 0.2, 100)
#z_arr = np.linspace(0.01, 10, 200)
#plt.plot(z_arr, g_z(z_arr, 0.03))
#plt.plot(e_arr, evans(e_arr * Ef))
#plt.plot(e_arr, PDFa(e_arr))
#plt.plot(e_arr, PDFa2(e_arr), ls='dashed', lw=2)
plt.plot(e_arr, muH(e_arr))
plt.plot(e_arr, muH2(e_arr), ls='dashed', lw=2)
#plt.plot(e_arr, MC(20000, e_arr))
plt.show()
