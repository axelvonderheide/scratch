from scipy.stats.distributions import norm, uniform # import normal distribution
import pylab as p  # import matplotlib with matlab interface
import numpy as np # import numpy package
from time import clock

if __name__ == '__main__':
    n_rv = 2 # number of random variables
    n_int = 5 # number of discretization points

    # set the mean and standard deviation of the two random variables
    m_la, std_la = 10.0, 1.0
    m_xi, std_xi = 1.0, 0.1

    # construct the normal distributions and get the methods
    # for the evaluation of the probability density functions
    g_la = norm(loc = m_la, scale = std_la)
    g_xi = norm(loc = m_xi, scale = std_xi)

    # discretize the range (-1,1) symmetrically with n_int points
    theta_arr = np.linspace(-(1.0 - 1.0 / n_int), 1.0 - 1.0 / n_int , n_int)
    # cover the random variable symmetrically around the mean 
    theta_la = m_la + 4 * std_la * theta_arr
    theta_xi = m_xi + 4 * std_xi * theta_arr
    # get hte size of the integration cell
    d_la = (8 * std_la) / n_int
    d_xi = (8 * std_xi) / n_int

    def Heaviside(x):
        ''' Heaviside function '''
        return (np.sign(x) + 1.0) / 2.0

    def q(e, la, xi):
        ''' Response function of a single fiber '''
        return la * e * Heaviside(xi - e)

    # prepare the sequence of the control strains in a numpy array
    e_arr = np.linspace(0, 2, 100)
    # define an array of the same size as e_arr
    mu_q_arr = np.zeros_like(e_arr)

    def loop_mu_q_e():
        # loop over the control variable (strain)
        for i, e in enumerate(e_arr):
            mu_q_e = 0.0  # interim variable for the summed stress contributions
            for la in theta_la: # loop over lambda range (array of values)
                for xi in theta_xi: # loop over xi range (array of values)
                    dG = g_la.pdf(la) * g_xi.pdf(xi) * d_la * d_xi
                    mu_q_e += q(e, la, xi) * dG
            mu_q_arr[ i ] = mu_q_e

    start_time = clock()
    #loop_mu_q_e()
    print 'loop-based: elapsed time', clock() - start_time

    dG_la = g_la.pdf(theta_la) * d_la
    dG_xi = g_xi.pdf(theta_xi) * d_xi

    dG_grid = dG_la[:, None] * dG_xi[None, :]

    def mu_q_e(e):
        ''' Summation / integration  over the random domain '''
        q_e_grid = q(e, theta_la[:, None], theta_xi[None, :])
        q_dG_grid = q_e_grid * dG_grid # element by element product of two (m,m) arrays
        return np.sum(q_dG_grid) # nsum has been imported at line 3 from numpy 

    p.subplot(121)
    mu_q_e_vct = np.vectorize(mu_q_e)
    start_time = clock()
    mu_q_arr = mu_q_e_vct(e_arr)
    print 'Regular grid of random variables: elapsed time', clock() - start_time
    p.plot(e_arr, mu_q_arr, color = 'blue', label = '$\\theta_j$ grid')

    def get_mu_q(dG, *theta):
        '''Generate an integrator method for the particular data type
        of dG and variables. 
        '''
        def mu_q(e):
            '''Template for the evaluation of the mean response.
            '''
            Q_dG = q(e, *theta)
            Q_dG *= dG # in-place multiplication
            return np.sum(Q_dG)
        return np.vectorize(mu_q)

    # Grid of constant probabilities
    pi_arr = np.linspace(0.5 / n_int, 1. - 0.5 / n_int, n_int)
    theta_la_ppf = g_la.ppf(pi_arr)
    theta_xi_ppf = g_xi.ppf(pi_arr)
    n_sim = n_int ** n_rv
    dG = 1.0 / n_sim
    mu_q_e_ppf = get_mu_q(dG, theta_la_ppf[:, np.newaxis], theta_xi_ppf[np.newaxis, :])
    start_time = clock()
    mu_q_arr = mu_q_e_ppf(e_arr)
    print 'Grid of constant probabilities: elapsed time', clock() - start_time
    p.plot(e_arr, mu_q_arr, color = 'green', label = '$\\pi_j$ grid')

    # Monte-Carlo implementation
    theta_la_rvs = g_la.rvs(n_sim)
    theta_xi_rvs = g_xi.rvs(n_sim)
    mu_q_e_rvs = get_mu_q(1.0 / n_sim, theta_la_rvs, theta_xi_rvs)
    start_time = clock()
    mu_q_arr = mu_q_e_rvs(e_arr)
    print 'Monte-Carlo: elapsed time', clock() - start_time
    p.plot(e_arr, mu_q_arr, color = 'red', label = 'Monte-Carlo')
    p.legend()
    p.xlabel('$\\varepsilon$', fontsize = 24)
    p.ylabel('$q$', fontsize = 24)

    ############################## Discretization grids ##################################
    p.subplot(122)
    expander = np.ones((n_int, n_int), dtype = int)
    p.plot((theta_la[np.newaxis, :] * expander).flatten(),
            (theta_xi[:, np.newaxis] * expander).flatten(),
            'k.', label = '$\\theta_j$ grid')
    p.plot((theta_la_ppf[np.newaxis, :] * expander).flatten(),
            (theta_xi_ppf[:, np.newaxis] * expander).flatten(),
            'o', color = 'gray', label = '$\\pi_j$ grid')
    p.plot(theta_la_rvs, theta_xi_rvs, 'kD', label = 'Monte-Carlo')
    p.ylabel('$\\theta_{\\xi}$', fontsize = 24)
    p.ylim(0.5, 1.5)
    p.xlim(5, 15)
    p.xlabel('$\\theta_{\lambda}$', fontsize = 24)
    p.legend()

    p.show()

