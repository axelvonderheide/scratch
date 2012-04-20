'''
Created on Apr 21, 2010

@author: rostislav
'''

from enthought.traits.api import HasTraits, Float, Int, Event, Array, Interface, \
    Tuple, Property, cached_property, Instance, Enum, DelegatesTo
from numpy import random, linspace, sum, sort, max, sign, nonzero, array, sqrt, mean
from scipy.integrate import trapz
from matplotlib import pyplot as plt
from scipy.stats import weibull_min, uniform
import time

def Heaviside( x ):
    return sign( sign( x ) + 1 )

class BundleStrength(HasTraits):

    # number of filaments in a bundle
    Nf = Float(10)
    # filament stiffness
    E = Float(1)
    # Weibull shape parameter for breaking strain
    xi_shape = Float(4.2)
    # Weibull scale parameter for breaking strain
    xi_scale = Float(0.06)
    # uniform distr location parameter for slack (lower bound)
    theta_loc = Float(0)
    # uniform scale parameter for slack (upper bound)
    theta_scale = Float(0.03)
    # uniform distr location parameter for E (lower bound)
    E_loc = Float(100)
    # uniform scale parameter for E (upper bound)
    E_scale = Float(200)
    # MC simulated array of breaking strains
    xi = Property(Array, depens_on = 'xi_shape, xi_scale')
    def _get_xi(self):
        # Nf random numbers from the interval (0,1)
        rn = random.rand(self.Nf)   
        xi_rand = weibull_min.ppf(rn, self.xi_shape, scale = self.xi_scale)
        return xi_rand
    
    # MC simulated array of slacks
    theta = Property(Array, depens_on = 'shape, scale')
    def _get_theta(self):
        # Nf random numbers from the interval (0,1)
        rn = random.rand(self.Nf)   
        theta_rand = uniform.ppf(rn, self.theta_loc, scale = self.theta_scale)
        return theta_rand
    
    # strains to plot
    e = linspace(0.,0.13,500)
    
    peaks = [0,0,0]
    
    def monte_carlo(self):
        xi = self.xi
        y = []
        # random slack array with Nf elements
        theta = self.theta 
        for eps in self.e:
            a = (eps - theta)*Heaviside(eps - theta)
            y.append(sum(a[a<xi]))
        y_montecarlo = array(y)*self.E/self.Nf
        self.peaks[0] = max(y_montecarlo)
        self.peaks[1] = max(sort(self.xi)*linspace(self.Nf,1,self.Nf)/self.Nf)
        return (self.e, y_montecarlo)

#    def analytical(self):
#        # analytical stress/strain realtionship for an infinite number of filaments
#        # as the actual strain multiplyied by the survival probability and integrated
#        # over random slack
#        y_analytical = []
#        thetas = linspace(self.theta_loc, self.theta_scale, 10000000)
#        step = thetas[1] - thetas[0]
#        PDF = uniform.pdf(thetas, loc = self.theta_loc, scale = self.theta_scale)
#        for eps in self.e:
#            integ_term = (eps - thetas)*Heaviside(eps - thetas)*\
#            (1-weibull_min.cdf((eps - thetas), self.xi_shape, scale = self.xi_scale))
#            time.clock()
#            y_analytical.append(trapz(integ_term*PDF, dx = step))
#            print time.clock()
#        self.peaks[2] = max(y_analytical)
#        return self.e, array(y_analytical)*self.E
    
    def analytical(self):
        # analytical stress/strain relationship for an infinite number of filaments
        # as the actual strain multiplyied by the survival probability and integrated
        # over random slack
        y_analytical = []
        thetas = linspace(self.theta_loc, self.theta_scale, 1000)
        CDF = uniform.cdf(thetas, loc = self.theta_loc, scale = self.theta_scale)
        for eps in self.e:
             integ_term = (eps - thetas)*Heaviside(eps - thetas)*\
             (1-weibull_min.cdf((eps - thetas), self.xi_shape, scale = self.xi_scale))
             time.clock()
             y_analytical.append(trapz(integ_term, CDF))
             print time.clock()
        self.peaks[2] = max(y_analytical)
        
        return self.e, array(y_analytical)*self.E
    
#    def asymptotic(self):
#        # analytical stress/strain realtionship for an infinite number of filaments
#        # as the actual strain multiplyied by the survival probability and integrated
#        # over random slack
#        y_asymptotic = []
#        thetas = linspace(self.theta_loc, self.theta_scale, 50)
#        E = linspace(self.E_loc, self.E_scale, 50)
#        for eps in self.e:
#             integ_term = E*(eps - thetas)*Heaviside(eps - thetas)*\
#             (1-weibull_min.cdf((eps - thetas), self.xi_shape, scale = self.xi_scale))
#             y_asymptotic.append(trapz(integ_term/self.theta_scale/self.E_scale, thetas,E))
#        self.peaks[2] = max(y_asymptotic)
#        return self.e, array(y_asymptotic)

    def get_peak(self):
        # peak force computed as:
        # 1: max of the simulated values
        print 'maximum MC simulated =', self.peaks[0]
        # 2: max of upwards sorted strengths divided by the number of survived filaments
        print 'maximum from sorted simulated values =', self.peaks[1]
        # 3: maximum analyticaly taking 1 - CDF as damage function
        # for bundle with infinite number of filaments
        print 'bundle analytical =', self.peaks[2]

if __name__ == '__main__':
    bs = BundleStrength()
    xdata, ydata = bs.analytical()
    plt.plot(bs.monte_carlo()[0], bs.monte_carlo()[1], 'r', xdata, ydata, linewidth = 2)
    #plt.plot(bs.analytical()[0], bs.analytical()[1], 'r.', bs.asymptotic()[0], bs.asymptotic()[1], linewidth = 2)
    plt.legend(('MonteCarlo ' + str(bs.Nf) + ' fil.','asymptotic'))
    plt.title('bundle stress strain relationship')
    bs.get_peak()
    plt.show()