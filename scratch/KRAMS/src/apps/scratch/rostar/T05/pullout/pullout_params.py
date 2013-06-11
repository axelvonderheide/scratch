'''
Created on Jun 24, 2010

@author: rostislav
'''

from enthought.traits.api import HasTraits, Array, Float
from matplotlib import pyplot as plt
from numpy import sqrt, linspace, random, tanh, array, hstack, loadtxt, \
    transpose, arange, frompyfunc, polyfit, poly1d, polyder, argmax, mean, var
from scipy.interpolate import interp1d
from scipy.optimize import leastsq


class PullOutParams(HasTraits):

    # clamp length [m]
    L = Float(0.05)
    # elasticity modulus of the bundle [Pa]
    E_mod = Float(140.e9)
    # cross-sectional area of the bundle [m2]
    A = Float(0.894e-6)
    # test data: u-P
    test_u = Array
    test_P = Array

    ###################################
    ######## constant friction ########
    ###################################

    # pull-out model with a constant friction bond law
    # params: friction
    def const(self, qf):
        return sqrt(2. * self.test_u * self.E_mod * self.A * qf)

    # error function - difference between test data and fitting function
    def const_error(self, qf):
        return self.test_P - self.const(qf)

    # returns the parameters of the fitted curve
    def const_params(self):
        iparams = 200.
        return leastsq(self.const_error, iparams)

    ###################################
    ######## bilinear bond law ########
    ###################################

    # pull-out model with bilinear bond law
    # params: friction, bond shear stiffness
    def bilinear (self, k, qf):
        w = sqrt(k / self.E_mod / self.A)
        mask = self.test_u > qf / k
        test_u_deb = self.test_u[mask]
        P_interp = interp1d(array([0., qf / k]), array([0., qf / w]))
        mask2 = self.test_u <= qf / k
        test_u_lin = self.test_u[mask2]
        P_lin = P_interp(test_u_lin)
        P_deb = sqrt((2. * test_u_deb * w ** 2 * self.E_mod * self.A * qf - qf ** 2)) / w
        return hstack((P_lin, P_deb))

    # error function - difference between test data and fitting function
    def bilinear_error(self, p):
        return self.test_P - self.bilinear(p[0], p[1])

    # returns the parameters of the fitted curve
    def bilinear_params(self):
        iparams = [1.e6, 100.]
        return leastsq(self.bilinear_error, iparams)


    ###################################
    ###### three params bond law ######
    ###################################

    # pull-out function with three bond law params
    # params: friction, max stress, bond shear stiffness
    def three_param_func (self, k, qf, qy):
        w = sqrt(k / self.E_mod / self.A)
        mask = self.test_u > qy / k
        test_u_deb = self.test_u[mask]
        P_interp = interp1d(array([0., qy / k]), array([0., qf / w]))
        mask2 = self.test_u <= qy / k
        test_u_lin = self.test_u[mask2]
        P_lin = P_interp(test_u_lin)
        P_deb = sqrt((2. * test_u_deb * w ** 2 * self.E_mod * self.A * qf - 2. * qf * qy + qf ** 2)) / w
        return hstack((P_lin, P_deb))

    # error function - difference between test data and fitting function
    def three_error(self, p):
        return self.test_P - self.three_param_func(p[0], p[1], p[2])

    # returns the parameters of the fitted curve
    def three_params(self):
        iparams = [23.e4, 150., 200.]
        return leastsq(self.three_error, iparams)

    ###########################################
    ################ plotting #################
    ###########################################
    
    def plot(self, P):
        plt.plot(self.test_u, self.test_P, color = 'black', label = 'test data')
        plt.plot(pull.test_u, P, color = 'red', label = 'fitted curve')
        plt.xlabel('displacement [m]')
        plt.ylabel('force [N]')
        plt.legend(loc = 'best')

if __name__ == '__main__':

    k_list = []
    i = 0
    for ii in linspace(0, 8, 9):
        i += 2
        po = PullOutParams()
        # import values from *.csv file
        values = transpose(loadtxt('Carbon 1 bar 1mm 1 Prozent.csv', delimiter = ';', skiprows = 4, \
                     usecols = (arange(20))))
        
        # strains of a perfectly clamped yarn P/(E*A) - forces from csv are in cN
        # therefore the eps is divided by 100
        eps_yarn = values[i + 1] / po.E_mod / po.A / 100.
        
        # strains from pull-out of the clamp
        # strains from csv are in percent therefore the eps is divided by 100
        # eps are the strains without the strains of the yarn free length
        eps = values[i] / 100.
        u_pull = (eps - eps_yarn) * 0.47
        P_pull = values[i + 1] / 100.

        # there is still the successive stiffness increase due to filament activation
        # this part is cut off by finding the inflection point and linearizing
        # the function from that point to zero force
        # first the test data are smoothed by a polynom of the order k

        def polynom():
            p = polyfit(u_pull, P_pull, 20)
            return poly1d(p)
           
        def poly_der():
            return polyder(polynom(), 1)
        
        
        # linearizing the curve from the first inflection point to zero force 
        # direction the second diff spline, prolonging it with a linear function
        # and shifting the whole curve to [0,0]
        idx = argmax(poly_der()(u_pull))
        u_inf = u_pull[idx]
        P_inf = P_pull[idx]
        x = u_inf - P_inf / poly_der()(u_pull)[idx]
        mask = u_pull > u_pull[idx]
        u = u_pull[mask]
        P = polynom()(u_pull)[mask]
        # number of points from the test data in the linear part
        points = len(u_pull) - len(u)
        xx = linspace(x, u_inf, points)
        yy = linspace(0, P_inf, points)
        uu = hstack((xx, u))
        uu -= uu[0]
        PP = hstack((yy, P))
        # constructing the PullOut instance and initializing the u and P 
        pull = PullOutParams(test_u = uu, test_P = PP)
        
        def constant():
            # extracting the bond shear stiffness
            # and friction of the bilinear bond law
            qf = pull.const_params()[0]
            # fitted curve with evaluated parameters
            P = pull.const(qf)
            
            ######################################
            ########## PLOTTING CONSTANT #########
            ######################################
            pull.plot(P)
            plt.title('pull out with constant frictional bond law')
            ytick = min(pull.test_P) + 0.05 * (max(pull.test_P) - min(pull.test_P))
            xtick = min(pull.test_u) + 0.7 * max(pull.test_u) - min(pull.test_u)
            plt.text(xtick, ytick, 'qf = %s' % (qf))
            plt.show()
    
        def bilinear():
            # extracting the bond shear stiffness
            # and friction of the bilinear bond law
            k, qf = pull.bilinear_params()[0]
            # fitted curve with evaluated parameters
            P = pull.bilinear(k, qf)
            k_list.append(k)
    
            ######################################
            ########## PLOTTING BILINEAR #########
            ######################################
    
            pull.plot(P)
            plt.title('pull out with bilinear bond law')
            ytick_qf = min(pull.test_P) + 0.05 * (max(pull.test_P) - min(pull.test_P))
            ytick_k = min(pull.test_P) + 0.12 * (max(pull.test_P) - min(pull.test_P))
            xtick = min(pull.test_u) + 0.6 * max(pull.test_u) - min(pull.test_u)
            plt.text(xtick, ytick_qf, 'qf = %s' % (qf))
            plt.text(xtick, ytick_k, 'k = %s' % (k))
            plt.show()
    
    
        def three():
            # extracting the bond shear stiffness, friction and
            # max bond stress of the three params bond law
            k, qf, qy = pull.three_params()[0]
            # fitted curve with evaluated parameters
            P = pull.three_param_func(k, qf, qy)
            k_list.append(k)
    
            ##########################################
            ########## PLOTTING THREE PARAMS #########
            ##########################################
    
            pull.plot(P)
            plt.title('pull out with three params bond law')
            ytick_k = min(pull.test_P) + 0.2 * (max(pull.test_P) - min(pull.test_P))
            ytick_qf = min(pull.test_P) + 0.12 * (max(pull.test_P) - min(pull.test_P))
            ytick_qy = min(pull.test_P) + 0.05 * (max(pull.test_P) - min(pull.test_P))
            xtick = min(pull.test_u) + 0.6 * max(pull.test_u) - min(pull.test_u)
            plt.text(xtick, ytick_qf, 'qf = %s' % (qf))
            plt.text(xtick, ytick_k, 'k = %s' % (k))
            plt.text(xtick, ytick_qy, 'qy = %s' % (qy))
            plt.show()
        
        def run():
            
            #constant()
            bilinear()
            #three()
        run()
    
    plt.hist(k_list, 7)
    plt.show() 
    print mean(k_list), sqrt(var(k_list)), k_list
