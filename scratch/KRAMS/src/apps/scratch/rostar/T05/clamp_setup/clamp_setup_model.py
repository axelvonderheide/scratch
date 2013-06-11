'''
Created on May 11, 2010

@author: rostislav
'''
from enthought.traits.api import HasTraits, Float, Property, Enum
from enthought.traits.ui.api import View, Item
from matplotlib import pyplot as plt
from numpy import linspace, hstack, sqrt, tanh, sinh, cosh, array
from scipy.optimize import brentq

''' module for plotting the strain, displacement and stress in a yarn held by two clamps '''



class ClampSetup(HasTraits):
    
    switch = Enum('displacement', 'strains', 'slip in clamp', modified = True)
    l1 = Float(0.025, modified = True, auto_set = False, enter_set = True)
    lt = Float(0.025, modified = True, auto_set = False, enter_set = True)
    L = Float(0.05, modified = True, auto_set = False, enter_set = True)
    k = Float(30000000., modified = True, auto_set = False, enter_set = True)
    E = Float(72.e9, modified = True, auto_set = False, enter_set = True)
    A = Float(0.89e-6, modified = True, auto_set = False, enter_set = True)
    Pt = Float(55., modified = True, auto_set = False, enter_set = True)
    Ps = Float(30., modified = True, auto_set = False, enter_set = True)
    
    w = Property(depends_on = 'E, A, k')
    def _get_w(self):
        return sqrt(self.k / self.E / self.A)
 
    P1 = Property(depends_on = 'Pt, Ps')
    def _get_P1(self):
        if self.Pt - self.Ps < 0.:
            print 'additional force Ps must not be larger than the applied test force Pt'
        else:
            return self.Pt - self.Ps

    def L1_residuum(self, Lt):
        P1 = self.P1
        Pt = self.Pt
        E = self.E
        A = self.A
        w = self.w
        L = self.L
        et = Pt / (E * A) / cosh(w * Lt)
        e1 = P1 / (E * A) / cosh(w * (L - Lt))
        return e1 - et

    def get_values(self):
        P1 = self.P1
        Pt = self.Pt
        E = self.E
        A = self.A
        w = self.w
        L = self.L
        l1 = self.l1
        lt = self.lt
        
        xdata = linspace(0., L, 200)
        e = Pt / (E * A) * cosh(w * xdata) / cosh(w * L)
        if e[0] * A * E > P1:
            ec = (Pt * w * sinh(w * xdata) - P1 * w * sinh(w * (xdata - L))) / (E * A * w * sinh(w * L))
            uc = (Pt * cosh(w * xdata) - P1 * cosh(w * (xdata - L))) / (E * A * w * sinh(w * L))
            uc = uc - uc[0] + P1 * l1 / A / E
            xe = hstack((0., xdata + l1, l1 + L + lt))
            xu = hstack((0., xdata + l1, l1 + L + lt))
            e = hstack((P1 / A / E, ec, Pt / A / E))
            u = hstack((0., uc, uc[-1] + Pt * lt / A / E))
            
        else:
            Lt = brentq(self.L1_residuum, 0, L)
            xdata1 = linspace(0, (L - Lt), 100)
            xdata2 = linspace(0, Lt, 100)
            e1 = P1 / (E * A) * cosh(w * xdata1) / cosh(w * (L - Lt))
            e1 = e1[::-1]
            e2 = Pt / (E * A) * cosh(w * xdata2) / cosh(w * Lt)
            e = hstack((P1 / A / E, e1, e2, Pt / A / E))
            
            u1 = P1 / (E * A * w) * sinh(w * xdata1) / cosh(w * (L - Lt))
            u1 = -u1[::-1]
            u1 = u1 - u1[0]
            u2 = Pt / (E * A * w) * sinh(w * xdata2) / cosh(w * Lt)
            u_left = hstack((0., u1 + P1 * l1 / A / E))
            u_right = hstack((u2 - u2[-1], Pt * lt / A / E))
            u = hstack((u_left, u_right + u_left[-1] - u_right[0]))
            
            xu = hstack((0., xdata1 + l1, xdata2 + xdata1[-1] + l1, l1 + L + lt))
            xe = hstack((0., hstack((xdata1 - xdata1[-1], xdata2)) + L - Lt + self.l1, self.l1 + L + self.lt))
        return xu, u, xe, e
    
    traits_view = View(Item('l1', label = 'fixation to clamping distance [m]'),
                       Item('lt', label = 'tested length [m]'),
                       Item('E', label = 'yarn elasticity modulus [N/m]'),
                       Item('A', label = 'yarn cross-section [m2]'),
                       Item('L', label = 'clamping length [m]'),
                       Item('k', label = 'shear stiffness [N/m]'),
                       Item('Pt', label = 'applied force [N]'),
                       Item('Ps', label = 'additional force [N]'),
                       Item('switch'),
                       ) 

if __name__ == '__main__':
    cs = ClampSetup()
    val = cs.get_values()
    plt.plot(val[0], val[1])
    plt.plot(val[2], val[3])
    plt.show()
