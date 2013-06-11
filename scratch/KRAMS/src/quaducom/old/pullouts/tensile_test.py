'''
Created on Mar 29, 2010

@author: rostislav
'''
""" defining the pullout model as a traited class """

from math import sqrt
from numpy import linspace, frompyfunc, array, sin, sign

from enthought.traits.api import HasTraits, Instance, Range, Interface,\
                                Array, on_trait_change, Property,\
                                cached_property, Bool, Float, Enum, DelegatesTo,\
                                on_trait_change, List, implements, Event
from enthought.traits.ui.api import View, Item, HSplit, VSplit, Group,\
                                HGroup, VGroup, Spring

from enthought.traits.ui.menu import OKButton
from parameters import Geometry, Material

def Heaviside(x):
    return sign(sign(x)+1)

class TensileTest( HasTraits ):
    '''Response of an elastic brittle filament with
    slack and delayed activation.
    '''
    
    geometry = Instance(Geometry)
    def _geometry_default(self):
        return Geometry()

    material = Instance(Material)
    def _material_default(self):
        return Material()
    
    lambd = DelegatesTo('geometry')
    
    theta = DelegatesTo('geometry')
    
    Ef = DelegatesTo('material')
    
    fu = DelegatesTo('material')
    
    Af = DelegatesTo('geometry')
    
    l = DelegatesTo('geometry')
    
    rf = DelegatesTo('geometry')

    def get_value(self):
        Af = self.Af
        Ef = self.Ef
        l = self.l
        lambd = self.lambd
        theta = self.theta 
        fu = self.fu
        # max displacement to plot
        u_max = 1.2*l*(fu / Ef * (1 + lambd) * (1 + theta) + theta*(1 + lambd))
        u_l = linspace(0,u_max,100)

        return u_l, Ef*Af*(u_l/l - theta*(1+lambd))/((1+theta)*(1+lambd))*\
            Heaviside(u_l/l - theta*(1+lambd))*\
            Heaviside(fu/Ef-(u_l/l - theta*(1+lambd))/((1+theta)*(1+lambd)))
            
    traits_view = View(Item('Ef'),
                       Item('fu'),
                       Item('l'),
                       Item('rf'),
                       Item('Af', style = 'readonly'),
                       Item('lambd'),
                       Item('theta'),
                       )

if __name__=="__main__":
    tt = TensileTest()
    tt.configure_traits()