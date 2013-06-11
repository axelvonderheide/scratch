'''
Created on 28.03.2012

@author: Axel
'''
from test_cb_resp_func import TESTCBShortFiber
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from enthought.traits.api import \
    HasTraits, Instance, on_trait_change, Int, Array, Tuple, List

from enthought.traits.api import HasTraits, Float, Property, \
                                cached_property, Range, Button
from enthought.traits.ui.api import View, Item, Tabbed, VGroup, \
                                VSplit, Group

from math import pi as Pi

import numpy as np

class Con_test(HasTraits):
    
    sim = Int(20, auto_set=False, enter_set=True,
               desc='simulations of one fiber', modified=True, param=True)
    
    nof = Int(100, auto_set=False, enter_set=True,
               desc='Number of fibers in crack', modified=True, param=True)
    
    E_f = Float(200e3, auto_set=False, enter_set=True,
               desc='steel modulus of elasticity [N/mm^2]', modified=True, param=True)
    
    E_m = Float(30e3, auto_set=False, enter_set=True,
               desc='matrix modulus of elasticity [N/mm^2]', modified=True, param=True)


    r_f = Float(0.075, auto_set=False, enter_set=True,
             desc='fiber radius[mm]', modified=True, param=True)
    

    f_f = Float(1., auto_set=False, enter_set=True, # [mm]
                 desc='snubbing coefficient', modified=True)
    
    V_f = Float(3., auto_set=False, enter_set=True,
             desc='volume fraction of steel fibers [%]', modified=True, param=True)

    f_length = Float(9., desc='in mm')


    A = Property(depends_on='height,width')
    @cached_property
    def _get_A(self):
        return self.height * self.width


    Ar = Property(depends_on='Vf')
    @cached_property
    def _get_Ar(self):
        return self.A * self.Vf / 100.

    Am = Property(depends_on='A_r, V_f')
    @cached_property
    def _get_Am(self):
        return self.A_m - self.A_r

    Kr = Property(depends_on='V_f, E_r')
    @cached_property
    def _get_Kr(self):
        return self.A_r * self.E_r

    Km = Property(depends_on='A_m, V_f, E_m')
    @cached_property
    def _get_Km(self):
        return self.A_m * self.E_m

    Kc = Property(depends_on='A_m, V_f, E_r, E_m')
    @cached_property
    def _get_Kc(self):
        return self.K_r + self.K_m
    
    nx = Float(200)

    w = Property(Array, depends_on='f_length')
    @cached_property
    def _get_w(self):
        '''discretizes crack-opening'''
        return np.linspace(0, self.f_length / 2, self.nx)
    
    def attr_fiber(self):
        #rand(nof)* definition range
        phi_arr = np.random.rand(self.nof) * Pi / 2
        le_arr = np.random.rand(self.nof) * self.f_length / 2
        phi_sh_arr = phi_arr.reshape(len(phi_arr), 1)
        le_sh_arr = le_arr.reshape(len(le_arr), 1)
        return phi_sh_arr, le_sh_arr
    
    def rigid_matrix_test(self):
        #sum of fiber-forces
        
        Cbsf = CBShortFiber()
        #fiber attr
        phi_arr_r = np.random.rand(self.nof) * Pi / 2
        le_arr_r = np.random.rand(self.nof) * self.f_length / 2
        phi_arr_r = phi_arr_r.reshape(len(phi_arr_r), 1)
        le_arr_r = le_arr_r.reshape(len(le_arr_r), 1)
        ####1#####
        resp_arr1 = Cbsf(self.w, 2, self.f_length, 2 * self.r_f, self.E_f, le_arr_r, phi_arr_r, self.f_f, 0, 0, 1e15)
        sum1 = np.sum(resp_arr1, axis=0) 
        ####2#####        add. simulations        ######
        resp_arr2 = Cbsf(self.w, 2, self.f_length, self.sim * 2 * self.r_f, self.E_f, le_arr_r, phi_arr_r, self.f_f, 0, 0, 1e15) / self.sim
        sum2 = np.sum(resp_arr2, axis=0) 
           
        return sum1, sum2
    
    
    
    def lin_ela_matr_test(self):
        #Input : Af,Ef,p,phi,le, A_m
        Cbsf = TESTCBShortFiber()
        ##############fiber attr#######################
        phi_arr = np.random.rand(10000) * Pi / 2
        le_arr = np.random.rand(10000) * self.f_length / 2
        phi_arr = phi_arr.reshape(len(phi_arr), 1)
        le_arr = le_arr.reshape(len(le_arr), 1)
        ###############################################
        ####1#####
        #################### w, tau, L_f, D_f, E_f, L_e, phi, f , nu
        resp_arr1 = Cbsf(self.w, 2, self.f_length, 2 * self.r_f, self.E_f, le_arr, phi_arr, self.f_f, 1)
        sum1 = np.sum(resp_arr1, axis=0) 
        ####2#####        add. simulations        ######
        resp_arr2 = Cbsf(self.w, 2, self.f_length, self.sim * 2 * self.r_f, self.E_f, le_arr, phi_arr, self.f_f, 1) / self.sim
        sum2 = np.sum(resp_arr2, axis=0) 
           
        return sum1, sum2
        
        
if __name__ == '__main__':
    from matplotlib import pyplot as plt
    a = Con_test()
    er1, er2 = a.rigid_matrix_test()
    ee1, ee2 = a.lin_ela_matr_test()
    plt.plot(a.w, er1, 'red')
    plt.plot(a.w, er2, 'blue')
    plt.plot(a.w, ee1, 'black')
    plt.plot(a.w, ee2, 'green')
    plt.show()
