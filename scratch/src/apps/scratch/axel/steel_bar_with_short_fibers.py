'''
Created on Aug 17, 2011

@author: rostar
'''
import numpy as np
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Range, Bool, on_trait_change, Int, Array, Tuple,Instance, \
    List
from enthought.traits.ui.api import  VGroup, View, Item
from math import pi as Pi
from quaducom.ctt import ICBM
from quaducom.ctt.homogenized_crack_bridges.steel_bar import SteelBar
#from quaducom.ctt.homogenized_crack_bridges.short_fibers_monte_carlo import ShortFibersMonteCarlo
from scipy.interpolate import interp1d
from apps.scratch.axel.short_fibers_monte_carlo import ShortFibersMonteCarlo

def H(x):
    return np.sign(np.sign(x) + 1.)

class SteelBarwithShortFibers(HasTraits):
    
    # applied force
    P = Float(modified = True) # [N]

    # closest crack from left
    Ll = Float(modified = True) # [mm]

    # closest crack from right
    Lr = Float(modified = True) # [mm]

    E_rod = Float(200e3, auto_set = False, enter_set = True,
               desc = 'steel modulus of elasticity [N/mm^2]', modified = True, param = True)
    
    E_fibers= Float(200e3, auto_set = False, enter_set = True,
               desc = 'steel modulus of elasticity [N/mm^2]', modified = True, param = True)
    
    sigma_max_rod= Float(500, auto_set = False, enter_set = True,
               desc = 'steel max [N/mm^2]', modified = True, param = True)

    Em = Float(30e3, auto_set = False, enter_set = True,
               desc = 'matrix modulus of elasticity [N/mm^2]', modified = True, param = True)

    tau_rod = Float(2.3, auto_set = False, enter_set = True,
               desc = 'sheer force per unit area [N/mm^2]', modified = True, param = True)
    
    tau_fibers = Float(2.3, auto_set = False, enter_set = True,
               desc = 'sheer force per unit area [N/mm^2]', modified = True, param = True)


    height = Float (30., desc = 'total specimen height [mm]',
                     auto_set = False, enter_set = True, modified = True)
    width = Float(30. , desc = 'total specimen width [mm]',
                   auto_set = False, enter_set = True, modified = True)

    Vf_fibers = Float(3., auto_set = False, enter_set = True,
             desc = 'volume fraction of steel fibers [%]', modified = True, param = True)

    r_fibers = Float(0.075, auto_set = False, enter_set = True,
             desc = 'fiber radius[mm]', modified = True, param = True)
    
    r_rod = Float(4, auto_set = False, enter_set = True,
             desc = 'rod radius[mm]', modified = True, param = True)

    f_fibers = Float(1., auto_set = False, enter_set = True, # [mm]
                 desc = 'snubbing coefficient', modified = True)
    
    
    discr = Float(400.)

    fiber_length = Float(9., desc = 'in mm')
    length = Float(400.)

    with_tau_distr = Bool(True , desc = 'enables Monte-Carlo simulation on every fiber"s tau', modified = True, param = True)
    weibull_tau_shape = Float(8.50666480278)
    weibull_tau_scale = Float(2.49019483074)
    with_f_distr = Bool(False , desc = 'enables Monte-Carlo simulation on every fiber"s snubbing coefficent')
    specimen_broken = Bool(False)

    Ac = Property(depends_on = 'height,width')
    @cached_property
    def _get_Ac(self):
        return self.height * self.width

    Ar = Property(depends_on = 'Vf')
    @cached_property
    def _get_Ar(self):
        return self.Ac * self.Vf / 100+self.r_rod**2*Pi

    Am = Property(depends_on = 'Ac, Vf')
    @cached_property
    def _get_Am(self):
        return self.Ac - self.Ar

    Kr = Property(depends_on = 'Vf, Er_rod')
    @cached_property
    def _get_Kr(self):
        return (self.Ar-self.Ac*self.Vf) * self.Er_rod+ self.Ac*self.Vf

    Km = Property(depends_on = 'Ac, Vf, Em')
    @cached_property
    def _get_Km(self):
        return self.Am * self.Em

    Kc = Property(depends_on = 'Ac, Vf, Er, Em')
    @cached_property
    def _get_Kc(self):
        return self.Kr + self.Km
    
    cbsfmc = Instance(ICBM)
    def _cbsfmc_default(self):
        return ShortFibersMonteCarlo( length = self.length, height= self.height , width=self.width )
    
    cbsb = Instance(ICBM)
    def _cbsb_default(self):
        return SteelBar(length = self.length)
    
    def get_eps_x_reinf (self, crack_x):
        return self.get_sigma_x_reinf(crack_x)
        
    def get_sigma_x_reinf(self, x):
        share_sf = self.share_of_force(x )
        self.cbsfmc.P=self.P *share_sf
        self.cbsb.P =self.P*(1-share_sf) 
        sigma_x_sfmc=self.cbsfmc.get_sigma_x_reinf(x)
        sigma_x_sb=self.cbsb.get_sigma_x_reinf(x)
        return sigma_x_sfmc
        
    
    def share_of_force( self , x ):
        crack_index=abs(x[0])
        self.cbsfmc.cbs_allocation(crack_index)
        func_sfmc=self.cbsfmc.invers_list[self.cbsfmc.crack_list.index(crack_index)]
        #adding both strength-crack_opening relationships.
        breaking_force_rod=self.sigma_max_rod*self.r_rod **2*Pi
        discr_P=np.linspace(0,breaking_force_rod/4,50)
        p_w_global=[]
        self.cbsb.Lr=self.Lr
        for i in discr_P:
                self.cbsb.P=i
                p= self.cbsb.P_w + func_sfmc(i)
                p_w_global.append(p)
        func_global = interp1d(discr_P,p_w_global)
        w_global = func_global(self.P)
        P_share_sfmc = self.cbsfmc.w_P_list[self.cbsfmc.crack_list.index(crack_index)](w_global)
        return P_share_sfmc/self.P
        
        
    traits_view = View(
                       VGroup(
                           Item('r', label = 'fiber radius', resizable = False, springy = True),
                           Item('Er', resizable = False, springy = False),
                           Item('Em', resizable = False, springy = False),
                           Item('tau', resizable = False, springy = False),
                           Item('height', label = 'specimen height', resizable = False, springy = False),
                           Item('width', label = 'specimen width', resizable = False, springy = False),
                           Item('f', label = 'snubbing coefficient f', resizable = False, springy = False),
                           Item('fiber_length', resizable = False, springy = False),
                           Item('with_tau_distr', label = 'with tau gauss distr', springy = False),
                           Item('weibull_tau_shape', label = 'mean tau', springy = False),
                           Item('weibull_tau_scale', label = 'Stdev tau', springy = False),
                           Item('with_f_distr', label = 'with f gauss distr', springy = False),
                           Item('mean_f', label = 'mean f', springy = False),
                           Item('stdev_f', label = 'Stdev f', springy = False),

                           springy = True,
                           label = 'CB parameters',
                           dock = 'tab',
                           id = 'cb.steel_bar.params',
                        ),
                            id = 'cb.steel_bar',
                            dock = 'fixed',
                            scrollable = True,
                            resizable = True,
                            height = 0.8, width = 0.8
                                   )
        
if __name__ == '__main__':
        from matplotlib import pyplot as plt
        sb = SteelBarwithShortFibers(Ll = 60., tau = 3.5 , P = 4000)
        #sb.configure_traits()
        x = np.linspace(0, 20, 200)
        #print x[35]
        x -= x[70]
    
        #print x
        #crack_x = x -
        #x = np.linspace( -10 , 10 , 100 )
    #    eps = sb.get_eps_x_reinf( x )
    #
    #    plt.plot( x, eps, lw = 2, color = 'black' )
    #    plt.show()
    
        import enthought.mayavi.mlab as m
        from stats.spirrid import orthogonalize
    
        resolution = 200
        profiles_list = []
        forces_arr = np.linspace(5, 620, resolution)
        for i, p in enumerate(forces_arr):
            sb.P = p
            profiles_list.append(sb.get_eps_x_reinf(x))
    
        profiles_arr = np.array(profiles_list)
    
        param_arr = orthogonalize([x, forces_arr])
        norm_param_arr = [ e / np.max(np.fabs(e)) for e in param_arr ]
    
        norm_profiles_arr = profiles_arr / np.max(np.fabs(profiles_arr))
        m.surf(norm_param_arr[0], norm_param_arr[1], norm_profiles_arr)
        m.show()
        
        
        
        
        
        
        
        
    
    
    
    
    