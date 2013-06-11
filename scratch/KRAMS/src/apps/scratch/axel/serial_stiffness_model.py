
'''
Created on 20.06.2011

@author: axel
'''
from enthought.traits.api import HasTraits, Float, Property, cached_property, \
    Range, Button, Instance, Enum, Bool, on_trait_change, Int, Event, Array, Tuple, \
    List
from enthought.traits.ui.api import  Tabbed, Group, VGroup, VSplit, View, Item#, \
   # View, Item, VGroup, HGroup, ModelView, HSplit, VSplit 
from enthought.traits.ui.menu import  OKButton#, CancelButton
from math import pi as Pi
from numpy import linspace, max, abs, array , sign, argmax #,\
    #mean, ones, round, zeros, diff, trapz ,frompyfunc ,hstack
from numpy import mean, reshape , sum as sumi
from numpy.random import  rand , randn
#from quaducom.crackbridge.crack_bridge import StressInFiberWithConstantFriction
from quaducom.ctt.scm_cuypers.analyt_scm_model import SCM
#from quaducom.ctt.scm_cuypers.reinf_cross_section import SimplyRatio, \
#   GridReinforcement
from quaducom.resp_func.cb_short_fiber import CBShortFiber
#from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
#from scipy.special import gamma
from scipy.stats import  norm #weibull_min
from stats.pdistrib.sin2x_distr import sin2x
from stats.random_field.gauss_1D import GaussRandomField
#from stats.spirrid.spirrid_nd import SPIRRID
#from traits.either_type import EitherType
from scipy.interpolate import interp1d
from numpy import trapz

# E MODUL 

def H( x ):
    return sign( sign( x ) + 1. )

class SerialStiffnessModel( SCM ):
    '''Global to be put in SCM'''


    length = Float( 1000.,
                    auto_set = False, enter_set = True, desc = 'total specimen length in [mm]', modified = True )

    l_rho = Float( 3., auto_set = False, enter_set = True, # [mm]
                 desc = 'autocorrelation length', modified = True )
    cracks = List


    nx = Int( 300, auto_set = False, enter_set = True,
                 desc = 'number of length discretization points', modified = True )

    applied_stress = Range( low = 1e-10, high = 50.0, value = 1.2,
                         auto_set = False, enter_set = True,
                 desc = 'current stress in [Mpa]', ctrl_param = True )
    result_arr = Tuple( Array, Array, Array )
    current_stress = Float
    E_f = 200e3



    x_arr = Property( Array, depends_on = 'length,nx' )
    @cached_property
    def _get_x_arr( self ):
        '''discretizes the specimen length'''
        return linspace( 0, self.length, self.nx )





    '''needs to be worked out'''

    sigma_m_ff = Property( depends_on = '+ctrl_param,+modified, reinf_ratio.+modified' )
    @cached_property
    def _get_sigma_m_ff( self ):
        '''stress in the matrix in an uncracked composite'''
        return self.applied_stress * self.E_m / self.E_c







    '''Steel fiber multiple cracking'''

    #Fiber values
    f = Float( 0.87, auto_set = False, enter_set = True, # [mm]
                 desc = 'snubbing coefficient', modified = True )
    fiber_length = Float( 17., desc = 'in mm' )
    fiber_volume_fraction = Float( 3., auto_set = False, enter_set = True,
                 desc = 'Volume Fraction', modified = True )

    r = 0.075
    tau = 1.67
    r = 0.075
    rho = 0.03

    #specimen dimensions
    height = Float ( 100., desc = 'total specimen height [mm]',
                     auto_set = False, enter_set = True, modified = True )
    width = Float( 100. , desc = 'total specimen width [mm]',
                   auto_set = False, enter_set = True, modified = True )

    #Crack bridge
    specimen_broken = Bool( False )
    cbs_pulloutlist = List
    weakest_cb = Float( 1e15 )
    n_of_f = Tuple( desc = 'No of fibers in specimen' )
    n_of_f = [[], [], []]

    w = Property ( Array , depends_on = 'fiber_length,tau,E_f,r' )
    @cached_property
    def _get_w( self ):
        return linspace( 0, 1.05 * self.fiber_length ** 2. * self.tau / ( self.E_f * self.r * 2 ), 100 )


    distr_of_fibers = Property ( Tuple )
    @cached_property
    def _get_distr_of_fibers( self ):
        #Giving mean,stdev of Fibers
        specimen_volume = self.length * self.width * self.height
        no_of_fibers_in_specimen = ( specimen_volume * self.fiber_volume_fraction / 100 ) / ( Pi * self.r ** 2 * self.fiber_length )
        prob_crackbridging_fiber = .5 * self.fiber_length / self.length
        mean = prob_crackbridging_fiber * no_of_fibers_in_specimen
        stdev = ( prob_crackbridging_fiber * no_of_fibers_in_specimen * ( 1 - prob_crackbridging_fiber ) ) ** 0.5
        return [mean, stdev]

    le_array = List

    distr_ppfs = Property( List )
    def _get_distr_ppfs( self ):
        no_of_fibers = norm( self.distr_of_fibers[0], self.distr_of_fibers[1] ).ppf( rand( 1 ) )
        phi_array = sin2x._ppf( rand( no_of_fibers ) ).reshape( no_of_fibers, 1 )
        le_array = ( rand( no_of_fibers ) * self.fiber_length / 2. )
        self.n_of_f[0].append( int( no_of_fibers ) )
        self.n_of_f[2].append( le_array )
        le_array = le_array.reshape( no_of_fibers , 1 )
        return phi_array, le_array

    all_resp_list = List
    invers_list = List
    def cbs_refresh( self ):
        isf = CBShortFiber()
        phi_array, le_array = self.distr_ppfs
        w = self.w.reshape( 1, 100 )
        all_resp = isf( w, self.tau, self.fiber_length, 2 * self.r, self.E_f, le_array, phi_array, self.f, 0, 0, 1e15 )
        resp = sumi( all_resp, axis = 0 )
        resp = resp[0:argmax( resp ) + 2]
        if resp.max() <= self.weakest_cb:
            self.weakest_cb = resp.max()
        self.cbs_pulloutlist.append( resp )
        function = interp1d( resp, self.w[0:argmax( resp ) + 2] )
        self.invers_list.append( function )
        self.all_resp_list.append( all_resp )





    u_list = List
    u = Property( Float )
    def _get_u( self ):
        if len( self.invers_list ) > 0:
            U = self.applied_stress * self.uncracked_length / self.E_c
            for i in range( len( self.invers_list ) ):
                U += self.invers_list[i]( self.applied_stress * self.width * self.height )
            return U
        else:
            return self.applied_stress * self.uncracked_length / self.E_c



    uncracked_length = Property( Float )
    def _get_uncracked_length( self ):
        sigma_x = array( len( self.x_arr ) * [self.sigma_m_ff] )
        return self.length * sum( self.sigma_m_x == sigma_x ) / self.nx



    applied_stress_list = List
    global_stiffness = List
    @on_trait_change( 'applied_stress' )
    def global_stiffness_evaluate( self ):
        #k_cbs = self.E_c * self.height * self.width / self.uncracked_length
        if self.applied_stress * self.width * self.height <= self.weakest_cb:
            '''for i in range ( len( self.cbs_pulloutlist ) ):
                index = abs( self.cbs_pulloutlist[i] - self.applied_stress * self.width * self.height ).argmin()
                k_cb = ( self.cbs_pulloutlist[i][ index + 1 ] - self.cbs_pulloutlist[i][index - 1] ) / ( self.w[2] - self.w[0] )
                print 'Steifigkeit', k_cb , 'Kraft', self.applied_stress * self.width * self.height
                k_cbs = 1 / ( 1 / k_cbs + 1 / k_cb )'''
            self.applied_stress_list.append( self.applied_stress )
            self.u_list.append( self.u )
        else: self.specimen_broken = True





    sigma_m_x = Property( Array, depends_on = 'cracks, +modified,cbs_refresh, +ctrl_param' )
    @cached_property
    def _get_sigma_m_x( self ):
        '''creates the matrix stress array given the current loading and an array of crack positions'''
        #self.cbs_refresh()
        sigma_x = array( len( self.x_arr ) * [self.sigma_m_ff] )
        if len( self.result_arr ) > 0:
            i = 0
            for crack in self.cracks:
                cr_arr = self.real_stress_curve( crack, self.n_of_f[2][i] )
                mask = sigma_x <= cr_arr
                sigma_x = sigma_x * mask + cr_arr * ( mask == False )
                i += 1
        return sigma_x



    def linear_stress_curve( self, crack_pos, no_of_fibers ):
        crack_x = self.x_arr - crack_pos
        return abs( self.tau * Pi * 2. * self.r * no_of_fibers / ( self.width * self.height ) * crack_x )

    def real_stress_curve ( self, crack_pos, le_array ):
        crack_x = ( self.x_arr - crack_pos )
        T = self.tau * Pi * 2. * self.r
        means = abs( self.active_fibers( crack_x , le_array ) * T * crack_x / ( self.width * self.height ) )
        return means

    def active_fibers ( self, crack_x, le_array ):
        crack_x = crack_x.reshape( len( crack_x ), 1 )
        matrix = abs( crack_x ) >= le_array
        y = sumi ( matrix, axis = 1 )
        return y




    #DON'T KNOW HOW TO UNWIRE THEM FROM STEEL FIBER METHODS

    @on_trait_change( '+modified, reinf_ratio.+modified' )
    def reset_history( self ):
        ''' if some params are changed, the sigma - eps history is not valid any more and is reseted'''
        self.cracks = []
        self.no_of_cracks = ( [0.], [0.] )
        self.sig_eps = ( [0.], [0.] )
        self.global_stiffness = []
        self.applied_stress_list = []
        self.crack_devel()
        self._get_distr_of_fibers()


    def launch( self, max, points ):
        self.reset_history()
        for stress in linspace( 0.01, max, points ):
            self.applied_stress = stress
            if self.applied_stress * self.width * self.height >= self.weakest_cb:
                break

    integ_list = List
    def crack_devel( self ):
        '''checks the matrix stress profile and compares it to the matrix CDF; adds new cracks'''
        x = self.x_arr

        # loading
        if self.sigma_m_ff > self.current_stress:
            self.current_stress = self.sigma_m_ff
            # find the points, where the strength is lower than the stress
            while sum( self.sigma_m_x >= self.random_field ) > 0.:
                cr_pos = argmax( self.sigma_m_x - self.random_field )
                self.cracks.append( self.x_arr[cr_pos] )
                self.cbs_refresh()

        # unloading           
        else:
            self.current_stress = self.sigma_m_ff

        eps_m = self.sigma_m_x / self.E_m
        eps_f = ( ( 1 + self.alpha ) * array( len( x ) * [self.applied_stress / self.E_c] ) - eps_m * self.alpha )
        self.integ_list.append( trapz( self.x_arr, eps_f ) )
        self.result_arr = ( x, eps_m , eps_f )


    no_of_cracks = Tuple( List, List )
    @on_trait_change( 'applied_stress' )
    def get_no_of_cracks( self ):
        self.no_of_cracks[0].append( self.current_stress )
        self.no_of_cracks[1].append( len( self.cracks ) )

    random_field = Property( Array, depends_on = 'l_rho, length, m, nsim,mean' )
    @cached_property
    def _get_random_field( self ):
        '''generates an array of random matrix strength'''
        rf = GaussRandomField( _l = self.l_rho, xgrid = self.x_arr , nsim = 1 , mean = 8 , stdev = 2.0 )
        return rf.random_field


    #Second vers.



    traits_view = View( 
                          Group( 
                            Tabbed( 
                              VGroup( 
                                   Group( Item( 'E_f', resizable = False,
                                               label = 'E-modulus',
                                               tooltip = "Young's modulus of the fiber" ),

                                         Item( 'sigma_fu', resizable = False,
                                               label = 'strength',
                                               help = "Strength of the fibers"
                                         ),
                                         label = 'fibers',
                                    ),

                                   Group( 
                                        Item( 'l_roh', resizable = False,
                                               label = 'autocorrelation length' ),
                                       Item( 'E_m', resizable = False,
                                             label = 'E-modulus',
                                             help = "Young's modulus of the matrix" ),
                                       Item( 'sigma_mu', resizable = False,
                                             label = 'strength',
                                             help = "Scale parameter of the matrix strength"
                                             'roughly corresponding to the mean strength' ),
                                       Item( 'm', resizable = False,
                                             label = 'Weibull-modulus',
                                             help = "Weibull modulus of the matrix strength distribution"
                                             'defining the scatter of the strength' ),
                                        label = 'matrix',
                                        scrollable = False,
                                        ),
                                   label = 'Components',
                                   dock = 'tab',
                                   id = 'scm.model.component_params',
                               ),
                                   VGroup( 
                                       Item( 'tau', resizable = False, springy = True ),
                                       Item( 'r', resizable = False, springy = False ),
                                       springy = True,
                                       label = 'Bond',
                                       dock = 'tab',
                                       id = 'scm.model.params',
                                    ),
                                    id = 'scm.model.allparams',
                                  ),
                               VGroup( 
                                   Item( 'reinf_ratio@', show_label = False, resizable = True ),
                                   label = 'Cross section parameters',
                                   dock = 'tab',
                                   id = 'scm.model.reinf_ratio',
                               ),
                               VGroup( 
                                    Item( 'orientation', label = 'fiber orientation' ),
                                    ),
                               id = 'scm.model.splitter',
                               springy = False,
                               layout = 'split',
                            ),
                            id = 'scm.model',
                            dock = 'fixed',
                            scrollable = True,
                            resizable = True,
                            buttons = [OKButton],
                            height = 0.8, width = 0.8
                                   )



