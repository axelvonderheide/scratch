'''
Created on Jun 23, 2010

@author: alexander
'''

from enthought.traits.api import \
    Float, Instance, Array, Int, Property, cached_property, on_trait_change, Bool, \
    HasTraits, File, Event

from enthought.traits.ui.api import \
    View, Item, FileEditor, HSplit, Group, VSplit, \
    Handler

from enthought.traits.ui.menu import \
    Action, CloseAction, HelpAction, Menu, \
    MenuBar, NoButtons, Separator, ToolBar

from enthought.pyface.api import ImageResource

from ibvpy.mats.mats_explore import MATSExplore
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore

from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like, \
                vstack, savetxt, hstack, argsort, fromstring, zeros_like, \
                copy, c_, newaxis, argmax, where, argsort, sqrt, frompyfunc, \
                max as ndmax, dot, fabs

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from scipy.optimize import brentq, newton, fsolve, brenth
from os.path import join
from ibvpy.core.tloop import TLoop, TLine
from ibvpy.core.scontext import SContext
from ibvpy.core.tstepper import TStepper

from promod.exdb.ex_run import ExRun

data_file_editor = FileEditor( filter = ['*.DAT'] )

from traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
import pickle
from copy import copy

from promod.simdb import SimDB
simdb = SimDB()

class SigFlCalib( HasTraits ):

    #--------------------------------------------------
    # Data source for calibration within simdb
    #--------------------------------------------------

    # @todo: get average value of all repetitions of the tensile test
    # @todo: check list of CCS fro consitency (must all be equal)
    # @todo: get ydata and ydata and calculate the average put in new
    #'mfn_line_array'

    ex_run = Instance( ExRun )

    composite_tensile_test = Property
    def _get_composite_tensile_test( self ):
        return self.ex_run.ex_type

    composite_cross_section = Property
    def _get_composite_cross_section( self ):
        return self.composite_tensile_test.ccs

    def get_data_tensile_test( self ):
        '''Use the data from the ExDB
        '''
        ctt = self.composite_tensile_test
        return ctt.eps_smooth, ctt.sig_c_smooth

    #---------------------------------------------------------------
    # PLOT OBJECT
    #-------------------------------------------------------------------
    figure = Instance( Figure )
    def _figure_default( self ):
        figure = Figure( facecolor = 'white' )
        figure.add_axes( [0.12, 0.13, 0.85, 0.74] )
        return figure


    #---------------------------------------------------------------
    # material properties textile reinforcement
    #-------------------------------------------------------------------

    # security factor 'gamma_tex' for the textile reinforcement
    #
#    gamma_tex = 1.5

    # reduction factor to drive the characteristic value from mean value
    # of the experiment (EN DIN 1990)
    #
#    beta = 0.81

    #---------------------------------------------------------------
    # material properties concrete
    #-------------------------------------------------------------------

    # characteristic compressive force [MPa]
    #
    f_ck = 60.0

    # concrete strength after 9 days
    #
#    f_ck = 49.0
    print 'f_ck', f_ck

    # security factor 'gamma_c' for hight strength concrete (>C55/67)
    #
#    gamma_c = 1.5 * ( 1 / ( 1.1 - f_ck / 500. ) )
#    print 'gamma_c', gamma_c

    # design value of the compressive force [MPa]
    #
#    f_cd = 0.85 * f_ck / gamma_c
#    print 'f_cd', f_cd

    # factor [-] to calculate the value of the resulting compressive
    # force, i.e. f_c = chi * fck / gamma_c
    # (for high strength concrete)
    #
    chi = 1.05 - f_ck / 500.
    print 'chi', chi

    # factor [-] to calculate the distance of the resulting compressive
    # force from the top, i.e. a = k * x
    #
    k = 1.05 - f_ck / 250.
    print 'k', k

    # compressive strain at the top at rupture [-]
    # value measured at the experiment (with DMS)
    # for 0: about 3.3
    # for 90: about 3.0
    # ultimute strain theoreticaly (Brockman): about 4.5
    # NOTE: strain was meassured at a distance of 5 cm
    #
    eps_c = 3.3 / 1000. # 3.3 promile
    print 'eps_c', eps_c

    #---------------------------------------------------------------
    # properties of the composite cross section
    #-------------------------------------------------------------------

    # height of the cross section [m]
    #
    thickness = 0.06
    print 'thickness', thickness

    # width of the cross section [m]
    #
    width = 0.20
    print 'width', width

    # total number of reinforcement layers [-]
    #
    n_layers = 12
    print 'n_layers', n_layers

    # spacing between the layers [m]
    #
    s_tex_z = thickness / ( n_layers + 1 )
    print 's_tex_z', s_tex_z

    # distance from the top of each reinforcement layer [m]:
    #
    z_i_arr = array( [ thickness - ( i + 1 ) * s_tex_z for i in range( n_layers ) ],
                     dtype = float )
    print 'z_i_arr', z_i_arr

    # iteration counter:
    #
    n = 1

    def get_lack_of_fit( self, u ):
        '''Return the difference between
        N_c (=compressive force of the compressive zone of the concrete)
        N_t (=total tensile force of the reinforcement layers)

        NOTE: eps_t (=tensile strain at the bottom [MPa]) is the search parameter
        to be found iteratively!
        '''
        print '\n'
        print '------------- iteration: %g ----------------------------' % ( self.n )

        eps_t, E_yarn = u
        # ------------------------------------
        # derived params depending on value for 'eps_t'
        # ------------------------------------

        # heights of the compressive zone:
        #
        x = abs( self.eps_c ) / ( abs( self.eps_c ) + abs( eps_t ) ) * self.thickness
        print 'x', x

        # strain at the height of each reinforcement layer [-]:
        #
        eps_i_arr = eps_t / ( self.thickness - x ) * ( self.z_i_arr - x )
        print 'eps_i_arr', eps_i_arr

        # use a ramp function to consider only positive strains
        eps_i_arr = ( fabs( eps_i_arr ) + eps_i_arr ) / 2.0
        print 'eps_i_arr', eps_i_arr

        # construct the constitutive law of the crack bridge - linear elastic 
        # with the search effective modulus of elasticity
        #
        xdata = array( [0., 1.] )
        ydata = array( [0., E_yarn ] )
        mfn_line_array_tensile_test = MFnLineArray( xdata = xdata, ydata = ydata )

        eps_u = 1.0
        print 'eps_u', eps_u

        # tensile stress at the height of each reinforcement layer [MPa]:
        #
        get_sig_i_arr = frompyfunc( mfn_line_array_tensile_test.get_value, 1, 1 )
        sig_i_arr = get_sig_i_arr( eps_i_arr )

        print 'sig_i_arr', sig_i_arr

        # tensile force of one reinforced composite layer [kN]:
        #
        n_rovings = 23   # for 0-direction (use 25 for 90-direction)
        A_roving = 0.461 #mm**2
        f_i_arr = sig_i_arr * n_rovings * A_roving / 1000.
        print 'f_i_arr', f_i_arr

        # sig_comp_i_arr [MPa]:
        #
        sig_comp_i_arr = f_i_arr / self.width / self.s_tex_z / 1000.
        print 'sig_comp_i_arr', sig_comp_i_arr

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (characteristic value)
        #
        N_tk = sum( f_i_arr )
        print 'N_tk', N_tk

        # total tensile force of all layers of the composite tensile zone [kN]:
        # (design value)
        #
#        N_td = N_tk / self.gamma_tex * self.beta
#        print 'N_td', N_td

        # distance [m] of the resulting compressive
        # force from the top, i.e. a = k * x / 2
        #
        a = self.k * x / 2.
        print 'a', a

#        k_exact = ( 1.74 * self.eps_c / 4.56 - ( self.eps_c / 4.56 ) ** 2 / ( 1 - 0.12 * self.eps_c / 4.56 ) )

        # total compressive force of the composite compressive zone [kN]:
        # (characteristic value)
        #
        N_ck = self.k * x * self.width * self.chi * self.f_ck * 1000.
        print 'N_ck', N_ck

        # total compressive force of the composite compressive zone [kN]:
        # (design value)
        #
#        N_cd = self.k * x * self.width * self.chi * self.f_cd * 1000.
#        print 'N_cd', N_cd

        # absolute error (equilibrium of sum N = 0):
        #
        d_N = N_ck - N_tk

        print 'd_N', d_N

        # resistance moment of one reinforced composite layer [kNm]:
        # (characteristic value)
        #

        print 'a * N_ck', a * N_ck
        print 'dot( f_i_arr, self.z_i_arr )', dot( f_i_arr, self.z_i_arr )
        M_Rk = dot( f_i_arr, self.z_i_arr ) - a * N_ck
        print 'M_Rk', M_Rk

        M_imposed = 3.49

        d_M = M_Rk - M_imposed

        self.n += 1

        return array( [ d_N, d_M ], dtype = float )

    def fit_response( self ):
        '''iterate 'eps_t' such that the lack of fit between the calculated
        normal forces in the tensile reinforcement and the compressive zone (concrete)
        is smaller then 'xtol' defined in function 'brentq'.
        NOTE: the method 'get_lack_of_fit' returns the relative error.
        '''

        # use scipy-functionality to get the iterated value of 'eps_t'
        # NOTE: get_lack_of_fit must have a sign change as a requirement
        # for the function call 'brentq' to work property.

        # The method brentq has optional arguments such as
        #   'xtol'    - absolut error (default value = 1.0e-12)
        #   'rtol'    - relative error (not supported at the time)
        #   'maxiter' - maximum numbers of iterations used
        #
        xtol = 1.0e-6

        #estimate does not converge!
#        u0 = array( [ 0.0001, 140000.0 ], dtype = float )

        # take as estimate for the solution for eps_t_new = 10 promile = 0.010!
        #
        u0 = array( [ 0.010, 140000.0 ], dtype = float )
        u_sol = fsolve( self.get_lack_of_fit, u0, xtol = xtol )

        # @todo: check if 'brenth' gives better fitting results; faster?
#            phi_new = brenth( self.get_lack_of_fit, 0., eps_t )
        eps_t_new, E_yarn = u_sol
        print 'eps_t_new', eps_t_new
        print 'E_yarn', E_yarn

if __name__ == '__main__':

    from promod.exdb.ex_run import ExRun
    import pylab as p

    #-------------
    # direction for 0 AND 90 direction:
    #-------------
    path = join( simdb.exdata_dir, 'tensile_tests',
                 'ZiE_2011-06-08_TT-12c-6cm-90-TU',
               )

    #-------------
    # tensile tests in 0-direction:
    #-------------
    tests = [ 'TT-12c-6cm-0-TU-V1.DAT']
#    tests = [ 'TT-12c-6cm-0-TU-V2.DAT']
#    tests = [ 'TT-12c-6cm-0-TU-V3.DAT']

    #-------------
    # tensile tests in 90-direction:
    #-------------
#    tests = [ 'TT-12c-6cm-90-TU-V1.DAT']
#    tests = [ 'TT-12c-6cm-90-TU-V2.DAT']
#    tests = [ 'TT-12c-6cm-90-TU-V3.DAT']

#    tests = [ 'TT-12c-6cm-90-TU-V1.DAT', 'TT-12c-6cm-90-TU-V2.DAT','TT-12c-6cm-90-TU-V3.DAT']

    for t in tests:
        ex_path = join( path, t )
        ex_run = ExRun( ex_path )
#        ex_run.ex_type._plot_c_smooth_stress_strain( p )
#        ex_run.ex_type._plot_c_smooth_stress_strain( p )

#            ex_run.ex_type._plot_force_edge_deflection( p )
#            ex_run.ex_type._plot_force_center_edge_deflection( p )

    # NOTE: works only if 'set_axes' in CTT is uncommented
    # @todo:
#    p.show()

    print '\n'
    print 'setup SigFlCalib'
    print '\n'

    sig_fl_calib = SigFlCalib( ex_run = ex_run,

                              )

#    print 'CTT', sig_fl_calib.composite_tensile_test
#    print 'CCS', sig_fl_calib.composite_tensile_test.ccs
#    print 'DTT', sig_fl_calib.get_data_tensile_test()
#    print 'DTT', sig_fl_calib.mfn_line_array_tensile_test.xdata

    sig_fl_calib.fit_response()
