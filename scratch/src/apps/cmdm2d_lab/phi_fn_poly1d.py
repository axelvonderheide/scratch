'''
Created on Sep 23, 2010

@author: alexander
'''
from enthought.traits.api import \
    Float, Instance, Array, Int, Property, cached_property, on_trait_change, Bool, \
    HasTraits, File, Event

from enthought.traits.ui.api import \
    View, Item, FileEditor, HSplit, Group, VSplit, CancelButton, OKButton, \
    Handler

from enthought.traits.ui.menu import \
    Action, CloseAction, HelpAction, Menu, \
    MenuBar, NoButtons, Separator, ToolBar

from enthought.pyface.api import ImageResource

from ibvpy.mats.mats_explore import MATSExplore
from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore

from numpy import copy, array, hstack, loadtxt, savetxt, poly1d, polyfit



from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from scipy.optimize import brentq, newton, fsolve, brenth
from os.path import join
from ibvpy.core.tloop import TLoop, TLine
from ibvpy.core.scontext import SContext
from ibvpy.core.tstepper import TStepper

from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import \
    MATS2DMicroplaneDamage, PhiFnGeneral, PhiFnStrainSoftening

from promod.exdb.ex_run import ExRun

data_file_editor = FileEditor( filter = ['*.DAT'] )

from traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
import pickle
from copy import copy

from promod.simdb import SimDB
simdb = SimDB()


matdata_dir = simdb.matdata_dir
ccs_unit_cell_dir = "CCSUnitCell"
file_name = 'PZ-0708-1_MAG-07-03_0.00300_all90.pickle'

pickle_file_path = join( matdata_dir, ccs_unit_cell_dir, file_name )
file = open( pickle_file_path, 'r' )
ccsuc = pickle.load( file )
df_list = ccsuc.damage_function_list
print 'df_list[0].damage_function.xdata', df_list[0].damage_function.xdata

xdata = df_list[0].damage_function.xdata
ydata = df_list[0].damage_function.ydata

phi_fn_poly = polyfit( xdata, ydata, 2 )
print 'phi_fn_poly', phi_fn_poly

#mats_eval = MATS2DMicroplaneDamage( 
#                                n_mp = 30,
#                                elastic_debug = False,
#                                stress_state = 'plane_stress',
#                                symmetrization = 'sum-type',
#                                model_version = 'compliance',
#                                phi_fn = PhiFnGeneral,
#                                )
#
##df = ccsuc.get_param( mats_eval, 'TT06-9u-V1.pickle' )
##df = ccsuc.get_param( mats_eval, 'TT06-9u-V2-all90.pickle' )
#df = ccsuc.get_param( mats_eval, 'TT06-9u-V3-all90.pickle' )
#print df
#
#phi_fn_poly = poly1d( xdata, ydata, 3 )
#print phi_fn_poly


#TT06-9u-V1.pickle
#TT06-9u-V2-all90.pickle
#TT06-9u-V3-all90.pickle

#xdate, ydata = trace.xdata, trace.ydata
#
#print xdata, xdata

#
#
#load( matdata / CCSUnitCell / PZ - 0708 - 1_MAG - 07 - 03_0.00300_all90.pickle



#    composite_tensile_test = Property
#    def _get_composite_tensile_test( self ):
#        return self.ex_run.ex_type
#
#    def get_target_data_exdb_tensile_test( self ):
#    '''Use the data from the ExDB
#    '''
#        ctt = self.composite_tensile_test
#        return ctt.eps_smooth, ctt.sig_c_smooth
#
#
#
#
#    mfn_line_array_target = Property( Instance( MFnLineArray ),
#                                      depends_on = 'ex_run' )
#    @cached_property
#    def _get_mfn_line_array_target( self ):
#        xdata, ydata = self.get_target_data_exdb_tensile_test()
#        return MFnLineArray( xdata = xdata, ydata = ydata )
