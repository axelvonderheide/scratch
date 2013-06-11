
from enthought.traits.api import \
    Float, HasTraits, Int, List, Array, Dict, String, WeakRef, Instance, Enum, Property, Any, Delegate, \
    cached_property, Constant, Button, on_trait_change, Bool

from ytta.spirrid.resp_func.resp_func import IResponseFunction 
from stats.pdistrib.pdistrib import IPDistrib, PDistrib
from ytta.spirrid.resp_mgr import ResponseManager
from ytta.spirrid.spirrid_integ import IIntegAlg, EquiDistant
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from enthought.traits.ui.api \
    import TabularEditor, View, Item, Group, HGroup, VGroup, HSplit, VSplit, Tabbed 
from enthought.traits.ui.menu \
    import OKButton, CancelButton
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from enthought.enable.component_editor import \
    ComponentEditor
from enthought.chaco.api import \
    Plot, AbstractPlotData, ArrayPlotData, \
    ArrayDataSource
from enthought.chaco.tools.api import \
    PanTool, ZoomTool
import math
from numpy import frompyfunc, trapz, linspace, \
    array, mgrid, ogrid, sqrt, max, pi,e, argmax, argmin, abs
import ytta.spirrid.resp_func as resp_func
from scipy.special import erf

#----------------------------------------------------------------------------------
#                                     RIDVariable
#----------------------------------------------------------------------------------
class RIDVariable(HasTraits):
    """
    Association between a random variable and distribution.
    """
    
    # should this variable be randomized
    
    randomized = Bool( False, randomization_changed = True )
    
    # name of the randomized variable (within the response function)
    #
    varname = String
    
    # --------------------------------------------

    # default view specification
    traits_view = View( HGroup( Item( 'varname', # , style = 'readonly', 
                                     show_label = False ),
                                     Item('randomized'),
                                     ),
                        resizable = True,
                        dock = 'tab',
                        id = 'rid_variable',
                        height=800 )

#-- Tabular Adapter Definition -------------------------------------------------

class RVAdapter ( TabularAdapter ):
    '''This adapter specialization lists the attributes of RIDVariable
    to be displayed in columns of the TabularEditor'''
    
    columns = [ ( 'Name',         'varname' ), 
                ( 'randomized',   'randomized' ) ]
                
    font                      = 'Courier 10'
    variable_alignment        = Constant( 'right' )
    
# -- Tabular Editor Definition -------------------------------------------------
# Edit the list of random variables
# 
rv_editor = TabularEditor(
    selected   = 'current_variable',
    adapter    = RVAdapter(),
    operations = [ 'move' ],
    auto_update = True
)


#from resp_func.fil_no_debonding import FilNoDebonding

#----------------------------------------------------------------------------------
#                                     SPIRRID                                     
#----------------------------------------------------------------------------------
class SPIRRID( HasTraits ):
    """
    Class for evaluating the Phoenix integral representing the strain-based
    fiber-bundle model.
    """
    # --------------------------------------------
    # Cached property for the response function
    # it depends on the response_function_type
    # so that it is constructed anew whenever
    # the type of the response function gets
    # changed.
    #
    response_function = Instance( IResponseFunction )
    def _response_function_default(self):
        return resp_func.FilStrengthDebonding()
                  
    variables = List()
    def _variables_default(self):
        return self.update_variables()
    
    def update_variables( self ):
        '''
        reset the RIDVariable list
        '''
        varset = []
        for key in ['Ef', 'Af', 'l', 'sigma_f','tau']:
            varset.append (
                    RIDVariable( varname = key ) )
        return varset

    # selected variable within the list of variables
    #
    current_variable = Instance( RIDVariable )
    def _current_variable_default(self):
        return self.variables[0]

    #-------------------------------------------------------
    # Methods performed prior to the evaluation
    #-------------------------------------------------------

    displ_force_stdev = Property( depends_on =  '+integration_changed'
                                                ',variables.+randomization_changed' )
#                                                ',variables.randomized' )
#                                                ',+response_function_changed,' )
    @cached_property
    def _get_displ_force_stdev(self):

        print 'displ_force_stdev_called'
        return 1
                

    # show the user interface
    traits_view = View( Group( VSplit( 
                                    HSplit(     
                                               VGroup(
                                                     Item('response_function',
                                                          style = 'custom',
                                                          show_label = False,
                                                          resizable = True), 
                                                          label = 'Response function',
                                                          dock = 'tab',
                                                          id = 'mine.spirrid.response_func'
                                                    ),
                                               VGroup(
                                                     Item('variables', show_label = False, editor = rv_editor,
                                                          resizable = True,
                                                          width = 200,
                                                          height = 200,
                                                          style = 'custom'
                                                          ),
                                                     label = 'Variables',
                                                     dock = 'tab',
                                                     id = 'mine.spirrid.variables'
                                                     ),                     
                                                id = 'mine.spirrid.resp_settings'
                                            ),
                                    HGroup(
                                           Item('current_variable', show_label = False, 
                                                style = 'custom', resizable = True),
                                                label = 'Random parameters',
                                                id = 'spirrid.current'
                                            ),
                                    id = 'mine.spirrid.resp_variables'       
                                    ),
                                label = 'Set params',
                                dock = 'tab',
                                id = 'spirrid.response+variables'
                                ),
                    
                        title = 'Simvisage.Spirrid',
                        resizable = True,
                        width = 1.0,
                        height = 1.0,
                        dock = 'tab',
                        buttons = [OKButton],
                        id = 'simvisage.spirrid_mine.dock') 

if __name__ == '__main__':


    global spirrid
    import cProfile

    spirrid = SPIRRID()
#    spirrid.response_function_type + 'FilEnergyDebonding'
#    spirrid._reset_response_function()
    
    if True:
        #cProfile.run('spirrid.displ_force_stdev', 'spirrid_tprof' )
        cProfile.run('spirrid.configure_traits()', 'spirrid_tprof' )
        
        import pstats
        p = pstats.Stats('spirrid_tprof')
        p.strip_dirs()
        print 'cumulative'
        p.sort_stats('cumulative').print_stats(50)
        print 'time'
        p.sort_stats('time').print_stats(50)
        #spirrid.configure_traits()

    else:
        spirrid.configure_traits()
    