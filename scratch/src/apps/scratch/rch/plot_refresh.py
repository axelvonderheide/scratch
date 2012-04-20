
from enthought.traits.api import \
    Float, HasTraits, Int, List, Array, Dict, String, WeakRef, Instance, Enum, Property, Any, Delegate, \
    cached_property, Constant, Button, on_trait_change, Bool

from stats.pdistrib.pdistrib import IPDistrib, PDistrib
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from enthought.traits.ui.api \
    import TabularEditor, View, Item, Group, HGroup, VGroup, HSplit, VSplit 
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
    array, mgrid, ogrid, sqrt, max, pi,e, argmax, argmin
from scipy.special import erf

#----------------------------------------------------------------------------------
#                                     SPIRRID                                     
#----------------------------------------------------------------------------------
class SPIRRID( HasTraits ):
    """
    Class for evaluating the Phoenix integral representing the strain-based
    fiber-bundle model.
    """
    plot_size_effect = Instance(Plot)
    def _plot_size_effect_default(self):
        p = Plot()
        p.tools.append(PanTool(p))
        p.overlays.append(ZoomTool(p))
        return p
        
    eval_size_effect = Button( 'Evaluate size effect curve' )
    def _eval_size_effect_fired(self):
        p = self.plot_size_effect
        N_array = linspace( 0, 100, 20 )
        means_array = N_array ** 2
    
        for name in p.plots.keys():
            p.delplot( name )
            
        p.datasources[ 'strength' ] = ArrayDataSource( means_array,
                                                        sort_order = 'none' )
        p.datasources[ 'cracks' ] = ArrayDataSource( N_array,
                                                        sort_order = 'none' )
#        p.datasources[ 'strngth' ] = ArrayDataSource( means_array + 1e-5,
#                                                        sort_order = 'none' )
        p.plot( ('cracks','strength'), name = 'size_effect1', color = 'red' )
        p.request_redraw()
#        p.plot( ('cracks','strngth'), name = 'size_effect2', color = 'red' )

    # show the user interface
    traits_view = View(                       Item('eval_size_effect', show_label = False),
                                              Item('plot_size_effect',                                            
                                                     editor=ComponentEditor(), 
                                                     show_label = False,
                                                     resizable = True
                                                     ),
                        title = 'Simvisage.Spirrid',
                        resizable = True,
                        width = 1.0,
                        height = 1.0,
                        dock = 'tab',
                        id = 'simvisage.spirrid.dock'
                        )

if __name__ == '__main__':


    spirrid = SPIRRID()
    spirrid.configure_traits()
    