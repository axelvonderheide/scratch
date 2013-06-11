
from enthought.traits.api import \
    Array, Bool, Callable, Enum, File, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
    Dict, Property, cached_property, WeakRef, Delegate, \
    ToolbarButton, on_trait_change

from enthought.traits.ui.api import \
    Item, View, HGroup, ListEditor, VGroup, VSplit, Group, HSplit

from enthought.traits.ui.menu import \
    NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
    MenuBar, Separator

from numpy import float_, zeros, arange, array, copy

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

#from rt_domain import RTraceDomainField

class RTrace( HasTraits ):
    name = Str( 'unnamed' )
    update_on = Enum( 'update', 'iteration' )
    clear_on = Enum( 'never', 'update' )
    save_on = Enum( None )

    #sctx = WeakRef( SContext )
    rmgr = WeakRef()

    # path to directory to store the data
    dir = Property
    def _get_dir( self ):
        return self.rmgr.dir

    # path to the file to store the data
    file = File

    def setup( self ):
        # In this case, the context must be a single point
        #
        pass

    refresh_button = ToolbarButton( 'Refresh',
                                   style = 'toolbar' )

    @on_trait_change( 'refresh_button' )
    def refresh( self, event = None ):
        self.redraw()

    def add_current_displ( self, sctx, U_k ):#TODO: to avoid class checking in rmngr - UGLY
        pass

    def register_mv_pipelines( self, e ):
        '''
        Eventually register pipeline components within the mayavi sceen.
        
        do nothing by default
        '''

