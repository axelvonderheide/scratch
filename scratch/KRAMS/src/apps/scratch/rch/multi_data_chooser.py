

from enthought.traits.api import \
    HasTraits, Instance, Button, Bool, on_trait_change, \
    Trait
from enthought.traits.ui.api import \
    View, Item
from enthought.chaco.api import \
    ArrayDataSource, Plot
from enthought.enable.component_editor import \
    ComponentEditor
from scipy.special import jn
from numpy import linspace

funcs = {'jn0': lambda x: jn(0,x), 
         'jn1': lambda x: jn(1,x),
         'jn2': lambda x: jn(2,x) }

class MultiDataChooser( HasTraits ):

    data_name = Trait( 'jn0', funcs )
    def _data_name_changed(self):
        # how to clear the plot?
        # I tried self.plot.delplot( name )
        # I tried self.plot.plots = {}
        # Even reconstruction of Plot does not show
        # any reaction
        self.plot = Plot()
        self.plot.datasources = self.ds
        self.plot.plot( ('x',self.data_name), color = 'red' )

    plot = Instance( Plot )

    def __init__(self,**kw):
        super(MultiDataChooser,self).__init__(**kw)
        x = linspace( 0,8,100 )
        
        # populate the plot with data sources
        self.ds = {}
        self.ds['x'] = ArrayDataSource( data = x )
        for key, fn in funcs.items():
            self.ds[key] = ArrayDataSource( data = fn( x ) )

        p = Plot()
        p.datasources = self.ds
        p.plot( ('x','jn0'), color = 'red' )
        self.plot = p

    traits_view = View( Item('data_name'),
                        Item('plot',show_label = False,
                             editor = ComponentEditor() ),
                        width = 0.8, 
                        height = 0.8, 
                        resizable = True )


MultiDataChooser().configure_traits()
