
from threading import Thread
from time import sleep
from enthought.traits.api import *
from enthought.traits.ui.api import View, Item, Group, VSplit, HSplit, \
     Handler

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, \
     Action
from mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
import math
import numpy
import wx

class TCHandler( Handler ):
    def setattr( self, info, object, name, value ):
        print 'setting attribute', name
        Handler.setattr( self, info, object, name, value )
        info.object._updated += 1
    def object__updated_changed( self, info ):
        print 'updated - changed'
        if info.initialized:
            info.ui.title += '*'
    def recalculate( self, info ):
        info.object.draw()

class MFnRange( HasTraits ):
    min = Float(0.0)
    max = Float(1.0)
    n = Int(100)
    _updated = Int

    view = View(
        Group(
        Item('min', label = 'minimum', resizable=True,
             tooltip = 'Lower bound of the x-axes'),
        Item('max', label = 'maximum', resizable=True,
             tooltip = 'Upper bound of the x-axes'),
        Item('n', label = 'number of points', resizable=True,
             tooltip = 'Number of points used for plotting'),
        orientation='horizontal'
        )
        #handler = TCHandler()
        )

class MFnSinFn( HasTraits ):
    a = Float(1.0, label='first parameter')
    b = Float(1.0, label='second parameter')

    def _a_changed( self, name ):
        print 'a changed'
    
    def _b_changed( self, name ):
        print 'b changed'
    
    def eval(self, x):
        return self.a * x * math.sin( self.b* x );

RecalcAction = Action( name = 'Recalculate', action = 'recalculate' )

class MFnEditor( HasTraits ):
    x_range = Instance( MFnRange )
    def _x_range_default(self):
        return MFnRange()
    
    figure = Instance(Figure)
    function = Instance( MFnSinFn )
    def _function_default(self):
        return MFnSinFn()
    
    _updated = Int

    def draw(self):
        dx = math.fabs( self.x_range.min-self.x_range.max) / self.x_range.n
        X = numpy.arange( self.x_range.min, self.x_range.max+dx, dx  )
        ufunc_get_val = numpy.frompyfunc( self.function.eval, 1, 1 )
        Y = numpy.array( ufunc_get_val( X ), 'float_' )
        #self.figure.axes[0].clear()
        self.figure.axes[0].plot( X, Y )
        self.figure.axes[0].set_xscale('log' ) #, subsx = [0, 1, 2, 3 ] )
        self.figure.axes[0].set_yscale('log' ) # , subsy = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] )

    def _figure_default(self):
        figure = Figure()
        figure.add_axes([0.05, 0.04, 0.9, 0.92])
        return figure

    def __call__( self, x ):
        self.eval(x)

    view = View(
        HSplit( Item('function', style="custom"),
                Group( Item('figure',  editor=MPLFigureEditor(),
                            dock='horizontal', height=400, width=600,
                            resizable=True),
                       Item('x_range', style="custom"),
                       show_labels=False,
                       help = 'This is the text associated with a view'
                       ),
                show_labels=False
                ),
        handler = TCHandler(),
        resizable=True,
        height=0.75, width=0.75,
        buttons= [OKButton,CancelButton,RecalcAction])


mfn = MFnEditor()
mfn.draw()
mfn.configure_traits(filename='my_file.dac')
