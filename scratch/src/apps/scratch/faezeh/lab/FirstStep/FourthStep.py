
# Imports:
from enthought.traits.api import HasTraits, Array, Float
    
from enthought.traits.ui.api import View, Item, Group, Handler
    
from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

from mfn_line_editor import MFnWTPlotItem

from numpy import linspace

from math import cos

import sys

class MFnWTHandler( Handler ):
    def open_data( self, info ):
        info.object.print_data()

    def save_file(self,info):
        sys.exit(0)

    def exit_file(self,info):
        sys.exit(0)

    
    def init_file(self, window):
        """ Creates a new action. """
        self._window = window
        self.name = "E&xit"

    def perform(self):
        """ Performs the action. """
        self._window.close()
        

class X2COSX(HasTraits):
    ''' Simple class emulating a function x**2 * cos(x)
    '''

    a = Float(1.)
    b = Float(2.)
    c = Float(3.)

    def _a_changed(self, old, new):
        self.refresh()
    def _b_changed(self, old, new):
        self.refresh()
    def _c_changed(self, old, new):
        self.refresh()

    X = Array
    X2cos = Array

    def refresh(self):
        self.X = linspace(0.,50.,100)
        self.X2cos = [ self(x) for x in self.X ]

    def __call__(self, x):
        return self.a * x**2 * cos(self.b*x + self.c)
    
    traits_view = View( Group( Item( name = 'a'), 
                               Item( name = 'b'),
                               Item( name = 'c' ),
                               orientation = 'horizontal' ),
                        MFnWTPlotItem("X", "X2cos", 
                                      #type_trait="line",
                                      type="line",
                                      
                                      # Basic axis and label properties
                                      show_label=False,
                                      resizable=True,
                                      orientation="h",
                                      x_label = "X",
                                      y_label = "Y",
                                      title = '',
                                      
                                      # Plot properties
                                      color = "blue",
                                      bgcolor = "white",
                                      
                                      # Specific to scatter plot
                                      marker = "circle",
                                      marker_size = 5,
                                      outline_color = "none",
                                      
                                      # Border, padding properties
                                      border_visible=True,
                                      border_width=0,
                                      padding_bg_color = "white"),
                        buttons= [OKButton,CancelButton],
                        menubar=MenuBar(Menu(Action(name="O&pen..", action="open_data"),
                                             Action(name="S&ave..", action="save_file"),
                                             Action(name="E&xit", action="exit_file"),
                                             name = 'File')),
                        handler = MFnWTHandler,                                             
                        resizable=True,
                        width=500, height=800)
    
import pickle
    
if __name__ == '__main__':
    fx = X2COSX()
    print fx(3.0)
    
    fx.a = 5.0
    
    filename = "myfile-fx.dat"
    file = open(filename,'w')
    pickle.dump(fx, file)
    file.close()
    fx = None
    
    file = open(filename,'r')
    fx2 = pickle.load(file)
    file.close()
    fx2.refresh()
    fx2.configure_traits()
    
