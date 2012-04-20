
from enthought.traits.api import Array, HasTraits

from ibvp_solve.mfn_line_editor import MFnWTPlotItem

from enthought.traits.ui.api import Item, View, Group

from enthought.chaco.api import \
     Plot, AbstractPlotData, ArrayPlotData

from numpy import array

from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
                                     MenuBar, Separator

from math import cos



class PhiDiagram( HasTraits ):

    xdata = Array
    ydata = Array
    
    
    
    traits_view = View( Group( Item("plot_type"),
                               orientation = 'horizontal' ),
                        MFnWTPlotItem("xdata", "ydata", 
                                      type_trait="plot_type",
                                      type="scatter",
                                      
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
                        menubar=MenuBar(Menu(Action(name="Print data",
                                                    action="print_data"),
                                             name="View")),
                        buttons= [OKButton,CancelButton],
                        resizable=True,
                        width=500, height=800)

    
    def __call__(self, x):
        return cos(x)
    
if __name__ == '__main__':    

    #phi_fn.configure_traits( view = 'xtraits_view' )

    #mf = MFnLineArray( xdata = [ 0.0, Pi / 6. - 0.01, Pi / 6. + 0.01, Pi ],
     #                   ydata = [ 1.0, 1.0, 0.0, 0.0] )
    

    fx = PhiDiagram()
    v = fx(60)
    print 'result=', v
    fx.configure_traits()
    
    