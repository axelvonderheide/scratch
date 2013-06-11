


#######################################################################################3
from enthought.traits.api import Array, Bool, Callable, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, ToolbarButton, Button, on_trait_change

from enthought.traits.ui.api import Item, View, Group, Handler


# Chaco imports
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
from enthought.enable.component_editor import ComponentEditor

from enthought.mayavi.core.source import Source

class MFnLineArray(HasTraits):

    # Public Traits
    xdata = Array
    ydata = Array

    traits_view = View( ChacoPlotItem("xdata", "ydata", 
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
                        resizable=True,
                        width=500, height=800)

if __name__ == '__main__':
    rm1 = MFnLineArray()
    rm1.configure_traits()
