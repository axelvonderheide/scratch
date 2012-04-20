"""
Example of how to use a DataView and bare renderers to create plots
"""

from numpy import linspace, sin, cos, array

from enthought.chaco.api import DataView, ArrayDataSource, ScatterPlot, LinePlot, LinearMapper, DataRange1D, FilledLinePlot
from enthought.chaco.tools.api import PanTool, ZoomTool
from enthought.enable.api import Window
from enthought.enable.example_support import DemoFrame, demo_main

class PlotFrame(DemoFrame):
    def _create_window(self):

        x = linspace(-5, 10, 500)
        
        y = 0.25 * cos(x)
        y2 = 0.5 * sin(2 * x)

        view = DataView(border_visible = True)

        x_ds = ArrayDataSource(x)
        y_ds = ArrayDataSource(y, sort_order="none")


        lineRed = FilledLinePlot(index = x_ds, value = y_ds,
                                    index_mapper = LinearMapper(range=view.index_range),
                                    value_mapper = LinearMapper(range=view.value_range),
                                    edge_color = "red",
                                    face_color = "paleturquoise",
                                    bgcolor = "white",
                                    border_visible = True)

        
#        lineRed = LinePlot(index = ArrayDataSource(x),
#                              value = ArrayDataSource(y),
#                              marker = "square",
#                              color = "red",
#                              outline_color = "transparent",
#                              index_mapper = LinearMapper(range=view.index_range),
#                              value_mapper = LinearMapper(range=view.value_range))



        lineBlue = LinePlot(index = lineRed.index,
                        value = ArrayDataSource(y2),
                        color = "blue",
                        index_mapper = LinearMapper(range=view.index_range),
                        value_mapper = LinearMapper(range=view.value_range))

        vertical_line = LinePlot(index = ArrayDataSource( array( [ 2.5, 2.5 ], dtype = float ) ),
                        value = ArrayDataSource( array( [ -0.4, 0.4 ], dtype = float ) ),
                        color = "green",
                        index_mapper = LinearMapper(range=view.index_range),
                        value_mapper = LinearMapper(range=view.value_range))


        # Add the plot's index and value datasources to the dataview's
        # ranges so that it can auto-scale and fit appropriately
        view.index_range.sources.append(vertical_line.index)
        view.value_range.sources.append(vertical_line.value)
        view.index_range.sources.append(lineRed.index)
        view.value_range.sources.append(lineRed.value)
        view.value_range.sources.append(lineBlue.value)

        # Add the renderers to the dataview.  The z-order is determined
        # by the order in which renderers are added.
        view.add(lineRed)
        view.add(lineBlue)
        view.add(vertical_line)

        view.tools.append(PanTool(view))
        view.overlays.append(ZoomTool(view))

        return Window(self, -1, component=view)

if __name__ == "__main__":
    demo_main(PlotFrame, size=(800,700), title="two graphs example")


