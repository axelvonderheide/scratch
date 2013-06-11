from enthought.chaco import api as chaco
from numpy import arange, pi, sin, cos
from enthought.enable.example_support import DemoFrame, demo_main
from enthought.enable.api import Window
from enthought.traits.api import false
from enthought.chaco.api import create_polar_plot

numpoints = 5000
low = 0
high = 2*pi
theta = arange(low, high, (high-low) / numpoints)

# Create the radius data
radius = cos(3*theta)

# Create a new polar plot with radius and theta data
myplot = create_polar_plot((radius,theta),color=(0.0,0.0,1.0,1), width=4.0)

# We now need to set the plot's size, and add a little padding for the axes.
# (Normally, when Chaco plots are placed inside WX windows, the bounds are
# set automatically by the window.)
myplot.padding = 50
myplot.bounds = [400,400]
myplot.outer_bounds = [600,600]

def main():
    # Now we create a canvas of the appropriate size and ask it to render
    # our component.  (If we wanted to display this plot in a window, we
    # would not need to create the graphics context ourselves; it would be
    # created for us by the window.)
    plot_gc = chaco.PlotGraphicsContext(myplot.outer_bounds)
    plot_gc.render_component(myplot)
    
    # Get the directory to save the image in
    import os, sys
    print 'Please enter a path in which to place generated plots.'
    print 'Press <ENTER> to generate in the current directory.'
    path = raw_input('Path: ').strip()
    
    if len(path) > 0 and not os.path.exists(path):
        print 'The given path does not exist.'
        sys.exit()
    
    # The file name to save the plot as
    file_name = "simple_polar.png"
        
    if not os.path.isabs(path):
        print 'Creating image: ' + os.path.join(os.getcwd(), path, file_name)
    else:
        print 'Creating image: ' + os.path.join(path, file_name)
    
    # Finally, we tell the graphics context to save itself to disk as an image.
    plot_gc.save(os.path.join(path, file_name))

if __name__ == '__main__':
    main()
