# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2010, Enthought
# License: BSD style

from enthought.mayavi import mlab
import numpy as np

mlab.clf()

# Number of lines
n_lines = 200#200
# Number of points per line
n_points = 100#100

# Create Example Coordinates
xyz = []
for ii in xrange( n_lines ):
    xyz.append( 
            np.cumsum( np.random.random( ( n_points, 3 ) ), axis=0 )
        )
xyz = np.array( xyz )
x, y, z = xyz.T
pts = mlab.pipeline.scalar_scatter( x.T, y.T, z.T )

connections = []
for i in range( n_lines ):
    connections.append( np.c_[np.arange( i * n_points, ( i + 1 ) * n_points - 1 ),
                                np.arange( i * n_points + 1, ( i + 1 ) * n_points )] )
pts.mlab_source.dataset.lines = np.vstack( connections )
# First option: tubes
if hasattr( mlab.pipeline, 'stripper' ):
    lines = mlab.pipeline.surface( 
                mlab.pipeline.tube( 
                    mlab.pipeline.stripper( 
                        pts,
                    ),
                    tube_sides=7, tube_radius=0.1,
                ),
            )
else:
    lines = mlab.pipeline.surface( 
                mlab.pipeline.tube( 
                    pts,

                    ),
                    tube_sides=7, tube_radius=0.1,
                )
mlab.show()





from enthought.mayavi import mlab
import numpy

# Create Example Coordinates
xyz = []
for ii in xrange( 50 ):
    xyz.append( numpy.array( [[0, 0, 0]] ) )
    for jj in xrange( 100 ):
        translation = xyz[ii][-1] + numpy.random.random( ( 1, 3 ) )
        xyz[ii] = numpy.concatenate( ( xyz[ii], translation ), axis=0 )

# Display (This is the part I need to optimize)
for ii in xrange( 50 ):
    mlab.plot3d( xyz[ii][:, 0], xyz[ii][:, 1], xyz[ii][:, 2], tube_sides=7, tube_radius=0.1 )
mlab.show()

import numpy as np

# The number of points per line
N = 300

# The scalar parameter for each line
t = np.linspace( -2 * np.pi, 2 * np.pi, N )

from enthought.mayavi import mlab
mlab.figure( 1, size=( 400, 400 ), bgcolor=( 0, 0, 0 ) )
mlab.clf()

# We create a list of positions and connections, each describing a line.
# We will collapse them in one array before plotting.
x = list()
y = list()
z = list()
s = list()
connections = list()

# The index of the current point in the total amount of points 
index = 0

# Create each line one after the other in a loop
for i in range( 50 ):
    x.append( np.sin( t ) )
    y.append( np.cos( ( 2 + .02 * i ) * t ) )
    z.append( np.cos( ( 3 + .02 * i ) * t ) )
    s.append( t )
    # This is the tricky part: in a line, each point is connected
    # to the one following it. We have to express this with the indices
    # of the final set of points once all lines have been combined
    # together, this is why we need to keep track of the total number of 
    # points already created (index)
    connections.append( np.vstack( 
                       [np.arange( index, index + N - 1.5 ),
                        np.arange( index + 1, index + N - .5 )]
                            ).T )    
    index += N

# Now collapse all positions, scalars and connections in big arrays
x = np.hstack( x )
y = np.hstack( y )
z = np.hstack( z )
s = np.hstack( s )
connections = np.vstack( connections )

# Create the points
src = mlab.pipeline.scalar_scatter( x, y, z, s )

# Connect them
src.mlab_source.dataset.lines = connections

# The stripper filter cleans up connected lines
lines = mlab.pipeline.stripper( src )

# Finally, display the set of lines
mlab.pipeline.surface( mlab.pipeline.tube( lines, tube_sides=7, tube_radius=0.1 ), opacity=.4, colormap='Accent' )

# And choose a nice view
mlab.view( 33.6, 106, 5.5, [0, 0, .05] )
mlab.roll( 125 )
mlab.show()



from enthought.mayavi import mlab
import numpy as np

mlab.clf()

# Number of lines
n_lines = 200
# Number of points per line
n_points = 100

# Create Example Coordinates
xyz = []
for ii in xrange( n_lines ):
    xyz.append( 
            np.cumsum( np.random.random( ( n_points, 3 ) ), axis=0 )
        )

fig = mlab.gcf()
fig.scene.disable_render = True
for this_xyz in xyz:
    x, y, z = this_xyz
    mlab.plot3d( x,
                y,
                z, tube_sides=7, tube_radius=0.1 )

fig.scene.disable_render = False

mlab.show()
