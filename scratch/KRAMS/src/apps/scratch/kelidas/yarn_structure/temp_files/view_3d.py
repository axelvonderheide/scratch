from enthought.traits.api import HasTraits, Float, Property, \
    cached_property, Event, Array, Instance, File, Int, Directory, Button, Enum
from enthought.traits.ui.api \
    import View, Item, FileEditor, DirectoryEditor, HistoryEditor
from enthought.traits.ui.menu import OKButton, CancelButton, Action, Menu, \
    MenuBar 
from numpy import \
    loadtxt, ones_like, vstack, c_, hstack, array, cumsum, \
    zeros_like, zeros, min, arange

#import wxversion

#wxversion.select( '2.8' )

from os.path import join
import os, re

from enthought.tvtk.api import tvtk
from enthought.mayavi.scripts import mayavi2

from enthought.mayavi import mlab
from yarn_structure import Data

class View_3d( HasTraits ):
    show_button = Button()
    data = Data() # Instance( Data )
    
    def _show_button_fired( self ):
        p, r = self.data.slice
    
        
        
        mlab.figure( 'Continuous fibers' )
        mlab.plot3d( p[:, 0], p[:, 1], p[:, 2], r,
                     tube_radius=20, tube_sides=20, colormap='Spectral' ) 
        
#        n_cut = self.data.bq.shape[1] - 1
#        for i in range( 0, len( p ) / n_cut ):  
#            mlab.plot3d( p[( i * n_cut ):( ( i + 1 ) * n_cut ), 0], p[( i * n_cut ):( ( i + 1 ) * n_cut ), 1], p[( i * n_cut ):( ( i + 1 ) * n_cut ), 2], r[( i * n_cut ):( ( i + 1 ) * n_cut )],
#                         tube_radius=20, tube_sides=20, colormap='Spectral' )
         
        #mlab.title( 'Continuous fibers' )
  
@mayavi2.standalone     
def Scatter():
    # Create some random points to view.
    data = Data()
    p, r = data.slice
    pd = tvtk.PolyData()
    pd.points = p
    verts = arange( 0, len( r ), 1 )
    verts.shape = ( len( r ), 1 )
    pd.verts = verts
    pd.point_data.scalars = r
    pd.point_data.scalars.name = 'scalars'

    # Now visualize it using mayavi2.
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface

    mayavi.new_scene()
    d = VTKDataSource()
    d.data = pd
    mayavi.add_source( d )
    mayavi.add_module( Outline() )
    s = Surface()
    mayavi.add_module( s )
    s.actor.property.set( representation='p', point_size=2 )

@mayavi2.standalone     
def Lines():
    # Create some random points to view.
    data = Data()
    p, r = data.slice
    pd = tvtk.PolyData()
    pd.points = p
    verts = arange( 0, len( r ), 1 )
    verts.shape = ( len( r ), 1 )
    pd.verts = verts
    pd.point_data.scalars = r
    pd.point_data.scalars.name = 'scalars'
    pd.lines = arange( 0, len( r ), 1 ).reshape( len( r ) / 15, 15 )
    
    # Now visualize it using mayavi2.
    from enthought.mayavi.sources.vtk_data_source import VTKDataSource
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.modules.iso_surface import IsoSurface
    from enthought.mayavi.modules.streamline import Streamline

    mayavi.new_scene()
    d = VTKDataSource()
    d.data = pd
    mayavi.add_source( d )
    o = Outline()
    mayavi.add_module( o )

    s = Surface() #enable_contours=True
    mayavi.add_module( s )
    s.actor.property.set( representation='p', point_size=2 )
    
def Tubes():
    from enthought.mayavi import mlab
    import numpy as np
    
    mlab.clf()
    data = Data()
    p, r = data.slice
    
    # display in nodes different variables e.g. bq 
    r = data.bq.flatten()
    x, y, z = p.T
    pts = mlab.pipeline.scalar_scatter( x.T, y.T, z.T , r )
    
    connections = arange( 0, len( r ), 1 ).reshape( len( r ) / 15, 15 ) 
#    for i in range( n_lines ):
#        connections.append( np.c_[np.arange( i * n_points, ( i + 1 ) * n_points - 1 ),
#                                    np.arange( i * n_points + 1, ( i + 1 ) * n_points )] )
    pts.mlab_source.dataset.lines = connections
    # First option: tubes
    if hasattr( mlab.pipeline, 'stripper' ):
        print 'stripper'
        lines = mlab.pipeline.surface( 
                    mlab.pipeline.tube( 
                        mlab.pipeline.stripper( 
                            pts,
                        ),
                        tube_sides=10, tube_radius=15,
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

        
if __name__ == '__main__':
    #Scatter()
    Tubes()
    #Lines()
    #data = View_3d()
    #data._show_button_fired()
    
