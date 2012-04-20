'''
Created on Apr 6, 2009

@author: jakub
'''
'''
Created on Mar 6, 2009

@author: jakub
'''

from enthought.tvtk.api import tvtk
from numpy import array

from enthought.mayavi.scripts import mayavi2
# Now view the data.
@mayavi2.standalone
def view():
    from enthought.mayavi.sources.api import VTKDataSource, VTKFileReader
    
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.modules.vectors import Vectors
    from enthought.mayavi.modules.api import TensorGlyph
    from enthought.mayavi.filters.api import WarpVector, ExtractTensorComponents

    mayavi.new_scene()
    # The single type one
    #src = VTKDataSource(data = '/home/jakub/simdata/spalling/sig_app9nodes_2.vtk')
    
    #src = VTKFileReader(base_file_name = '/home/jakub/simdata/spalling/nodes_0.vtk')
    src = VTKFileReader()
    VTKFileReader()
    mayavi.add_source(src) 
    #mayavi.open_vtk('/home/jakub/simdata/spalling/sig_app9nodes_2.vtk')
    #tg = TensorGlyph()
    #tg.glyph.glyph.scale_factor = 100.
    #mayavi.add_module(tg)
    #mayavi.add_module(Surface())
    #mayavi.add_module(Vectors())
    #wv = WarpVector()
    #tc = ExtractTensorComponents()
    #mayavi.add_filter(wv)
    #mayavi.add_filter(tc)
    mayavi.add_module(Surface())
    src.base_file_name = '/home/jakub/simdata/spalling/eps_app1nodes_0.vtk'

if __name__ == '__main__':
    view()
