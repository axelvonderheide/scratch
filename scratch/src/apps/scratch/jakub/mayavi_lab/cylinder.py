'''
Created on Apr 16, 2009

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
    from enthought.mayavi.sources.api import VTKDataSource, ParametricSurface
    from enthought.mayavi.modules.outline import Outline
    from enthought.mayavi.modules.surface import Surface
    from enthought.mayavi.modules.vectors import Vectors
    from enthought.mayavi.modules.api import TensorGlyph
    from enthought.mayavi.filters.api import WarpVector, ExtractTensorComponents

    mayavi.new_scene()
    mesh = tvtk.Cylinder()
    # The single type one
    src = ParametricSurface(#function = 'boy'
                            parametric_function = mesh)
    #src = VTKDataSource(data = mesh)
    mayavi.add_source(src)
    mayavi.add_module(Surface())
    
 

if __name__ == '__main__':
    view()
