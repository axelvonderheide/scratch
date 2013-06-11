'''
Created on Apr 28, 2009

@author: jakub
'''
import numpy as np
a = np.random.random((4, 4))
from enthought.mayavi.api import Engine
e = Engine()
e.start()
s = e.new_scene()
from enthought.mayavi.sources.api import ArraySource
src = ArraySource(scalar_data=a)
e.add_source(src)
from enthought.mayavi.filters.api import WarpScalar, PolyDataNormals
warp = WarpScalar()
e.add_filter(warp, obj=src)
normals = PolyDataNormals()
e.add_filter(normals, obj=warp)
from enthought.mayavi.modules.api import Surface
surf = Surface()
e.add_module(surf, obj=normals)