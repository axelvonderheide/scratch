'''
Created on Aug 5, 2009

@author: jakub
'''

from enthought.mayavi import mlab

@mlab.show
def image():
   src = mlab.pipeline.open('eps_app1nodes_0.vtk')
   tc = mlab.pipeline.extract_tensor_components(src)
   mlab.pipeline.surface(tc)

if __name__ == '__main__':
    image()
