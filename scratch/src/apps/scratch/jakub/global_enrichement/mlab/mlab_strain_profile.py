'''
Created on Apr 30, 2010

@author: jakub
'''
from enthought.mayavi import mlab
from os import chdir, getcwd

@mlab.show
def image():
    #print getcwd()
    fig = mlab.figure()
    fig.scene.background = (1.,1.,1.)
    fig.scene.jpeg_quality = 100
    
    src = mlab.pipeline.open('/home/jakub/simdata/global_enrichement/eps_f0nodes_1.vtk')
    tc = mlab.pipeline.extract_tensor_components(src)
    tc.filter.scalar_mode = 'component' 
    ws = mlab.pipeline.warp_scalar(tc, warp_scale = 2.)
    #mlab.pipeline.surface(ws)
    cp = mlab.pipeline.cut_plane(ws)
    
    cp.filters[0].widget.enabled = False
    cp.filters[0].plane.normal = (0.,1.,0.)
    su = mlab.pipeline.surface(cp)
    su.actor.property.color = (0.75294117647058822, 0.75294117647058822, 0.75294117647058822)
    su.actor.property.line_width = 4.
    su.actor.mapper.scalar_visibility = False
    
    #print cp.filters
    src2 = mlab.pipeline.open('/home/jakub/simdata/global_enrichement/eps_m0nodes_1.vtk')
    tc2 = mlab.pipeline.extract_tensor_components(src2)
    tc2.filter.scalar_mode = 'component' 
    ws2 = mlab.pipeline.warp_scalar(tc2, warp_scale = 2.)
    #mlab.pipeline.surface(ws)
    cp2 = mlab.pipeline.cut_plane(ws2)
    
    cp2.filters[0].widget.enabled = False
    cp2.filters[0].plane.normal = (0.,1.,0.)
    
    su2 = mlab.pipeline.surface(cp2)
    su2.actor.property.color = (0.,0.,0.)
    su2.actor.property.line_width = 4.
    su2.actor.mapper.scalar_visibility = False
    ax = mlab.axes(color = (0.,0.,0.),
                   xlabel = 'x',
                   ylabel = '',
                   zlabel = 'strain',
                   extent = [0., 3.1, 0., 0., 0., 2.],
                   ranges = [1., 0., 0., 0., 0., 1.],
                   y_axis_visibility = False)
    ax.title_text_property.color = (0.0, 0.0, 0.0)
    ax.label_text_property.color = (0.0, 0.0, 0.0)
    ax.title_text_property.bold = False
    ax.label_text_property.bold = False
    mlab.view(90., 90., 6.0, 'auto')
    
    #mlab.roll(125)
    
    #mlab.savefig('quad5.png')
    
if __name__ == '__main__':
    image()