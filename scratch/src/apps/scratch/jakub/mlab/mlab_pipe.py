'''
Created on Aug 15, 2009

@author: jakub
'''

from enthought.traits.api import \
    HasTraits, Bool, Str, Enum, Int
    
from enthought.mayavi import mlab
from enthought.mayavi.sources.api import VTKFileReader, VTKDataSource

class MVUnstructuredGrid( HasTraits ):
 
    warp = Bool(False)
    idx = Int()
        
    name = Str
    position = Enum('nodes','int_pnts')

    @mlab.show  
    def rebuild_pipeline(self, pd, max_idx):
        self.src = VTKDataSource( name = self.name, data = pd )
        self.msrc = mlab.pipeline.add_dataset(self.src)
        #self.src = VTKFileReader(base_file_name = self.position+'_0.vtk')
        self.max_idx = max_idx
        self.idx = max_idx
        #src = mlab.pipeline.open('u1nodes_0.vtk')
        if self.warp:
            self.src.point_vectors_name = 'u_'+str(self.idx)
            self.msrc = mlab.pipeline.warp_vector(self.msrc)
    
        name = self.name +'_'+ str(self.idx)
        if name in self.src._point_tensors_list:
            self.src.point_tensors_name = name
            self.msrc = mlab.pipeline.extract_tensor_components(self.msrc)
        elif name in self.src._point_vectors_list:
            self.src.point_vectors_name = name
            self.msrc = mlab.pipeline.extract_vector_components(self.msrc)
        elif name in self.src._point_scalars_list:
            self.src.point_scalars_name = name

        surf = mlab.pipeline.surface(self.msrc)   
        surf.actor.property.point_size = 5

    def _update_fields(self):
        if self.warp:
            self.src.point_vectors_name = 'u_'+str(self.idx)
        name = self.name +'_'+ str(self.idx)
        if name in self.src._point_tensors_list:
            self.src.point_tensors_name = name
        elif name in self.src._point_vectors_list:
            self.src.point_vectors_name = name
        elif name in self.src._point_scalars_list:
            self.src.point_scalars_name = name

    def prev_ts(self):
        if self.idx > 0:
            self.idx -= 1  
            self._update_fields()   
    
    def next_ts(self):
        if self.idx < self.max_idx:
            self.idx += 1
            self._update_fields()   
