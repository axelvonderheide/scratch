'''
Created on Apr 30, 2009

@author: jakub
'''
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property


class FERefinementLevelGrid(HasTraits):
    
    dof_offset = Int(0)
    
    s_list = List
    
    subgrids = Property(List)
#    def _set_subgrids(list):
#        for i,grid in enumerate(list):
#            grid.order = i
#        return list
#    
    def _get_subgrids(self):
        for i,grid in enumerate(self.s_list):
            grid.order = i
            if i == 0:
                grid.prev_grid = None
            else:
                grid.prev_grid = self.s_list[i-1]
        return self.s_list
    
    n_dofs = Property()
    #no need to be cached? summing small number of integers
    def _get_n_dofs(self):
        n_dofs = 0
        for grid in self.subgrids:
            n_dofs += grid.n_dofs
        return n_dofs
    
class IFEGrid(HasTraits):
    pass    
    
class FEGrid(IFEGrid):
    
    order = Int(0)
    prev_grid = Instance(IFEGrid)
    a = Int(1)
    
    dof_offset = Property(Int, depends_on = 'prev_grid.n_dofs,prev_grid.dof_offset')
    @cached_property
    def _get_dof_offset(self):
        if self.prev_grid:
            return self.prev_grid.n_dofs+self.prev_grid.dof_offset
        else:
            return 0

    n_dofs = Property(Int,depends_on = 'a' )
    @cached_property
    def _get_n_dofs(self):
        return self.a
        
if __name__ == '__main__':
    g1 = FEGrid()
    g2 = FEGrid()
    g3 = FEGrid()
    
    i_grid =  FERefinementLevelGrid()
    i_grid.s_list = [g1,g2,g3]
    print "n_dofs", i_grid.n_dofs
    print "order ", g3.order
    print "offset ", g3.dof_offset
    
    g1.a = 5
    print "AFTER ",g1.n_dofs
    print "n_dofs", i_grid.n_dofs
    print "order ", g3.order
    print "offset ", g3.dof_offset
    