from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict

from numpy import array

class LILMatrix(HasTraits):
    '''
    class that mimics the sparse.lil_matrix from scipy 
    '''
    ndofs = Int()
    rows = Property(Array)
    @cached_property
    def _get_rows(self):
        return array(self.ndofs,1)
    
    data = Property(Array)
    @cached_property
    def _get_data(self):
        return array(self.ndofs,1)
    
    def set_item(self, ix, matrix):
        '''
        Sets the the submatrix to given indices
        @param ix:
        @param matrix:
        '''
        for i in ix:
            self.rows[i]= ix 
        self.data[ix] = matrix
        
if __name__ == "__main__":
    
    A = array([[1,2],[3,4]])
    ix = [3,4]
    mtx = LILMatrix(ndofs = 10)
    mtx.set_item(ix, A)
    mtx.configure_traits()
        