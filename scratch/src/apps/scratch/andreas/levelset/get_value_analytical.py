'''
Created on Sep 15, 2009

@author: Andreas
'''
from numpy import dstack, array, zeros_like, mgrid
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict, \
     Class
     
from enthought.mayavi.mlab import surf, colorbar, show, mesh
from time import time

class LevelSet(HasTraits):
    
    out_grid = Array
    scalars = Array
    
    def get_value(self, grid):
        t_start = time()
        # f1=y    f2=y+1    fg=x
        x=dstack(grid)[0][0]
        y=dstack(grid)[0][1]
        #x, y = grid
        def f1(x,y):
            return y
        def f2(x,y):
            return y+1
        def fg(x,y,xg):
            return x-xg
        result=zeros_like(x)
        fgpos=(fg(x,y,0)>=0)    
        for i in range(0,len(grid)):
            if fgpos[i]==True:
                result[i]=f1(x[i],y[i])
            elif fgpos[i]==False:
                result[i]=f2(x[i],y[i])
        self.out_grid = grid
        self.scalars = result 
        print "Elapsed time     : %8.2f sec" % (time () - t_start)
        #return grid, result
        
    @show
    def gsurf(self): 
        """Test contour_surf on regularly spaced co-ordinates like MayaVi."""
        x, y = self.out_grid.T
        s = surf(x, y, self.scalars, warp_scale = 0.)
        colorbar(s)
        return s
#test   
if __name__ == '__main__':
    #print 'distance ',get_value(array([[2,2],[3,-3],[-2,-2],[-2,1]]))
    ls = LevelSet()
    print 'distance ',ls.get_value(array([[2,2],[3,-3],[-2,-2],[-2,1]]))
    #print 'distance ',ls.get_value(mgrid[-10.:10.:0.25, -10.:10.0:0.25])
    print ls.out_grid
    ls.gsurf()