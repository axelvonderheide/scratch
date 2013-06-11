'''
Created on Sep 15, 2009

@author: Andreas
'''

#===============================================================================
# 
#    This class may be used to estimate a level set using signed distance function. The following functions are used within the class.
# 
#    set_grid()    set_fns()     _eval_fn()     get_value()    show()
#    
#    Functions and limits have to be implemented manually using python syntax. will therefore be generated with the set_fns() command.
#    some examples for syntax:
#  
#        >=    greater or equal to    
#        ==    equal to     
#        <=    smaller or equal to     
#        !=    unequal to   (this may be helpful for modeling cones(SDF for circles))
#        |     for or               
#        &     for and
#        ;     for separating function command line from limits command
# 
#    The limits have to include the hole domain, so that a function is defined everywhere, for each grid point, 
#    if not the level set of this point will be zero. 
#        
#    some examples: 
#    
#    linear and constant function         
#    
#    LevelSet.set_fns(['y+x;(x>0)','y;(x<0)']
#    
#    two circle
#    
#    a='((x**2+y**2)/abs(x**2+y**2)**0.5)-5.0;((x**2+y**2)**0.5<=7.5) & ((x**2+y**2)>0.0)'
#    b='-5;((x==0)&(y==0))'
#    c='(-(x**2+y**2)/abs(x**2+y**2)**0.5)+10;(x**2+y**2)**0.5>=7.5'
#    ls.set_fns([a,b,c])                        
#    
#    The grid may be imported or generated using set_grid.
#    Display of LevelSet is only available for regularly spaced grid. 
#
#===============================================================================

# Import of functions
from numpy import dstack, array, zeros_like, mgrid, where
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, Dict, \
     Class  
from enthought.mayavi.mlab import surf, colorbar, show, mesh
from time import time

class LevelSet(HasTraits):
    
    out_grid = Array                #grid 
    scalars = Array                 #scalar distance to center line
    x_elements = Int(100.0)         #number of elements in x-direction used by set_grid
    y_elements = Int(100.0)         #number of elements in y-direction used by set_grid
    x_limit = array([-10.0,10.0])   #limits in x-direction used by set_grid
    y_limit = array([-10.0,10.0])   #limits in y-direction used by set_grid
    fns = List                      #definition of function
    lim= List                       #definition of limits

    def set_grid(self):
        """ defines a grid within a certain area with x_e elements in x and y_e elements in y direction. 
            values x_elements, y_elements, x_limits, y_limits should be defined in advance before starting set_grid"""        
        x_e=self.x_elements
        y_e=self.y_elements
        x_l=self.x_limit
        y_l=self.y_limit      
        a=(x_l[1]-x_l[0])/x_e
        b=(y_l[1]-y_l[0])/y_e
        return mgrid[x_l[0]:x_l[1]:a, y_l[0]:y_l[1]:b]
  
  
    def set_fns(self, fn_list):
        """ set_fns evaluates the functions and limits from a string and will be evaluated by python. General python  syntax has to be used.
            the following form has to be followed: set_fns('function#1,limit#1';...;'function#n,limit#n')"""
        self.fns = []
        self.lim = []
        for fn in fn_list.split(';'):
            print fn
            afn, limits = fn.split('@')
            self.fns.append(afn)
            self.lim.append(limits) 
            

    def _eval_fn(self,x,y,fn, limits):#J:methods begining with underscore - private methods, not for user, helper methods used just internaly in the class
        """evaluates the value for each function, as defined in set_fns will be 0 if not within the defined domain"""
        return where(eval(limits),eval(fn),False)          


    def get_value(self, grid):
        """ get_value() evaluates the value for all function over, as defined in set_fns, 
            will be 0 if grid point is not defined in any function. Time for calculation is printed."""
        self.out_grid = grid
        t_start = time()
        x,y=grid
        self.scalars=zeros_like(x)
        for fn, limit in zip(self.fns,self.lim):
            self.scalars=where(self.scalars==0,self._eval_fn(x,y, fn, limit),self.scalars)
        print "Elapsed time:   %8.6f sec" % (time () - t_start) 

        
    @show
    def show(self): 
        """Test contour_surf on regularly spaced co-ordinates like MayaVi."""
        x, y = self.out_grid
        s = surf(x, y, self.scalars, warp_scale = 0.)
        return s
    
#test   
if __name__ == '__main__':
    ls = LevelSet()
    ls.x_elements=int(397.0)
    ls.y_elements=int(397.0)
    ls.x_limit=([-0.125,0.125])
    ls.y_limit=([-0.010,0.050])
    a='(y-0.040) @ ((x>=-0.0578)&(x<=0.0578)&(y>=0.020))'
    b='((y-0.029)**2+(x-0.0578)**2)**0.5-0.011 @ (x<=0.1400)&(x>=0.0578)&(y>=0.020)'
    c='((y-0.029)**2+(x+0.0578)**2)**0.5-0.011 @ (x>=-0.1400)&(x<=-0.0578)&(y>=0.020)'
    d='-(((y-0.011)**2+(x-(0.0578+0.01264911))**2)**0.5-0.011) @ (x>=0)&(x<=(0.0578+0.01264911))&(y<=0.020)'
    e='-(((y-0.011)**2+(x+0.0578+0.01264911)**2)**0.5-0.011) @ (x<=0)&(x>=-(0.0578+0.01264911))&(y<=0.020)'
    f='y @ ((x<=-(0.0578+0.01264911))|(x>=(0.0578+0.01264911)))&(y<=0.02)'
    ls.set_fns(a+";"+b+";"+c+";"+d+";"+e+";"+f)
    ls.get_value(ls.set_grid())
    ls.show()