from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

from numpy import array, sort, where


def signum(x):
    return x/abs(x)

class LevelSet2D(HasTraits): 
    
    def get_2d_intersections(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError
    
    def get_value(X_pnt):
        raise NotImplementedError
    
class Plane(LevelSet2D):
    
    A = Float(0.)
    B = Float(0.)
    C = Float(0.)
    
    def get_2d_intersections(self, xmin, xmax, ymin, ymax):
        if self.A == 0. and ymin< (signum(self.B)*(-self.C))< ymax:
            return [xmin,-self.C],[xmax,-self.C]
        elif self.B == 0. and xmin< (signum(self.A)*(-self.C))< xmax:
            return [-self.C,ymin],[-self.C,ymax]
        elif self.A == 0. and self.B == 0.:
            print "A and B are zero, this is forbidden!"
        else:
            xval = array([(-self.B * ymin -self.C)/self.A,\
                         (-self.B * ymax -self.C)/self.A])
            xm = where(xmin < xval,xval,xmin)
            xcoord = where(xm < xmax,xm,xmax)
            
            yval = array([(-self.A * xmin -self.C)/self.B,\
                          (-self.A * xmax -self.C)/self.B])
            ym = where(ymin<yval, yval,ymin)
            ycoord = where(ym<ymax, ym,ymax)
            if (self.A * self.B) > 0:
                return [xcoord[0],ycoord[1]],[xcoord[1],ycoord[0]]
            elif (self.A * self.B) < 0:
                return [xcoord[0],ycoord[0]],[xcoord[1],ycoord[1]]            
        
    def get_value(self,X_pnt):
        return self.A * X_pnt[0]+self.B * X_pnt[1] +self.C
    
    
if __name__ == '__main__':
    ls = Plane(A = 0., B= 1., C = - 0.5 )
    
    #ls.configure_traits()
    X_pnt = [1.,0.,0.]
    print "val ",ls.get_value(X_pnt)
    print "intersect ",ls.get_2d_intersections(0., 1., 0., 1.)