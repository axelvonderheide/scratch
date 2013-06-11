from math import e, sin
from enthought.traits.api import \
    HasTraits, Float, Int, List, Array, Interface, \
    Str, Tuple, Property, cached_property, \
    Any, Enum, implements, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, HSplit, VGroup, Group

from numpy import frompyfunc, linspace, array, log as ln
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

class IntervalDividing( HasTraits) :
    """ Enter the interval <a,b> and maximum absolute error for terminating the process. """
    """ The root of the function f(x) must be defined at this interval. """
    a = Float(0, auto_set = False, enter_set = True)
    b = Float(2, auto_set = False, enter_set = True)
    error = Float(0.0001, auto_set = False, enter_set = True)
#    function = Expression(e**x -2*x - x**x*sin(e**x), auto_set = False, enter_set = True)

    
    def __init__(self,**kw):
        super(IntervalDividing,self).__init__(**kw)
        self._replot()
    
    @on_trait_change('a,b,error,f')
    def _replot(self):
        xdata, ydata = self.get_plot_data()
        self.plot_data.xdata = xdata
        self.plot_data.ydata = ydata
   
    """ define the function f(x) """
         
    def f(self,x):
        y = e**x -2*x - x**x*sin(e**x)
        return y   

    def get_plot_data(self):
        xdata = linspace(self.a,self.b,500)
        farray = frompyfunc( self.f, 1, 1 )
        ydata = array( farray( xdata ), dtype=float )
        return xdata, ydata
 
    results = Property( Float, depends_on = 'f, a, b, error' )        
    @cached_property
    def _get_results(self):
        """dividing the interval <a,b>,
        returns aprox x, error estimation, No. of steps ..."""
        int = [self.a, self.b]
        if self.f(self.a) * self.f(self.b) > 0:
            print "None or more than 1 roots in selected interval"
        else:
            while abs(int[0] - int[1])/2 > self.error:
                if self.f(int[0]) * self.f((int[0] + int[1])*0.5) < 0:
                    int.insert(1,(int[0]+int[1])*0.5), int.pop()
                else:
                    int.insert(1,(int[0]+int[1])*0.5), int.pop(0) 
            return [(int[0] + int[1])/2,
                    abs(int[0] - int[1])/2,
                    (ln(2)-ln( (abs( int[0] - int[1] ) /2)/(self.b-self.a)))/ln(2)]


    value = Property( Float, depends_on = 'f, a, b, error' )
    @cached_property
    def _get_value(self):
        return self.results[0]

    err = Property( Float, depends_on = 'f, a, b, error' )
    @cached_property
    def _get_err(self):
        return self.results[1]

    steps = Property( Float, depends_on = 'f, a, b, error' )
    @cached_property
    def _get_steps(self):
        return self.results[2]


    
    plot_data = Instance(MFnLineArray)
    def _plot_data_default(self):
        return MFnLineArray()
    
    trait_view = View(HSplit(VGroup(Group(
                                   Item('a'),
                                   Item('b'),
                                   Item('error'),
                                   Item('value', style='readonly', label = 'approx value of x'),
                                   Item('err', style='readonly', label = 'error estimation'),
                                   Item('steps', style='readonly', label = 'No. of steps'),
                                   label = 'select interval <A,B> and max error',
                                   ),
                            ),
                            Item('plot_data',
                            style = 'custom',
                            show_label = False,
                            resizable = True,
                            )
                        ),                  
                       resizable = True,
                       width = 0.6,
                       height = 0.6,
                       )
                       
    

    
         
if __name__ == '__main__':
    idev = IntervalDividing()
    idev.configure_traits()