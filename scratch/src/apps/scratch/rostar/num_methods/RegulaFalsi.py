from math import e
from enthought.traits.api import \
    HasTraits, Float, Int, List, Array, Interface, \
    String, Tuple, Property, cached_property, \
    Any, Enum, implements, Instance, on_trait_change

from enthought.traits.ui.api import View, Item, HSplit, VGroup, Group

from numpy import frompyfunc, linspace, array, log as ln
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

class RegulaFalsi( HasTraits) :
    """ Enter the interval <a,b> and maximum absolute error for terminating the process. """
    """ The root of the function f(x) must be defined at this interval. """
    a = Float(0, auto_set = False, enter_set = True)
    b = Float(2, auto_set = False, enter_set = True)
    error = Float(0.0001, auto_set = False, enter_set = True)
    
    
    def __init__(self,**kw):
        super(RegulaFalsi,self).__init__(**kw)
        self._replot()
    
    @on_trait_change('a,b,error,f')
    def _replot(self):
        xdata, ydata = self.get_plot_data()
        self.plot_data.xdata = xdata
        self.plot_data.ydata = ydata
   
    """ define the function f(x) """
         
    def f(self,x):
        y = e**x - 2*x - 2
        return y   

    def get_plot_data(self):
        xdata = linspace(self.a,self.b,50)
        farray = frompyfunc( self.f, 1, 1 )
        ydata = array( farray( xdata ), dtype=float )
        return xdata, ydata
 
    results = Property( Float, depends_on = 'f, a, b, error' )        
    @cached_property
    def _get_results(self):
        """ returns approximation of x using Regula Falsi ..."""
        a = self.a
        b = self.b
        int = [a, b]
        f = self.f
        error = self.error
        step = []
        if f(a) * f(b) > 0:
            print "None or more than 1 root in selected interval"
        else:
            while abs(f(int[0] - f(int[0])*(int[1]-int[0])/(f(int[1])-f(int[0])))) > error:
                int.insert(0,int[0] - f(int[0])*(int[1]-int[0])/f((int[1])-f(int[0]))), int.pop(1)
                step.append('st')
            return [int[0], len(step)]
 
    value = Property( Float, depends_on = 'f, a, b, error' )
    @cached_property
    def _get_value(self):
        return self.results[0]

    steps = Property( Float, depends_on = 'f, a, b, error' )
    @cached_property
    def _get_steps(self):
        return self.results[1] 
    
    
    plot_data = Instance(MFnLineArray)
    def _plot_data_default(self):
        return MFnLineArray()
    
    trait_view = View(HSplit(VGroup(Group(
                                   Item('a'),
                                   Item('b'),
                                   Item('error'),
                                   Item('value', style='readonly',
                                        label = 'approx value of x'),
                                   Item('error', style='readonly',
                                        label = 'error estimation'),
                                   Item('steps', style='readonly',
                                        label = 'No. of steps'),
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
    idev = RegulaFalsi()
    idev.configure_traits()