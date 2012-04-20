

from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate, Button, \
     Interface

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

import enthought.units as units
from enthought.units.angle import degree, radian

from math import pi as Pi, cos, sin, exp


class MFnPolar180( HasTraits ):
    
    ralpha_range = Range( 0., 90., 90.,
                  label = 'Reinforcement range',
                  auto_set = False)
    ralpha = Range( 0., 180., 0.,
                  label = 'Reinforcement direction',
                  auto_set = False)
    salpha = Range( 0., 90., 45.,
                  label = 'Transition range',
                  auto_set = False)
    res_integ = Range(0.,1.,0.4, label = 'Residual integrity')
    fn = Instance( MFnLineArray )
    def _fn_default(self):
        xdata, ydata = self._get_fn_data()
        return MFnLineArray( xdata = xdata, ydata = ydata ) 

    @on_trait_change('ralpha_range,ralpha,salpha,res_integ')
    def _set_fn_data(self):
        self.fn.xdata, self.fn.ydata = self._get_fn_data()
    
    def _get_fn_data(self):
        ralpha       = units.convert( self.ralpha, degree, radian )
        ralpha_range = units.convert( self.ralpha_range, degree, radian )
        salpha       = units.convert( self.salpha, degree, radian )

        r_min = ralpha - ralpha_range / 2.
        r_max = ralpha + ralpha_range / 2.

        if r_min < 0:  s_min = r_min + Pi
        else:          s_min = r_min
        if r_max > Pi: s_max = r_max - Pi
        else:          s_max = r_max
        
        step_size = salpha
        rinteg = self.res_integ
        if s_min < s_max:
            self.single_interval = True
            return ([0.0, s_min-(step_size), s_min,
                    s_max, s_max+(step_size), Pi], 
                    [ 0.0, 0.0, rinteg, rinteg, 0.0, 0.0 ] )
        else:
            self.single_interval = False
            return ([0.0, s_max, s_max+(step_size),
                    s_min-(step_size), s_min, Pi],
                    [ rinteg, rinteg, 0.0, 0.0, rinteg, rinteg ])

    traits_view = View(Item('ralpha'),
                       Item('ralpha_range'),
                       Item('salpha'),
                       Item('res_integ'),
                       Item('fn@', show_label = False ),
                       resizable = True)
            
    get_value = Delegate('fn')
    
    def __call__(self, x):
        return self.get_value(x)

if __name__ == '__main__':
    mfn = MFnPolar180()
    mfn.configure_traits()
