'''
Created on May 19, 2011

@author: kelidas
'''


from quaducom.ymb_micro.ymb_data import YMBCutData, YMBSegmentData, YMBSlider, YMBSource, yarn_list, var_dict
from numpy import linspace, array, searchsorted

from enthought.traits.api import HasTraits, Property, cached_property, Event, \
    Instance, Int, on_trait_change, Bool, Str, Tuple, List, Array, Enum, Trait, Range
from enthought.traits.ui.api import Item, View, Group, HGroup, HSplit
import matplotlib.pyplot as plt

source = YMBSource( yarn_type = 'MAG' )
data = YMBSegmentData( source = source, cf_limit = 1 )


slack = ( data.bond_free_length - data.bond_free_length_x ) / data.bond_free_length_x
slack_m = slack[:, 7]
slack_f = slack.flatten()

#plt.hist( slack_m[slack_m > -1], bins = 100 )
plt.hist( slack_f[slack_f > -1], bins = 1000 )
plt.show()








