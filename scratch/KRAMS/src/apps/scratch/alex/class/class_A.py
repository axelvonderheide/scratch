'''
Created on Feb 17, 2010

@author: alexander
'''


from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from traits.editors.mpl_figure_editor import MPLFigureEditor
from matplotlib.figure import Figure
    
import os

import csv

from numpy import array, fabs, where, copy, ones

from scipy.io import read_array
from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

#-- Tabular Adapter Definition -------------------------------------------------

from string import replace
from os.path import exists



class reinforcement( HasTraits ):
    type = Str('MAG')
    a_tex = 0
    
    def set_a( self ):
        if self.type == '2D':
            self.a_tex = 10
    
        
if __name__ == '__main__':
    reinf = reinforcement(type = '2D')
    print reinf.type
    print reinf.a_tex
    
    