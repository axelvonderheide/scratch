'''
Created on Apr 22, 2010

@author: alexander
'''
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button, Bool, Enum, Event, implements, \
    DelegatesTo

from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, VGroup, \
    TableEditor, EnumEditor, Handler, FileEditor, VSplit, Group, \
    InstanceEditor, HGroup, Spring
    
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

from numpy import \
    array, fabs, where, copy, ones, linspace, ones_like, hstack

from scipy.io import read_array
from numpy import \
    loadtxt, argmax, polyfit, poly1d, frompyfunc, dot

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

from enthought.traits.ui.api import TextEditor 
float_editor_9g = TextEditor( evaluate=float, format_str="%.9g", enter_set=True, auto_set=False )
float_editor_0f = TextEditor( evaluate=float, format_str="%.0f", enter_set=True, auto_set=False )
float_editor_3f = TextEditor( evaluate=float, format_str="%.3f", enter_set=True, auto_set=False )
float_editor_4f = TextEditor( evaluate=float, format_str="%.4f", enter_set=True, auto_set=False )

#-- Tabular Adapter Definition -------------------------------------------------

from string import replace
from os.path import exists

#-----------------------------------------------------------------------------------
# ExDesignReader
#-----------------------------------------------------------------------------------
from enthought.traits.ui.file_dialog  \
    import open_file, FileInfo, TextInfo, ImageInfo

from enthought.traits.ui.api \
    import View, Item, TabularEditor
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter

from mathkit.array.smoothing import smooth

from promod.matdb.trc.textile_cross_section \
    import TextileCrossSection
    
from promod.matdb.trc.textile_fabric_layout \
    import TextileFabricLayout 
  
from promod.matdb.trc.concrete_mixture \
    import ConcreteMixture
    
from promod.matdb.trc.composite_cross_section \
    import CompositeCrossSection
  

def get_unit( trait ):
    '''Extract the string stored for 'unit' in the trait metadata 
    '''
    dict = trait.__dict__
    return dict['_metadata']['unit'] 


def convert_unit( unit_old, unit_new  ):
    '''Returns a factor for converting 'unit_old' to 'unit_new' 
    unit_old  - given unit to be converted
    unit_new - converted unit (target unit) 
    '''
    m_dict = {'m':1, 'dm':10, 'cm':100, 'mm':1000,}
    return m_dict[ unit_new ]/m_dict[ unit_old ] 


class A( HasTraits ):
    '''Read the data from the directory
    '''
    b = Int(3, unit = 'm')
#    u = b.__dict__
#    unit_old = u['_metadata']['unit']
#    print 'u', u['_metadata']['unit']
    u = get_unit(b)
    f = convert_unit(u,'mm')
    print 'u', u
    print 'f', f

    
if __name__ == '__main__':
    A = A()
    print 'b', A.__getstate__()
    A.configure_traits()

