
from enthought.traits.api import \
    HasTraits, Directory, List, Int, Float, Any, \
    on_trait_change, File, Constant, Instance, Trait, \
    Array, Str, Property, cached_property, WeakRef, \
    Dict, Button
    
from enthought.util.home_directory import \
    get_home_directory

from enthought.traits.ui.api import \
    View, Item, DirectoryEditor, TabularEditor, HSplit, Tabbed, VGroup, HSplit,\
    TableEditor, Group, ListEditor, VSplit
    
from enthought.traits.ui.table_column import \
    ObjectColumn
    
from enthought.traits.ui.menu import \
    OKButton, CancelButton
    
from enthought.traits.ui.tabular_adapter \
    import TabularAdapter
    
from numpy import array, loadtxt, arange, sqrt, zeros, arctan, sin, cos, ones_like,\
                vstack, savetxt, hstack, argsort, fromstring    

from math import pi
from string import split 
import os

#import csvPixel.py

from numpy import array, ones_like, arange

from scipy.io import read_array
from numpy import loadtxt, argmax, polyfit, poly1d, frompyfunc, sqrt

from enthought.traits.ui.table_filter \
    import EvalFilterTemplate, MenuFilterTemplate, RuleFilterTemplate, \
           EvalTableFilter

from enthought.pyface.dock.dock_sizer \
    import DockSizer, DockControl, DockRegion, DockStyle, DockSplitter, \
           no_dock_info, clear_window, features


class AsTbpFile( HasTraits ):
    '''
    Represent a single test specifying the design parameters.
    and access to the measured data.
    '''

  #  data_file = Str
    Names = []
    Units = []
    description = Str
    
    # open input file containing raw data:
    file = open( 'as_tbp.txt','r')
    
    n = file.read().split()
    
    lenght = len(n)
    i = 0
    while i < lenght:
        if n[i] == '#BEGINCHANNELHEADER':
#            desc = n[i+1]
#            j = 0
#            while j < len(desc):
#                description.append(desc[4+j])
#            print description
            Names.append(n[i+1])
            Units.append(n[i+3])
#            print 'N', Names
#            print n[i+3]
        i = i +1
  
    print 'Names =', Names
    print 'Units = ', Units

    print Names[0]
    
 #   filestr = file.read()
#    if filestr == '#BEGINCHANNELHEADER'
#    print 'str', filestr
    
    
#    lines = file.readlines()
#    for line in lines:
#        pair = line.split()
#        print pair
#        
#        lenght = len(pair)
#        if lenght > 0: 
#          if pair[0] == '#BEGINCHANNELHEADER' 
            
#            print '0', pair
#        if lenght >= 2:
#            print '2', pair 
          #print 'pair0 =', pair[0]
#          if  lenght < 2:
#            print 'pair1 =', pair[1]   
            
        
        
       # if xval == '200':
#        print 'pair0 =', pair[0]
#        lenght = len(pair)
#        print a
          #print 'pair1 =', pair[1]
        
        
#        if xval[0] == '200' and xval <> '' :
#          print xval[1]


#        Names = line
 #       if Names <> '':
 #         print 'N', Names[0], Names[1], Names[2]
#        if a_line == '200':  
 #         print '2'


    traits_view = View(
                       # Item( 'data_file', style = 'readonly' ),
                        Tabbed(
                            Item( 'Names' , show_label = False ),
                            Item( 'Units' , show_label = False),
                            scrollable = True,
                         ),
                      resizable = True,
                      scrollable = True,
                      height = 200,
                      width = 1000
                      )
   
isf = AsTbpFile() # data_file = 'as_tbp.txt' )
isf.configure_traits()
