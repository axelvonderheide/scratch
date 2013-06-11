#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Created on May 25, 2009 by Rostislav Chudoba
#
#-------------------------------------------------------------------------------

from enthought.traits.api \
    import HasTraits, HasPrivateTraits, Color, Str, Int, Float, Enum, List, \
           Bool, Instance, Any, Font, Event, Property, Interface, \
           on_trait_change, cached_property, implements, Color, Tuple
           
#-------------------------------------------------------------------------------
#  'MFnPlotAdapter' class:
#-------------------------------------------------------------------------------

class MFnPlotAdapter ( HasPrivateTraits ):
    
    """ The base class for adapting function implementation in order to plot
    it in a plotting toolkit (either Chaco or Matplotlib)
    """
    
    # label on the x axis (only used when var_x empty)
    #
    label_x = Str( 'displacement' )
    
    # label on the y axis (only used when var_y empty)
    label_y = Str( 'force' )

    # title of the diagram
    #
    title = Str('FORCE vs DISPLACEMENT')
    
    # color of the title
    title_color = Str('black')
    
    # when specified take the label of the axis from the value 
    # of the variable trait in object
    #
    var_x = Str('')
    
    # when specified take the label of the axis from the value 
    # of the variable trait in object
    #
    var_y = Str('')
    
    # label for line in legend
    #
    line_label = Str('line1')
    
    # @todo - further attributes to be made available
    
    # limits of number positions for switching into scientific notation 1eXX
    scilimits = Tuple(-3.,4.)
    
    # Plot properties
    
    # line properties can be also specified as a tuple
    # of RGB or RGBA values from 0-1
    line_color = Str( "black" )
    
    # linewidth
    linewidth = Float(2.0)
    
    # linestyle, e.g. 'dashed'
    line_style = Str('solid')
    
    # color of the plotting background
    bgcolor = Str( "white" )

  
    # Border, padding properties
    border_visible = Bool( False )
    
    border_width =  Int( 0 )
    
    padding_bg_color = Str( 'white' )
    
    # labels in legend
    
    legend_labels = Tuple('label1','label2','label3')
                