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
# Thanks for using Simvisage open source!
#
# Created on Nov 18, 2011 by: matthias

from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Bool, File, Array
    
import numpy as np

class Cnstr_Entity(HasTraits):
    '''
     Main Constrain-Object 
    '''
    dtype = None
    # left-hand side coefficients of the constraint equations 
    cnstr_lhs = Array(value = [], dtype = float)
    # right-hand side values of the constraint equations
    cnstr_rhs = Array(value = [], dtype = float)
    
    
class Cnstr_Line(Cnstr_Entity):
    '''
     Constrains projected on a line
    '''
    dtype = 'line'
    # directional Vector
    a_dir_v = Array(value =[1.0, 0.0, 0.0], dtype = float)
    
    # Base Vector
    base_v = Array(value = [0.0, 0.0, 0.0], dtype = float)
    
    
class Cnstr_Plane(Cnstr_Entity):
    '''
     Constrains projected on a plane
    '''
    dtype = 'plane'
    # directonal Vectors
    #a_dir_v = Array(value =[1.0, 0.0, 0.0], dtype = float)
    #b_dir_v = Array(value =[0.0, 1.0, 0.0], dtype = float)
    
    # normal of the plane
    normal = Array(value = [0.0, 0.0, 1.0], dtype = float)
    
    # Base Vector
    base_v = Array(value = [0.0, 0.0, 0.0], dtype = float)
        