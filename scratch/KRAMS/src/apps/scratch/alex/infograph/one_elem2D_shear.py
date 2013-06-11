'''
Example of a tensile test using a one element discretization
'''

from enthought.traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from ibvpy.api import \
    TStepper as TS, TLoop, TLine, \
    IBVPSolve as IS, DOTSEval, \
    RTraceGraph, RTraceDomainListField,  \
    BCDof, BCDofGroup,  BCSlice
    
from ibvpy.fets.fets2D.fets2D4q import \
    FETS2D4Q
    
from ibvpy.fets.fets2D.fets2D9q import \
    FETS2D9Q
    
from ibvpy.mesh.fe_grid import \
    FEGrid
    
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic
    
from mathkit.geo.geo_ndgrid import \
    GeoNDGrid
    
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint
    
from numpy import \
    sin, cos, c_, arange, hstack, array, loadtxt

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimModel, SimPStudy, SimPar, SimOut, ISimModel

class SimQuadPlate( HasTraits ):

    implements( ISimModel )
    
    shape = Int( 10,
                   ps_levels = (1, 40, 4 ) )
    
    edge_length  = Float( 1.0,
                          ps_levels = (1, 10, 9) )

    nu = Float( 0.2,
                ps_levels = (0.0, 0.5, 10 ) )

    def get_sim_outputs( self ):
        return [ SimOut( name = 'u_upper_right', unit = 'm' ) ]
        
    def eval(self):
        
        mats = MATS2DElastic( E = 35000., nu = self.nu, stress_state  = "plane_strain" )
        
        fets_eval = FETS2D4Q( mats_eval = mats ) 
        
        domain = FEGrid( coord_max = ( self.edge_length,  self.edge_length, 0.), 
                         shape   = ( int( self.shape ), int( self.shape ) ),
                         fets_eval = fets_eval )

        upper_right_corner_dof = domain[-1,-1,-1,-1].dofs
        
        ts = TS(
                sdomain = domain,
        
        #        # simple shear test: clamped at left side loaded at right side
        #        bcond_list = [BCSlice( var = 'u', value = 0., dims = [0,1], slice = domain[0, 0, 0, :] ),
        #                      BCSlice( var = 'u', value = 0., dims = [0]  , slice = domain[0, 0,-1, :] ),
        #                      BCSlice( var = 'f', value = 1.0, dims = [1] , slice = domain[0, 0,-1, :] )],
        
                # shear test: fixed at 000, one support in x-direction at left top; load at top right
                bcond_list = [BCSlice( var = 'u', value = 0., dims = [0,1], slice = domain[0, 0, 0, 0] ),
                              BCSlice( var = 'u', value = 0., dims = [0]  , slice = domain[0, 0, 0,-1] ),
                              BCSlice( var = 'f', value = 1.0, dims = [1] , slice = domain[0, 0,-1,-1] )],

                rtrace_list = [ 
                                
                             RTraceDomainListField(name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True),
                             RTraceDomainListField(name = 'Stress' ,
                                            var = 'sig_app', idx = 0, warp = True, 
                                            record_on = 'update'),
                            ]             
                )
         
        # Add the time-loop control
        tloop = TLoop( tstepper = ts,
                       tline  = TLine( min = 0.0,  step = 1., max = 1.0 ) )                   
        
        tloop.eval()
        return tloop

    def peval(self):
        u = self.eval() 
        u_upper_right = u[ upper_right_corner_dof ][0,0,0]
        return array( [ u_upper_right ], dtype = 'float_' )

if __name__ == '__main__':

    sim_ps = SimArray( model = SimQuadPlate() )
    sim_ps.configure_traits()