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
    IBVModel, DOTSEval, \
    RTraceGraph, RTraceDomainListField,  \
    BCDof, BCDofGroup,  BCSlice
    
from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets2D.fets2D4q import \
    FETS2D4Q
from ibvpy.fets.fets2D.fets2D9q import \
    FETS2D9Q
from ibvpy.fets.fets2D.fets2D4q9u import \
    FETS2D4Q9U

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
    SimPStudy, SimPar, SimOut, ISimModel

class SimQuadPlate( IBVModel ):

    implements( ISimModel )
    
    shape = Int( 10,
                   ps_levels = (1, 40, 4 ) )
    
    edge_length  = Float( 1.0,
                          ps_levels = (1, 10, 9) )

    nu = Float( 0.2,
                ps_levels = (0.0, 0.5, 10 ) )

    fets = Instance( FETSEval, 
                     ps_levels = [ 'fe_type_linear', 'fe_type_quadratic' ] )
    def _fets_default(self):
        return self.fe_type_linear

    mats = Instance( MATS2DElastic )
    def _mats_default(self):
        return MATS2DElastic( E = 35000., nu = self.nu, stress_state  = "plane_strain" )
    
    fe_type_linear = Instance( FETSEval )    
    def _fe_type_linear_default(self):
        return FETS2D4Q( mats_eval = self.mats ) 
        
    fe_type_quadratic = Instance( FETSEval )    
    def _fe_type_quadratic_default(self):
        return FETS2D4Q9U( mats_eval = self.mats ) 

    def _tloop_default(self):
        return self._new_tloop()

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()
            
        u_upper_right_x = U[ self.upper_right_corner_dof ][0,0,0]
        u_upper_right_y = U[ self.upper_right_corner_dof ][0,0,1]
        return array( [ u_upper_right_x, u_upper_right_y ], 
                        dtype = 'float_' )

    tloop = Property( depends_on = '+ps_levels' )
    def _get_tloop(self):
        
        fets_eval = self.fets
        fets_eval.mats_eval.nu = self.nu
        
        domain = FEGrid( coord_max = ( self.edge_length,  self.edge_length, 0.), 
                         shape   = ( int( self.shape ), int( self.shape ) ),
                         fets_eval = fets_eval )

        self.upper_right_corner_dof = domain[-1,-1,-1,-1].dofs
        
        elem_length = self.edge_length / float( self.shape )
        line_load = 2.0

        force_bc = []
        if self.fets == self.fe_type_linear:
            nodal_load = 1./2. * line_load * elem_length        
            force_bc = [ BCSlice( var = 'f', value = nodal_load, dims = [0] , slice = domain[-1, :,-1, :] ) ]
        elif self.fets == self.fe_type_quadratic:
            nodal_corner_elem_load = 1./6. * line_load * elem_length
            nodal_edge_elem_load = 4./6. * line_load * elem_length
            elem_corner_nodes = domain[-1, :,-1, [0,-1] ]
            elem_edge_nodes   = domain[-1, :,-1, 1:-1 ]
            force_bc = [BCSlice( var = 'f', value = nodal_corner_elem_load, dims = [0] , slice = elem_corner_nodes ),
                        BCSlice( var = 'f', value = nodal_edge_elem_load, dims = [0] , slice = elem_edge_nodes )]            
        else:
            raise ValueError, 'Element type not found'
        
        ts = TS(
                sdomain = domain,
                # 
                bcond_list = [ BCSlice( var = 'u', value = 0., dims = [0,1], slice = domain[0, :, 0, :] ) ] + force_bc,
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
        
        return tloop

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_upper_right_x', unit = 'm' ),
                 SimOut( name = 'u_upper_right_y', unit = 'm' ) ]

if __name__ == '__main__':

    sim_mod = SimQuadPlate()

    ui = False
    
    if ui:
        sim_mod.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_mod )
        app.main()
    
    else:
        sim_ps = SimArray( model = sim_mod )
        sim_ps.configure_traits()