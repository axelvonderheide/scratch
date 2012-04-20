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
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U

from ibvpy.mesh.fe_grid import \
    FEGrid
    
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic
    
from mathkit.geo.geo_ndgrid import \
    GeoNDGrid
    
from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint
    
from numpy import \
    sin, cos, c_, arange, hstack, array, loadtxt, ix_

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

class SimQuadPlate( IBVModel ):

    implements( ISimModel )
    
    shape_x = Int( 20,
                     ps_levels = (1, 40, 4 ) )
    
    shape_y = Int( 2,
                     ps_levels = (1, 4, 4 ) )

    height  = Float( 1.0 )

    length = Float( 15.0 )

    nu = Float( 0.2,
                ps_levels = (0.0, 0.5, 3 ) )

    fets = Instance( FETSEval, 
                     ps_levels = [ 'fe_type_linear', 
                                  'fe_type_quadratic',
                                  'fe_type_quadratic_lag' ], transient = True )
    def _fets_default(self):
        #return self.fe_type_linear
        return self.fe_type_quadratic_lag

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_upper_right_x', unit = 'm' ),
                 SimOut( name = 'u_upper_right_y', unit = 'm' ) ]


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
    
class SimQuadPlate2D( SimQuadPlate ):

    mats = Instance( MATS2DElastic, transient = True )
    def _mats_default(self):
        return MATS2DElastic( E = 35000., nu = self.nu, stress_state  = "plane_stress" )
    
    fe_type_linear = Instance( FETSEval, transient = True )    
    def _fe_type_linear_default(self):
        return FETS2D4Q( mats_eval = self.mats ) 

    fe_type_quadratic = Instance( FETSEval, transient = True )    
    def _fe_type_quadratic_default(self):
        return FETS2D9Q( mats_eval = self.mats ) 
        
    fe_type_quadratic_lag = Instance( FETSEval, transient = True )    
    def _fe_type_quadratic_lag_default(self):
        return FETS2D4Q9U( mats_eval = self.mats ) 
        
    @cached_property
    def _get_tloop(self):
        
        fets_eval = self.fets
        fets_eval.mats_eval.nu = self.nu
        
        domain = FEGrid( coord_max = ( self.length,  self.height, 0.), 
                         shape   = ( int( self.shape_x ), int( self.shape_y ) ),
                         fets_eval = fets_eval )

        self.upper_right_corner_dof = domain[-1,-1,-1,-1].dofs

        line_load = -1.0

        old = False
        if old:
            elem_length = self.length / float( self.shape_x )
    
            force_bc = []
            if self.fets == self.fe_type_linear:
                nodal_load = 1./2. * line_load * elem_length
                force_bc = [ BCSlice( var = 'f', value = nodal_load, dims = [1] , slice = domain[:,-1,:,-1] ) ]
            elif self.fets == self.fe_type_quadratic:
                nodal_corner_elem_load = 1./6. * line_load * elem_length
                nodal_edge_elem_load = 4./6. * line_load * elem_length
                elem_corner_nodes = domain[:, -1, [0,-1],-1 ]
                elem_edge_nodes   = domain[:, -1, 1:-1  ,-1 ]
                force_bc = [BCSlice( var = 'f', value = nodal_corner_elem_load, dims = [1] , slice = elem_corner_nodes ),
                            BCSlice( var = 'f', value = nodal_edge_elem_load, dims = [1] , slice = elem_edge_nodes )]            
            else:
                raise ValueError, 'Element type not found'
        else:
            force_bc = [ BCSlice( var = 'f', value = line_load, dims = [1] , slice = domain[:,-1,:,-1] ) ]
        
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

class SimQuadPlate3D( SimQuadPlate ):

    mats = Instance( MATS3DElastic, transient = True )
    def _mats_default(self):
        return MATS3DElastic( E = 35000., nu = self.nu )
    
    fe_type_linear = Instance( FETSEval, transient = True )    
    def _fe_type_linear_default(self):
        return FETS3D8H( mats_eval = self.mats ) 
        
    fe_type_quadratic = Instance( FETSEval, transient = True )    
    def _fe_type_quadratic_default(self):
        return FETS3D8H20U( mats_eval = self.mats ) 

    fe_type_quadratic_lag = Instance( FETSEval, transient = True )    
    def _fe_type_quadratic_lag_default(self):
        return FETS3D8H27U( mats_eval = self.mats ) 
        
    @cached_property
    def _get_tloop(self):

        fets_eval = self.fets
        fets_eval.mats_eval.nu = self.nu        
        
        domain = FEGrid( coord_max = ( self.length,  self.height, 1.), 
                         shape   = ( self.shape_x, self.shape_y, 1 ),
                         fets_eval = fets_eval )

        self.upper_right_corner_dof = domain[-1,-1,-1,-1,-1,-1].dofs

        surface_load = - 1.0

        if False: # to be removed
            # NOTE: the depth of the beam, i.e. the length in z-direction, is set to "1"
            elem_length = self.length / float( self.shape_x )
            force_bc = []
            if self.fets == self.fe_type_linear:
                nodal_load = 1./4. * surface_load * elem_length * 1. 
                upper_surface = domain[:,-1,:,:,-1,:]
                force_bc = [ BCSlice( var = 'f', value = nodal_load, dims = [1] , slice = upper_surface ) ]
            elif self.fets == self.fe_type_quadratic:
                nodal_corner_elem_load =  -1./12. * surface_load * elem_length * 1.
                nodal_edge_elem_load   =   1./3. * surface_load * elem_length * 1.
    
                ci = ix_( [0,-1], [0,-1] )
                elem_corner_nodes     = domain[:, -1, :,  ci[0], -1, ci[1] ]
                elem_edge_nodes_part1 = domain[:, -1, :,  1:-1 , -1, [0,-1] ]
                elem_edge_nodes_part2 = domain[:, -1, :, [0,-1], -1,  1:-1  ]
                force_bc = [BCSlice( var = 'f', value = nodal_corner_elem_load, dims = [1] , slice = elem_corner_nodes ),
                            BCSlice( var = 'f', value = nodal_edge_elem_load,   dims = [1] , slice = elem_edge_nodes_part1 ),
                            BCSlice( var = 'f', value = nodal_edge_elem_load,   dims = [1] , slice = elem_edge_nodes_part2 ),            
                            ]
            elif self.fets == self.fe_type_quadratic_lag:
                nodal_corner_elem_load =  1./36. * surface_load * elem_length * 1.
                nodal_edge_elem_load   =  4./36. * surface_load * elem_length * 1.
                nodal_middle_elem_load = 16./36. * surface_load * elem_length * 1.
    
                ci = ix_( [0,-1], [0,-1] )
                elem_corner_nodes     = domain[:, -1, :,  ci[0], -1, ci[1] ]
                elem_edge_nodes_part1 = domain[:, -1, :,  1:-1 , -1, [0,-1] ]
                elem_edge_nodes_part2 = domain[:, -1, :, [0,-1], -1,  1:-1  ]
                elem_middle_nodes     = domain[:, -1, :,  1:-1 , -1,  1:-1  ]
                force_bc = [BCSlice( var = 'f', value = nodal_corner_elem_load, dims = [1] , slice = elem_corner_nodes ),
                            BCSlice( var = 'f', value = nodal_edge_elem_load,   dims = [1] , slice = elem_edge_nodes_part1 ),
                            BCSlice( var = 'f', value = nodal_edge_elem_load,   dims = [1] , slice = elem_edge_nodes_part2 ),            
                            BCSlice( var = 'f', value = nodal_middle_elem_load, dims = [1] , slice = elem_middle_nodes )]
                
            else:
                raise ValueError, 'Element type not found'
        else:
            force_bc = [ BCSlice( var = 'f', value = surface_load, 
                                  dims = [1] , slice = domain[:,:,-1,:,:,-1] ) ]
        
        self.ts = TS(
                sdomain = domain,
                # 
                bcond_list = [ BCSlice( var = 'u', value = 0., dims = [0], slice = domain[0, :, :, 0, :, :] ),
                              BCSlice( var = 'u', value = 0., dims = [1], slice = domain[0, 0, :, 0, 0, :] ),                            
                              BCSlice( var = 'u', value = 0., dims = [2], slice = domain[0, 0, 0, 0, 0, 0] )                               
                              ] + force_bc,
                rtrace_list = [ 
                                
                             RTraceDomainListField(name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True),
                             RTraceDomainListField(name = 'Stress' ,
                                            var = 'sig_app', idx = 0, warp = True, 
                                            record_on = 'update'),
                            ]             
                )
         
        # Add the time-loop control
        tloop = TLoop( tstepper = self.ts,
                       tline  = TLine( min = 0.0,  step = 1., max = 1.0 ) )                   
        
        return tloop

if __name__ == '__main__':

    sim_model = SimQuadPlate3D()

    do = 'ps'

    if do == 'eval':

        print 'default response', sim_model.peval()

    if do == 'ui':
        sim_model.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_model )
        app.main()
    
    elif do == 'ps':
                
        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()

    elif do == 'pickle':
        
        import pickle
        filename = '/tmp/sim.pickle'
        file = open(filename,'w')
        pickle.dump( sim_model, file )
        file.close()
        file = open(filename,'r')
        sm = pickle.load( file )
        file.close()
