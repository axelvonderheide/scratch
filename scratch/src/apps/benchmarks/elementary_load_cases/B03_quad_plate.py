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
    sin, cos, c_, arange, hstack, array, loadtxt

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

class SimQuadPlate( IBVModel ):

    implements( ISimModel )

    # discretization in x,y-direction:
    shape_xy = Int( 10,
                     ps_levels = (4, 20, 3 ) )
    
    # discretization in z-direction:
    shape_z = Int( 1,
                     ps_levels = (1, 2, 2 ) )


    length    = Float( 1.25 ) # edge length of the quadratic plate

    thickness = Float( 0.03 )
    
    nu = Float( 0.2 )

    fets = Instance( FETSEval, 
                     ps_levels = [ 'fe_linear', 
                                  'fe_quad_serendipity',
                                  'fe_quad_lagrange' ] )
    def _fets_default(self):
        return self.fe_quad_serendipity

    mats = Instance( MATS3DElastic )
    def _mats_default(self):
        return MATS3DElastic( E = 34000., nu = self.nu )
    
    fe_linear = Instance( FETSEval )    
    def _fe_linear_default(self):
        return FETS3D8H( mats_eval = self.mats ) 
        
    fe_quad_serendipity = Instance( FETSEval )    
    def _fe_quad_serendipity_default(self):
        return FETS3D8H20U( mats_eval = self.mats ) 

    fe_quad_lagrange = Instance( FETSEval )    
    def _fe_quad_lagrange_default(self):
        return FETS3D8H27U( mats_eval = self.mats ) 

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()
            
        u_center_top_z = U[ self.center_top_dof ][0,0,2]
        return array( [ u_center_top_z ], 
                        dtype = 'float_' )

    tloop = Property( depends_on = '+ps_levels' )
    @cached_property
    def _get_tloop(self):
 
        fets_eval = self.fets
        fets_eval.mats_eval.nu = self.nu
        
        # only a quarter of the plate is simulated due to symmetry:
        domain = FEGrid( coord_max = ( self.length / 2,  self.length / 2, self.thickness), 
                         shape   = ( self.shape_xy, self.shape_xy, self.shape_z ),
                         fets_eval = fets_eval )

        self.center_top_dof = domain[-1,-1,-1,-1,-1,-1].dofs

        # [MN] apply 1 kN at the center of the plate, i.e. 1/4 kN at the symmetric model)
        nodal_load = 0.25e-3 

        bc_symplane_yz  = BCSlice( var = 'u', value = 0., dims = [0], slice = domain[-1,:,:,-1,:,:] )
        bc_symplane_xz  = BCSlice( var = 'u', value = 0., dims = [1], slice = domain[:,-1,:,:,-1,:] )
        bc_support_000  = BCSlice( var = 'u', value = 0., dims = [0,1,2], slice = domain[0,0,0,0,0,0] )
        bc_center_load  = BCSlice( var = 'f', value = -nodal_load, dims = [2], slice = domain[-1,-1,-1,-1,-1,-1] )

        ts = TS(
                sdomain = domain,
                bcond_list =  [bc_symplane_yz, bc_symplane_xz, bc_support_000, 
                               bc_center_load
                               ],                                                  
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
        return [ SimOut( name = 'u_center_top_z', unit = 'm' ) ]

if __name__ == '__main__':

    sim_model = SimQuadPlate()

    ui = True
    
    if ui:
        print 'eval', sim_model.peval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_model )
        app.main()
    
    else:
                
        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()
        
        