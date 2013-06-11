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
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U    
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U    
from ibvpy.fets.fets2D5.fets2D58h import \
    FETS2D58H
from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic
        
from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic
    
from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage
    
from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear,\
    PhiFnStrainSoftening
            
from numpy import \
    sin, cos, c_, arange, hstack, array, loadtxt, ix_, unique

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

class WaistedShape( HasTraits ):
    '''
    Geo tranformation for a waisted tensile specimen.
    
    the upper face can vary between the value ly0 and ly1
    other dimensions are constant.
    
    dx is a  coefficient of a quadratic term defining the density
    of the mesh in the x direction. 
    dx = 0 means uniform mesh
    dx > 0 makes the mesh finer at the left hand side
    '''
    ly0 = Float(  1.0, enter_set = True , auto_set = False )
    ly1 = Float(  1.0, enter_set = True , auto_set = False )
    lx  = Float(  1.0, enter_set = True , auto_set = False )
    #lz  = Float( 0.03, enter_set = True , auto_set = False )
    lz  = Float( 1.0, enter_set = True , auto_set = False )
    dx  = Float( 0.0, enter_set = True , auto_set = False )
    
    CS = Property
    def _get_CS(self):
        return self.ly1 * self.lz
    
    def __call__( self, points ):
    
        ly0 = self.ly0
        ly1 = self.ly1
        lx = self.lx
        lz = self.lz
        # number of local grid points for each coordinate direction
        # values must range between 0 and 1
        xi, yi, zi = points[:,0], points[:,1], points[:,2]

        xi2 = xi * ( 1 + self.dx * ( 4 * xi**2 - 4 * xi ) )
        y = ( ly0 + ( ly1 - ly0 ) / lx * xi2 ) * yi
        x = lx * xi2
        z = lz * zi
        return c_[ x, y, z ]

class SimTensileTestWaisted( IBVModel ):
    '''Tensile test with wasted shape at the supports
    '''
    implements( ISimModel )

    # discretization in x,y-direction:
    shape_x = Int( 1 )
    
    shape_y = Int( 1 )

    # discretization in z-direction:
    shape_z = Int( 1 ) 
    
    geo_transform = Instance( WaistedShape )
    def _geo_transform_default(self):
        return WaistedShape()
        
    nu = Float( 0.25 )

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()
            
        u = U[ self.upper_right_dof ]
        return array( [ u ], 
                        dtype = 'float_' )

    tloop = Property() #  depends_on = '+ps_levels' )
    @cached_property
    def _get_tloop(self):
 
        u_max = 0.005
        u_max = 1.0
        rho = 0.03
        E_f = 70000
        E_m = 24000  
        E_c = E_m * ( 1 - rho ) + E_f * rho

        phi_fn = PhiFnStrainHardeningLinear( E_m = E_m, E_f = E_f, rho = rho,
                                             sigma_0 = 5.0,  
                                             alpha = 0.0, beta = 0.0, Elimit = 0.01 )

        # @todo: test elastic case: compare values of MATS2D5 with no damage (phi_fn = 1 (const.)).
        Epp = 100000.
        Efp = 1. 
        phi_fn = PhiFnStrainSoftening( Epp = Epp, Efp = Efp )

        mats = MATS2D5MicroplaneDamage( E = 1., 
                                        nu = 0.25,
                                        n_mp = 10,
                                        symmetrization = 'sum-type',
                                        model_version = 'compliance',
                                        phi_fn = phi_fn )         

#        # @todo: test elastic case: compare values of 2D5 and 3D
         # this is ok! FETS2D5 seams to work properly!
#        mats = MATS3DElastic( E = 1., nu = 0.25 )         

        fets =  FETS2D58H20U( mats_eval  = mats ) 
#        fets =  FETS2D58H( mats_eval  = mats ) 

        # only a quarter of the plate is simulated due to symmetry:
        domain = FEGrid( coord_max = ( 1, 1, 1 ), 
                         shape   = ( self.shape_x, self.shape_y, self.shape_z ),
                         geo_transform = self.geo_transform,
                         fets_eval = fets )

        self.center_top_dof = domain[-1,-1,-1,-1,-1,-1].dofs

        # bc:
        bc_symplane_yz_0 = BCSlice( var = 'u', value = 0., dims = [0], slice = domain[0,:,:,0,:,:])

        # apply displacement in x-direction (tension xx):
        bc_plane_yz_1_x    = BCSlice( var = 'u', value = u_max, dims = [0], slice = domain[-1,:,:,-1,:,:] )
        # apply displacement in y-direction (shear xy):
        bc_plane_yz_1_y    = BCSlice( var = 'u', value = u_max, dims = [1], slice = domain[-1,:,:,-1,:,:] )
        # apply displacement in y-direction (shear xz):
        bc_plane_yz_1_z    = BCSlice( var = 'u', value = -u_max, dims = [2], slice = domain[-1,:,:,-1,:,:] )

        bc_symplane_xz   = BCSlice( var = 'u', value = 0., dims = [1], slice = domain[:,0,:,:,0,:])
        bc_symplane_xy   = BCSlice( var = 'u', value = 0., dims = [2], slice = domain[:,:,0,:,:,0])

        self.upper_left_dof  = domain[  0, -1, -1,  0, -1, -1 ].dofs[0,0,0]
        self.upper_right_dof = domain[ -1, -1, -1, -1, -1, -1 ].dofs[0,0,0]

        self.right_dofs = unique( domain[-1,:,:,-1,:,:].dofs[:,:,0].flatten() )

        self.f_u_diagram = RTraceGraph( name = 'force - displ',
                                   var_y = 'F_int', idx_y_arr = self.right_dofs,
                                   var_x = 'U_k', idx_x = self.upper_right_dof,
                                   record_on = 'update',
                                   transform_x = 'x',
                                   transform_y = 'y' )   
        ts = TS(
                sdomain = domain,
                bcond_list =  [ bc_symplane_yz_0, 
                                bc_plane_yz_1_x , 
#                                bc_plane_yz_1_y , 
                                bc_plane_yz_1_z , 
                                bc_symplane_xz, 
                                bc_symplane_xy ],                                                  
                                
                rtrace_list = [ 
                               self.f_u_diagram, 
                             RTraceDomainListField(name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True),
                             RTraceDomainListField(name = 'Stress' ,
                                            var = 'sig_app', idx = 0, warp = True, 
                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'IStress' ,
#                                            position = 'int_pnts',
#                                            var = 'sig_app', idx = 0, warp = True, 
#                                            record_on = 'update'),
                             RTraceDomainListField(name = 'Strain' ,
                                            var = 'eps_app', idx = 0, warp = True, 
                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'IStrain' ,
#                                            position = 'int_pnts',
#                                            var = 'eps_app', idx = 0, warp = True, 
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'Damage' ,
#                                            var = 'omega_mtx', idx = 0, warp = True, 
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'IDamage' ,
#                                            position = 'int_pnts',
#                                            var = 'omega_mtx', idx = 0, warp = True, 
#                                            record_on = 'update'),                              
                                            ]             
                
                )
         
        # Add the time-loop control
        tloop = TLoop( tstepper = ts, tolerance = 1e-4, KMAX = 200, RESETMAX = 0, 
                       tline  = TLine( min = 0.0, step = 1.0, max = 1.0 ) )                   
#                       tline  = TLine( min = 0.0, step = 0.03, max = 1.0 ) )                   
        
        return tloop

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_right', unit = 'm' ) ]

if __name__ == '__main__':

    sim_model = SimTensileTestWaisted( geo_transform = WaistedShape( dx = 0.0 ) )

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()
        #sim_model.f_u_diagram.configure_traits()
    
    elif do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_model )
        app.main()
    
    elif do == 'ps':
        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()
        