
'''
Created on Jul 30, 2010

@author: abach
'''

from enthought.traits.api import \
    HasTraits, Float, Array, implements, Property, cached_property, Instance, \
    Int, List, Bool, HasTraits

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
    TLine, BCDof, BCDofGroup, BCSlice, IBVModel


from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic


from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets2D.fets2D4q import \
    FETS2D4Q
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U
from ibvpy.mesh.fe_grid import \
    FEGrid
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray

from numpy import \
    array, tensordot, dot, zeros, c_, ix_, shape, \
    cos, sin, arctan, where, abs, all, any, diag

from ibvpy.mats.mats3D.mats3D_tensor import map3d_sig_eng_to_mtx
from math import sqrt, asin, acos, pi as Pi

# Interpolation
from scipy.interpolate import Rbf

from simiter.sim_pstudy import ISimModel, SimOut, SimPStudy


class hinge_forces( IBVModel ):

    implements( ISimModel )

    fe_grid_1 = Property( Instance( FEGrid ) )
    @cached_property
    def _get_fe_grid_1( self ):
        return FEGrid( coord_min = ( 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0 ),
                       shape = ( 5, 5 ),
                       fets_eval = FETS2D4Q( mats_eval = MATS2DElastic( E = 10e3, nu = 0.2, ) ) )
#                                                                        initial_strain = self.temperature_strain_z ) ) )


    initial_strain_roof = Bool( False )
    initial_strain_col = Bool( False )
    alpha = Float( 1e-5 )
    t_up = Float( 10. )
    t_lo = Float( 10. )
    def temperature_strain_z( self, X_pnt, x_pnt ):
        alpha, t_up, t_lo = self.alpha, self.t_up, self.t_lo
        delta_t = t_lo + ( t_up - t_lo ) * x_pnt[1]
        epsilon_0 = alpha * delta_t
        # return the initial volumetric strain tensor with n_dims  
        return diag( [ epsilon_0 for i in range( 2 ) ] )


    def geo_transfom_grid2( self, points ):
        xi, yi = points[:, 0], points[:, 1]
        return c_[xi + 1.0, yi]

    fe_grid_2 = Property( Instance( FEGrid ) )
    @cached_property
    def _get_fe_grid_2( self ):
        return FEGrid( coord_min = ( 0.0, 0.0 ),
                       coord_max = ( 1.0, 1.0 ),
                       geo_transform = self.geo_transfom_grid2,
                       shape = ( 5, 5 ),
                       fets_eval = FETS2D4Q( mats_eval =
                                             MATS2DElastic( E = 10e3,
                                                            nu = 0.2, ) ) )
                                                            #initial_strain = self.temperature_strain_z ) ) )

    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        F_int = self.tloop.tstepper.F_int
        F_ext = self.tloop.tstepper.F_ext

        F_int_slice_l_1 = F_int[self.slice_edge_left_1]
        F_int_slice_r_1 = F_int[self.slice_edge_right_1]

        F_ext_slice_l_1 = F_ext[self.slice_edge_left_1]
        F_ext_slice_r_1 = F_ext[self.slice_edge_right_1]

        print"element1"
        print"internal-left", F_int_slice_l_1
        print"internal-right", F_int_slice_r_1

        print"external-left", F_ext_slice_l_1
        print"external-right", F_ext_slice_r_1

        F_int_slice_l_2 = F_int[self.slice_edge_left_2]
        F_int_slice_r_2 = F_int[self.slice_edge_right_2]

        F_ext_slice_l_2 = F_ext[self.slice_edge_left_2]
        F_ext_slice_r_2 = F_ext[self.slice_edge_right_2]

        print"element2"
        print"internal-left", F_int_slice_l_2
        print"internal-right", F_int_slice_r_2

        print"external-left", F_ext_slice_l_2
        print"external-right", F_ext_slice_r_2

        from numpy import sum
        print"summe", sum( F_int_slice_l_2[:, :-1, 0] ) + F_int_slice_l_2[-1, -1, 0]

#        print "shape F_int", shape( F_int )
#        print "slice_edge", self.slice_edge
#        from numpy import sum
#        print "F_int_slice", F_int_slice
#        print "F_int sum", sum( F_int_slice[:, :, 0] )





        return array( 1.0 )

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [  SimOut( name = 'U_z', unit = 'm' ), ]


    sig_app = Property( Instance( RTraceDomainListField ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_sig_app( self ):
        return RTraceDomainListField( name = 'sig_app' ,
#                                      position = 'int_pnts',
                                      var = 'sig_app',
                                      record_on = 'update', )

    u = Property( Instance( RTraceDomainListField ), depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_u( self ):
        return RTraceDomainListField( name = 'displacement' ,
                                      var = 'u', warp = True,
                                      record_on = 'update', )

    rtrace_list = List
    def _rtrace_list_default( self ):
        return [ self.sig_app, self.u]

    tline = Instance( TLine )
    def _tline_default( self ):
        return TLine( min = 0.0, step = 1.0, max = 1.0 )

    tloop = Property( depends_on = '+ps_levels, +input' )
    @cached_property
    def _get_tloop( self ):

        fe_grid_1 = self.fe_grid_1
        fe_grid_2 = self.fe_grid_2

        self.slice_edge_left_1 = self.fe_grid_1[0, :, 0, :].dofs
        self.slice_edge_right_1 = self.fe_grid_1[-1, :, -1, :].dofs
        self.slice_edge_left_2 = self.fe_grid_1[0, :, 0, :].dofs
        self.slice_edge_right_2 = self.fe_grid_1[-1, :, -1, :].dofs


        bc_force_list = [BCSlice( var = 'f'  , dims = [0],
                                  slice = fe_grid_2[-1, :, -1, :],
                                  value = 1.0 )]

        bc_disp_list = [BCSlice( var = 'u'  , dims = [0],
                                 slice = fe_grid_1[0, : , 0, :],
                                 value = 0 ),
                        BCSlice( var = 'u'  , dims = [1],
                                 slice = fe_grid_1[0, 0, 0, 0],
                                 value = 0 ),
#                        BCSlice( var = 'u'  , dims = [0],
#                                 slice = fe_grid_2[-1, :, -1, :],
#                                 value = 0 ),
#                        BCSlice( var = 'u'  , dims = [1],
#                                 slice = fe_grid_2[:, 0, :, 0],
#                                 value = 0 ),
                        ]


        bc_link_list = [BCSlice( var = 'u'  , dims = [0, 1],
                             slice = fe_grid_1[-1 , : , -1 , : ],
                             link_slice = fe_grid_2[0 , : , 0 , : ], link_coeffs = [1.0],
                             value = 0.0 ) ]

        ts = TS( sdomain = [fe_grid_1, fe_grid_2],
                 dof_resultants = True,
                 bcond_list = bc_force_list +
                              bc_disp_list +
                              bc_link_list,
                 rtrace_list = self.rtrace_list
               )

        tloop = TLoop( tstepper = ts,
                       tolerance = 1e-4,
                       tline = self.tline )

        return tloop

if __name__ == '__main__':

    sim_model = hinge_forces ()
#    print sim_model

    do = 'ui'

    if do == 'eval':
        print 'eval', sim_model.peval()

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
        file = open( filename, 'w' )
        pickle.dump( sim_model, file )
        file.close()
        file = open( filename, 'r' )
        sm = pickle.load( file )
        file.close()

