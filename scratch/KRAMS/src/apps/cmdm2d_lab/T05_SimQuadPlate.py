
from enthought.traits.api import \
    Array, Bool, Enum, Float, HasTraits, \
    Instance, Int, Trait, Str, Enum, \
    Callable, List, TraitDict, Any, Range, \
    Delegate, Event, on_trait_change, Button, \
    Interface, implements, Property, cached_property

from ibvpy.api import \
    TStepper as TS, TLoop, TLine, \
    IBVModel, DOTSEval, \
    RTraceGraph, RTraceDomainListField, \
    BCDof, BCDofGroup, BCSlice

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.fets.fets3D.fets3D8h20u import \
    FETS3D8H20U
from ibvpy.fets.fets3D.fets3D8h27u import \
    FETS3D8H27U
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mats.mats2D5.mats2D5_cmdm.mats2D5_cmdm import \
    MATS2D5MicroplaneDamage

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.matsXD.matsXD_cmdm.matsXD_cmdm_phi_fn import \
    PhiFnGeneralExtended, \
    PhiFnGeneral, PhiFnStrainHardening, PhiFnStrainHardeningLinear

from mathkit.geo.geo_ndgrid import \
    GeoNDGrid

from mathkit.mfn.mfn_ndgrid.mfn_ndgrid import \
    MFnNDGrid, GridPoint

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

from numpy import \
    sin, cos, c_, arange, hstack, array, max, frompyfunc, linspace

from time import time
from os.path import join

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from simiter.sim_pstudy import \
    SimPStudy, SimOut, ISimModel

from promod.exdb.ex_run_view import \
    ExRunView

from promod.matdb.trc.ccs_unit_cell import \
    CCSUnitCell, DamageFunctionEntry

from promod.simdb import \
    SimDB

from promod.simdb.simdb_class import \
    SimDBClass, SimDBClassExt

simdb = SimDB()

from pickle import dump, load


class SimQuadPlate( IBVModel ):
    '''Plate test prepared for parametric study.
    '''

    input_change = Event
    @on_trait_change( '+input,ccs_unit_cell.input_change' )
    def _set_input_change( self ):
        self.input_change = True

    implements( ISimModel )

    #-----------------
    # discretization:
    #-----------------
    #
    # discretization in x,y-direction:
    shape_xy = Int( 8, input = True,
                      ps_levels = ( 8, 12, 3 ) )

    # discretization in z-direction:
    shape_z = Int( 2, input = True,
                      ps_levels = ( 1, 3, 3 ) )

    #-----------------
    # geometry:
    #-----------------
    #
    # edge length of the quadratic plate (entire length without symmetry)
    length = Float( 1.25, input = True
                       )

    thickness = Float( 0.03, input = True
                       )

    # vary the failure strain in PhiFnGeneralExtended:
    factor_eps_fail = Float( 1.0, input = True,
                             ps_levels = ( 1.0, 1.2, 3 ) )

    #-----------------
    # composite cross section unit cell:
    #-----------------
    #
    ccs_unit_cell_key = Enum( CCSUnitCell.db.keys(),
                              simdb = True, input = True,
                              auto_set = False, enter_set = True )

    ccs_unit_cell_ref = Property( Instance( SimDBClass ),
                                  depends_on = 'ccs_unit_cell_key' )
    @cached_property
    def _get_ccs_unit_cell_ref( self ):
        return CCSUnitCell.db[ self.ccs_unit_cell_key ]

    #-----------------
    # damage function:
    #-----------------
    #
    material_model = Str( input = True )
    def _material_model_default( self ):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup  
        return self.ccs_unit_cell_ref.damage_function_list[0].material_model

    calibration_test = Str( input = True
                            )
    def _calibration_test_default( self ):
        # return the material model key of the first DamageFunctionEntry
        # This is necessary to avoid an ValueError at setup  
        return self.ccs_unit_cell_ref.damage_function_list[0].calibration_test

    damage_function = Property( Instance( MFnLineArray ),
                                depends_on = 'input_change' )
    @cached_property
    def _get_damage_function( self ):
        return self.ccs_unit_cell_ref.get_param( self.material_model, self.calibration_test )

    #-----------------
    # phi function extended:
    #-----------------
    #
    phi_fn = Property( Instance( PhiFnGeneralExtended ),
                       depends_on = 'input_change,+ps_levels' )
    @cached_property
    def _get_phi_fn( self ):
        return PhiFnGeneralExtended( mfn = self.damage_function,
                                     factor_eps_fail = self.factor_eps_fail )

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int( 28, #input = True
                )

    # composite E-modulus 
    #
    E_c = Property( Float, depends_on = 'input_change' )
    @cached_property
    def _get_E_c( self ):
        return self.ccs_unit_cell_ref.get_E_c_time( self.age )

    # Poisson's ratio 
    #
    nu = Property( Float, depends_on = 'input_change' )
    @cached_property
    def _get_nu( self ):
        return self.ccs_unit_cell_ref.nu

    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    mats_eval = Property( Instance( MATS2D5MicroplaneDamage ),
                          depends_on = 'input_change' )
    @cached_property
    def _get_mats_eval( self ):
        mats_eval = MATS2D5MicroplaneDamage( 
                                E = self.E_c,
                                nu = self.nu,
                                n_mp = 30,
                                symmetrization = 'sum-type',
                                model_version = 'compliance',
                                phi_fn = self.phi_fn )

        return mats_eval

    #-----------------
    # fets:
    #-----------------

    # use quadratic serendipity elements
    #
    fets = Property( Instance( FETSEval ),
                     depends_on = 'input_change' )
    @cached_property
    def _get_fets( self ):
        return FETS2D58H20U( mats_eval = self.mats_eval )

    def peval( self ):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        self.f_w_diagram_center.refresh()
        F_max = max( self.f_w_diagram_center.trace.ydata )

        u_center_top_z = U[ self.center_top_dofs ][0, 0, 2]
        return array( [ u_center_top_z, F_max ],
                        dtype = 'float_' )

    tloop = Property( depends_on = 'input_change' )
    @cached_property
    def _get_tloop( self ):

        fets_eval = self.fets
        #fets_eval.mats_eval.nu = self.nu

        # only a quarter of the plate is simulated due to symmetry:
        domain = FEGrid( coord_max = ( self.length / 2, self.length / 2, self.thickness ),
                         shape = ( self.shape_xy, self.shape_xy, self.shape_z ),
                         fets_eval = fets_eval )

        self.center_top_dofs = domain[-1, -1, -1, -1, -1, -1].dofs
        self.edge_top_dofs = domain[ 0, -1, -1, 0, -1, -1].dofs

        # if shape_xy is an even number:
        if self.shape_xy % 2 == 0:
            x_idx_center_edge = ( self.shape_xy / 2 ) - 1
            self.center_edge_top_dofs = domain[x_idx_center_edge, -1, -1, -1, -1, -1].dofs
        # if shape_xy is an odd number:
        else:
            # get the midside node of the center-edge-element
            x_idx_center_edge = ( self.shape_xy - 1 ) / 2
            # valid only for quadratic elements
            self.center_edge_top_dofs = domain[x_idx_center_edge, -1, -1, 1, -1, -1].dofs

        #--------------------------------------------------------------
        # boundary conditions for the symmetry and the single support
        #--------------------------------------------------------------
        bc_symplane_yz = BCSlice( var = 'u', value = 0., dims = [0], slice = domain[-1, :, :, -1, :, :] )
        bc_symplane_xz = BCSlice( var = 'u', value = 0., dims = [1], slice = domain[:, -1, :, :, -1, :] )
        bc_support_000 = BCSlice( var = 'u', value = 0., dims = [2], slice = domain[0, 0, 0, 0, 0, 0] )

        #--------------------------------------------------------------
        # loading
        #--------------------------------------------------------------
        # w_max = center displacement:
        w_max = -0.08 # [m]

        #--------------------------------------------------------------
        # var 1: nodal displacement applied at top at the center of the plate
        #--------------------------------------------------------------
        bc_center_w = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-1, -1, -1, -1, -1, -1] )

        #--------------------------------------------------------------
        # var 2: apply nodal displacements at top nodes of center element 
        #--------------------------------------------------------------
        bc_center_w_elem = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-1, -1, -1, :, :, -1] )

        #--------------------------------------------------------------
        # var 3: displacement applied along the edge of the 2nd element row (from the middle) for 
        # shape_xy = 10 an an edge length of L_quarter = 0.625m the distance from the center equals 0.125m
        #--------------------------------------------------------------
        bc_center_w_xline = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-2:, -2, -1, :, 0, -1] )
        bc_center_w_yline = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-2, -2:, -1, 0, :, -1] )


        #--------------------------------------------------------------
        # ts 
        #--------------------------------------------------------------

        # force-displacement-diagramm 
        # 
        center_dof = self.center_top_dofs[0, 0, 2]
        self.f_w_diagram_center = RTraceGraph( name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_dof,
                                       var_y = 'F_int', idx_y = 2,
                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       transform_y = '4 * 1000 * y' )

        center_edge_dof = self.center_edge_top_dofs[0, 0, 2]
        self.f_w_diagram_center_edge = RTraceGraph( name = 'displacement (center/edge) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_edge_dof,
                                       var_y = 'F_int', idx_y = 2,
                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       transform_y = '4 * 1000 * y' )

        edge_dof = self.edge_top_dofs[0, 0, 2]
        self.f_w_diagram_edge = RTraceGraph( name = 'displacement (edge) - reaction 2',
                                       var_x = 'U_k'  , idx_x = edge_dof,
                                       var_y = 'F_int', idx_y = 2,
                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       transform_y = '4 * 1000 * y' )

        ts = TS( 
                sdomain = domain,
                bcond_list = [bc_symplane_yz, bc_symplane_xz, bc_support_000,
                               # var 1:
                               bc_center_w,
#                               # var 2:
#                               bc_center_w_elem,
#                               # var 3:
#                               bc_center_w_xline, bc_center_w_yline
                              ],
                rtrace_list = [
                             self.f_w_diagram_center,
                             self.f_w_diagram_center_edge,
                             self.f_w_diagram_edge,
                             RTraceDomainListField( name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True ),

#                             RTraceDomainListField(name = 'Stress' ,
#                                            var = 'sig_app', idx = 0, warp = True, 
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'Strain' ,
#                                        var = 'eps_app', idx = 0, warp = True, 
#                                        record_on = 'update'),
                             RTraceDomainListField( name = 'Damage' ,
                                        var = 'omega_mtx', idx = 0, warp = True,
                                        record_on = 'update' ),

#                             RTraceDomainListField(name = 'IStress' ,
#                                            position = 'int_pnts',
#                                            var = 'sig_app', idx = 0, 
#                                            record_on = 'update'),
#                             RTraceDomainListField(name = 'IStrain' ,
#                                            position = 'int_pnts',
#                                            var = 'eps_app', idx = 0, 
#                                            record_on = 'update'),
                              ]
                )

        # Add the time-loop control
        tloop = TLoop( tstepper = ts,

#                       # allow only a low tolerance 
#                       #
#                       KMAX = 50, 
#                       tolerance = 5e-4, 

                       # allow a high tolerance 
                       #
                       KMAX = 20,
                       tolerance = 0.001,

                       RESETMAX = 0,
                       debug = False,
                       tline = TLine( min = 0.0, step = 0.05, max = 1.0 ) )

        return tloop

    def get_sim_outputs( self ):
        '''
        Specifies the results and their order returned by the model
        evaluation.
        '''
        return [ SimOut( name = 'u_center_top_z', unit = 'm' ),
                 SimOut( name = 'F_max', unit = 'kN' ) ]


if __name__ == '__main__':

    sim_model = SimQuadPlate( ccs_unit_cell_key = 'FIL-10-09_2D-02-06a_0.00273_90_0',
                              calibration_test = 'TT11-10a-average',
                              age = 28 )

    do = 'validation'
#    do = 'show_last_results'

    if do == 'ui':
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp( ibv_resource = sim_model )
        app.main()

    if do == 'pstudy':
        sim_ps = SimPStudy( sim_model = sim_model )
        sim_ps.configure_traits()

    if do == 'validation':

        from promod.exdb.ex_run import ExRun
        import pylab as p

        # PT-9a
#        path = join( simdb.exdata_dir, 'plate_tests', 'PT-9a' )
#        tests = [ 'PT01-9a.DAT', 'PT02-9a.DAT' , 'PT09-9a.DAT' ]

        # PT-10a
        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]

        for t in tests:
            ex_path = join( path, t )
            ex_run = ExRun( ex_path )
            ex_run.ex_type._plot_smoothed_force_deflection_center( p )

        sim_model.tloop.eval()

        sim_model.f_w_diagram_center.refresh()
        pickle_path = 'pickle_files'
        file_name = 'f_w_diagram_c_average_step0-05_KMAX20_tol0-001_nelems12-12-4.pickle'
        pickle_file_path = join( pickle_path, file_name )
        file = open( pickle_file_path, 'w' )
        dump( sim_model.f_w_diagram_center.trace, file )
        file.close()

        sim_model.f_w_diagram_center.trace.mpl_plot( p, color = 'red' )

        p.show()

    if do == 'show_last_results':
        from promod.exdb.ex_run import ExRun
        import pylab as p

        pickle_path = 'pickle_files'

        ### shape ###:

        # plate elem(12,12,4); w=0.08; step = 0.05; KMAX = 20, tol = 0.001
        file_name = 'f_w_diagram_c_average_step0-05_KMAX20_tol0-001_nelems12-12-4.pickle'
        pickle_file_path = join( pickle_path, file_name )
        file = open( pickle_file_path, 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'blue' )

#        # plate elem(8,8,4); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems8-8-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

        # plate elem(10,10,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_nelems10-10-2.pickle'
        pickle_file_path = join( pickle_path, file_name )
        file = open( pickle_file_path, 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'black' )

        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
        pickle_file_path = join( pickle_path, file_name )
        file = open( pickle_file_path, 'r' )
        trace = load( file )
        p.plot( trace.xdata, trace.ydata, color = 'red' )


        ### eps_factor ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001, factor_eps_fail = 1.4
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_epsfactor_1-4.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001, factor_eps_fail = 1.2
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001_epsfactor_1-2.pickle'
#        pickle_file_path = join( pickle_path, file_name )
        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        ### average - V1 (c) ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_V1_step0-1_KMAX20_tol0-001.pickle'
        pickle_file_path = join( pickle_path, file_name )
        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )


        ### c ce e - average ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_ce_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_e_average_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )


        ### c ce e - V1 ###:

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_c_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )

#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_ce_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_e_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )


        ### step ###:

#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-02_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'red' )
#
#        # plate elem(8,8,2); w=0.08; step = 0.1; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-1_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'black' )


        ### tolerance ###:

#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 50, tol = 0.0005
#        file_name = 'f_w_diagram_V1_step0-02_KMAX50_tol0-0005.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'green' )
#        
#        # plate elem(8,8,2); w=0.08; step = 0.02; KMAX = 20, tol = 0.001
#        file_name = 'f_w_diagram_V1_step0-02_KMAX20_tol0-001.pickle'
#        pickle_file_path = join( pickle_path, file_name )
#        file = open( pickle_file_path, 'r' )
#        trace = load( file )
#        p.plot( trace.xdata, trace.ydata, color = 'yellow' )


        path = join( simdb.exdata_dir, 'plate_tests', 'PT-10a' )
        tests = [ 'PT10-10a.DAT', 'PT11-10a.DAT' , 'PT12-10a.DAT' ]
#        tests = [ 'PT10-10a.DAT' ]

        for t in tests:
            ex_path = join( path, t )
            ex_run = ExRun( ex_path )
            ex_run.ex_type._plot_smoothed_force_deflection_center( p )
#            ex_run.ex_type._plot_force_edge_deflection( p )
#            ex_run.ex_type._plot_force_center_edge_deflection( p )

        p.show()
