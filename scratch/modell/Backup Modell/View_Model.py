
from enthought.traits.api import HasTraits, Int, List, Float, Str, Property, Range, Array, Bool, Enum
from etsproxy.traits.ui.api import ModelView, View, Item, Label, Group, HGroup
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
from etsproxy.traits.api import Instance, cached_property, Button
import numpy as np
from spirrid.rv import RV
from scipy.optimize import brentq, fminbound
from scipy.integrate import cumtrapz
from mathkit.mfn.mfn_line.mfn_line import MFnLineArray
import time
from scipy.optimize import minimize
from scipy.stats import exponweib
from hom_CB_elastic_mtrx import CompositeCrackBridge
from hom_CB_elastic_mtrx_view import CompositeCrackBridgeView
from matplotlib import pyplot as plt
from reinforcement import ContinuousFibers, ShortFibers
from stats.pdistrib.weibull_fibers_composite_distr import WeibullFibers
from traitsui.api import CheckListEditor
from spirrid import make_ogrid as orthogonalize
from mayavi import mlab
from mayavi.api import Engine
from numpy import array
import os

# from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_py_loop import CompositeCrackBridgeLoop
class Plot_app( HasTraits ):
    # MATERIAL PARAMETERS
    # reinf1
    x_name3D = Enum( 'V_f1', 'V_f2', 'V_f3', 'V_f4', 'r1', 'r2', 'r3', 'r4', 'E_f1', 'E_f2', 'E_f3', 'E_f4', 'lf3', 'lf4', 'Ll', 'Lr' )
    y_name3D = Enum( 'V_f3', 'V_f1', 'V_f2', 'V_f4', 'r1', 'r2', 'r3', 'r4', 'E_f1', 'E_f2', 'E_f3', 'E_f4', 'lf3', 'lf4', 'Ll', 'Lr' )
    Ll = Float( 100. )
    Lr = Float( 100. )
    discr_amin = Int( 50 )
    ###############
    standard_configurations_1 = Enum( 'custom', 'Carbon', 'PolyEthylen', 'AR_Glass' )
    standard_configurations_2 = Enum( 'custom', 'Carbon', 'PolyEthylen', 'AR_Glass' )
    standard_configurations_3 = Enum( 'custom', 'Steel', 'PolyEthylen', 'AR_Glass' )
    standard_configurations_4 = Enum( 'custom', 'Steel', 'PolyEthylen', 'AR_Glass' )
    load1 = Bool( False )
    load2 = Bool( False )
    load3 = Bool( False )
    load4 = Bool( False )
    
    def config_loader1( self ):
        if self.standard_configurations_1 == 'custom':
            pass
        else:
            exec '{}{}{}'.format( 'self.', self.standard_configurations_1, '_c1()' )

    def config_loader2( self ):
        if self.standard_configurations_2 == 'custom':
            pass
        else:
            exec '{}{}{}'.format( 'self.', self.standard_configurations_2, '_c2()' )

    def config_loader3( self ):
        if self.standard_configurations_3 == 'custom':
            pass
        else:
            exec '{}{}{}'.format( 'self.', self.standard_configurations_3, '_c3()' )

    def config_loader4( self ):
        if self.standard_configurations_4 == 'custom':
            pass
        else:
            exec '{}{}{}'.format( 'self.', self.standard_configurations_4, '_c4()' ) 
    
    
    
    
    def Carbon_c1( self ):
        self.r1 = 3.45e-3
        self.tau1 = ['weibull_min', 0.006, .23, .03 ]
        self.V_f1 = 0.01
        self.E_f1 = 240e3
        self.xi_shape1 = 5.
        self.xi_sv01 = 0.0026 
        self.n_int1 = 50 
        
    def Carbon_c2( self ):
        self.r2 = 3.45e-3
        self.tau2 = ['weibull_min', 0.006, .23, .03 ]
        self.V_f2 = 0.01
        self.E_f2 = 240e3
        self.xi_shape2 = 5.
        self.xi_sv02 = 0.0026 
        self.n_int2 = 50 
        
    def PolyEthylen_c1( self ):
        pass
    def PolyEthylen_c2( self ):
        pass
    def AR_Glass_c1( self ):
        pass
    def AR_Glass_c2( self ):
        pass
    
    def AR_Glass_c3( self ):
        pass
    
    def AR_Glass_c4( self ):
        pass
    
    def PolyEthylen_c3( self ):
        self.r3 = 3.8e-3
        self.tau3 = 1.02
        self.lf3 = 12.7 
        self.snub3 = 0.7
        self.V_f3 = 1
        self.E_f3 = 120e3
        self.xi3 = 100. 
        self.n_int3 = 50 
    
    def PolyEthylen_c4( self ):
        self.r4 = 3.8e-3
        self.tau4 = 1.02
        self.lf4 = 12.7 
        self.snub4 = 0.7
        self.V_f4 = 1
        self.E_f4 = 120e3
        self.xi4 = 100. 
        self.n_int4 = 50 
    
    def Steel_c3( self ):
        self.r3 = 0.3 
        self.tau3 = 1.76 
        self.lf3 = 17. 
        self.snub3 = .03 
        self.V_f3 = 0.02 
        self.E_f3 = 200e3 
        self.xi3 = 100. 
        self.n_int3 = 50 
    
    def Steel_c4( self ):
        self.r4 = 0.3 
        self.tau4 = 1.76 
        self.lf4 = 17. 
        self.snub4 = .03 
        self.V_f4 = 0.02 
        self.E_f4 = 200e3 
        self.xi4 = 100. 
        self.n_int4 = 50 
    
        
    ###############
    # reinf1
    r1 = Float( 0.00345 )
    tau1 = List
    def _tau1_default( self ):
        return ['weibull_min', 0.006, .23, .03 ]
    
    tau1rv = Property( depends_on = 'tau1' )
    def _get_tau1rv( self ):
         if len( self.tau1 ) == 3:
                    return RV( self.tau1[0], loc = np.float( self.tau1[1] ), scale = np.float( self.tau1[2] ) )
         else:
                    return RV( self.tau1[0], loc = np.float( self.tau1[1] ), shape = np.float( self.tau1[2] ) , scale = np.float( self.tau1[3] ) )
    
    V_f1 = Float( 0.03 )
    E_f1 = Float( 240e3 )
    xi_shape1 = Float( 4. )
    xi_sv01 = Float( 0.0026 )
    n_int1 = Int( 50 )
    ###############
    # reinf2
    r2 = Float( 0.00345 )
    tau2 = List
    def _tau2_default( self ): 
        return ['weibull_min', 0.006, .23, .03 ]
    
    tau2rv = Property( depends_on = 'tau1' )
    def _get_tau2rv( self ):
         if len( self.tau1 ) == 3:
                    return RV( self.tau2[0], loc = np.float( self.tau2[1] ), scale = np.float( self.tau2[2] ) )
         else:
                    return RV( self.tau2[0], loc = np.float( self.tau2[1] ), shape = np.float( self.tau2[2] ) , scale = np.float( self.tau2[3] ) )
    
    
    V_f2 = Float( 0.03 )
    E_f2 = Float( 180e3 )
    xi_shape2 = Float( 4. )
    xi_sv02 = Float( 0.0025 )
    n_int2 = Int( 50 )
###########################
    # reinf3
    r3 = Float( 0.3 )
    tau3 = Float( 1.76 )
    lf3 = Float( 17. )
    snub3 = Float( .03 )
    phi3 = List
    def _phi3_default( self ):
        return ['sin2x', 0., 1.]
    V_f3 = Float( 0.02 )
    E_f3 = Float( 200e3 )
    xi3 = Float( 100. )
    n_int3 = Int( 50 )
    # reinf4
    r4 = Float( 0.3 )
    tau4 = Float( 1.76 )
    lf4 = Float( 17. )
    snub4 = Float( .03 )
    phi4 = List
    def _phi4_default( self ):
        return ['sin2x', 0., 1.]
    V_f4 = Float( 0.02 )
    E_f4 = Float( 180e3 )
    xi4 = Float( 100. )
    n_int4 = Int( 50 )

    model = Instance( CompositeCrackBridge ) 
    def _model_default( self ):
        return CompositeCrackBridge( E_m = 25e3,
                                 reinforcement_lst = [self.reinf1, self.reinf3],
                                 Ll = 1000.,
                                 Lr = 1000.,
                                 discr_amin = 70 )
    ccb_view = Instance( CompositeCrackBridgeView )
    def _ccb_view_default( self ):
        return CompositeCrackBridgeView( model = self.model )
    
    
    
    def set_model( self ):
        # TODO to be integrated in model- still only works with cont and sf fibers
        self.reinf2_bool
        reinf_ls = []
        if self.reinf1_bool:
            reinf_ls.append( self.reinf1 )
        if self.reinf2_bool:
            reinf_ls.append( self.reinf2 )
        if self.reinf3_bool:
            reinf_ls.append( self.reinf3 )
        if self.reinf4_bool:
            reinf_ls.append( self.reinf4 )
        self.ccb_view.model.reinforcement_lst = reinf_ls
        self.ccb_view.model.Ll = self.Ll
        self.ccb_view.model.Lr = self.Lr
        self.ccb_view.model.discr_amin = self.discr_amin  
    
    
    reinf1_bool = Bool( True )
    reinf2_bool = Bool( False )
    reinf3_bool = Bool( True )
    reinf4_bool = Bool( False )
    
    reinf1 = Property( depends_on = 'material_ok' )
    def _get_reinf1( self ):
        return ContinuousFibers( r = self.r1,
                                tau = self.tau1rv,
                                V_f = self.V_f1,
                                E_f = self.E_f1,
                                xi = WeibullFibers( shape = self.xi_shape1, sV0 = self.xi_sv01 ),
                                n_int = self.n_int1,
                                label = 'carbon' )
        
    reinf2 = Property( depends_on = 'material_ok' )
    def _get_reinf2( self ):
        return ContinuousFibers( r = self.r2,
                                tau = RV( self.tau2[0], loc = np.float( self.tau2[1] ), scale = np.float( self.tau2[2] ) ),
                                V_f = self.V_f2,
                                E_f = self.E_f2,
                                xi = WeibullFibers( shape = self.xi_shape2, sV0 = self.xi_sv02 ),
                                n_int = self.n_int2,
                                label = 'carbon' )

    
    reinf3 = Property( depends_on = 'material_ok' )
    
    def _get_reinf3( self ):
        return ShortFibers( r = self.r3,
                                tau = self.tau3,
                                lf = self.lf3,
                                snub = self.snub3,
                                phi = RV( ( self.phi3[0] ), loc = np.float( self.phi3[1] ), scale = np.float( self.phi3[2] ) ),
                                V_f = self.V_f3,
                                E_f = self.E_f3,
                                xi = self.xi3 ,
                                n_int = self.n_int3,
                                label = 'Short Fibers' )
    
    reinf4 = Property( depends_on = 'material_ok' )
    def _get_reinf4( self ):
        return ShortFibers( r = self.r4,
                                tau = self.tau4,
                                lf = self.lf4,
                                snub = self.snub4,
                                phi = RV( ( self.phi4[0] ), loc = np.float( self.phi4[1] ), scale = np.float( self.phi4[2] ) ),
                                V_f = self.V_f4,
                                E_f = self.E_f4,
                                xi = self.xi4 ,
                                n_int = self.n_int4,
                                label = 'Short Fibers' )

    
    
    w_profile = Float( 1. )
    def profile( self ):
        w = self.w_profile
        self.ccb_view.model.w = w
        plt.plot( self.ccb_view.x_arr, self.ccb_view.epsm_arr, color = 'blue', lw = 2 )  # , label='w_eval=' + str(ccb_view.w_evaluated) + ' w_ctrl=' + str(ccb_view.model.w))
        plt.plot( self.ccb_view.x_arr, self.ccb_view.mu_epsf_arr, color = 'red', lw = 2 )
        plt.xlabel( 'position [mm]' )
        plt.ylabel( 'strain' )
        plt.show( 'block' )
   
    choice3D = Enum( 'sigma', 'w' )

    renamed_choice3D = Property( depends_on = 'choice3D' )
    @cached_property
    def _get_renamed_choice3D( self ):
        return ['sigma', 'w' ].index( self.choice3D )
    
    ls_xfrom_3D = Float( 0 )
    ls_xtil_3D = Float( 0.1 )
    ls_xdp_3D = Int( 10 )
    ls_yfrom_3D = Float( 0 )
    ls_ytil_3D = Float( 0.1 )
    ls_ydp_3D = Int( 10 )
    plot3d_saveto = Str( 'name_it' )
    save3D = Bool( False )
        
        
    def save_plot3D( self ):
        # createFolderNamed 
        try:
            os.chdir( 'DATA' )
            os.makedirs( self.plot3d_saveto )
            os.chdir( self.plot3d_saveto )
            x_axis = self.compute3D_plot[0]
            y_axis = self.compute3D_plot[1]
            z_axis = self.compute3D_plot[2] 
            np.save( 'x', x_axis )
            np.save( 'y', y_axis )
            np.save( 'z', z_axis )
            os.chdir( os.pardir )
            os.chdir( os.pardir )
        except OSError:
            os.chdir( 'DATA1' )
            os.makedirs( self.plot3d_saveto )
            os.chdir( self.plot3d_saveto )
            x_axis = self.compute3D_plot[0]
            y_axis = self.compute3D_plot[1]
            z_axis = self.compute3D_plot[2] 
            np.save( 'x', x_axis )
            np.save( 'y', y_axis )
            np.save( 'z', z_axis )
            os.chdir( os.pardir )
            os.chdir( os.pardir )
            # let exception propagate if we just can't
            # cd into the specified directory

    compute3D_plot = Property()
    @cached_property
    def _get_compute3D_plot( self ):
        # lib_list = self.library_read( 'plot3D' )
        # if lib_list[0] == True:
        #    return lib_list[1], lib_list[2], lib_list[3]
        # else: 
            choice = self.renamed_choice3D
            z_axis_y = []
            z_axis = []
            x_string_exec = '{}{}{}'.format( 'self.', self.x_name3D, '=x' )
            y_string_exec = '{}{}{}'.format( 'self.', self.y_name3D, '=y' )
            x_ls = np.linspace( self.ls_xfrom_3D, self.ls_xtil_3D, self.ls_xdp_3D )
            y_ls = np.linspace( self.ls_yfrom_3D, self.ls_ytil_3D, self.ls_ydp_3D )
            # method = [self.ccb_view.sigma_c_max[0], self.ccb_view.sigma_c_max[1]]
            for i, x in enumerate( x_ls ):
                exec x_string_exec
                for y in y_ls:
                    exec y_string_exec
                    self._material_ok_fired()
                    self.set_model()
                    z_i = self.ccb_view.sigma_c_max[choice]
                    z_axis_y.append( z_i )
                z_axis.append( z_axis_y )
                z_axis_y = []
            e_arr = orthogonalize( [x_ls, y_ls] )
            x_axis = e_arr[0]
            y_axis = e_arr[1]
            return x_axis, y_axis , np.array( z_axis )
    
    
    load_plot_bool = Bool( False )
    Button_load_plot3D = Button( 'launch' )
    def _Button_load_plot3D_fired( self ):
        self.load_3D_plot()
    
    
    autowarp_bool = Bool( True )
    z_scale = Float( 2. )
    y_scale = Float( 1. )
    x_scale = Float( 1. )
    
    def load_3D_plot( self ):
        config = 0
        os.chdir( 'DATA' )
        # os.chdir( self.load_folder3D )
        os.chdir( self.dropdownls3D )
        x = np.load( 'x.npy' )
        y = np.load( 'y.npy' )
        z = np.load( 'z.npy' )
        print 'min_x', np.min( x )
        print 'max_x', np.max( x )
        print 'min_y', np.min( y )
        print 'max_y', np.max( y )
        print 'min_z', np.min( z )
        print 'max_z', np.max( z )
        # print x_axis, y_axis, z_axis
        engine = Engine()
        engine.start()
        if len( engine.scenes ) == 0:
            engine.new_scene  # print os.curdir
            # os.chdir( '.' )()
        if self.autowarp_bool:
            x = x / x[-1]
            y = y / y[0][-1]
            z = z / z[-1] * self.z_scale
        
        # print x, y, z
        mlab.surf( x, y , z, representation = 'wireframe', line_width = 10 )  # ,warp_scale = "auto"
        if 'config.py' in os.listdir( os.curdir ):
            execfile( 'config.py' )
            os.chdir( os.pardir )
            os.chdir( os.pardir )
        else:
            os.chdir( os.pardir )
            os.chdir( os.pardir )
            surface = engine.scenes[0].children[0].children[0].children[0].children[0].children[0]
            surface.actor.mapper.scalar_range = np.array( [ 6.97602671, 8.8533387 ] )
            surface.actor.mapper.scalar_visibility = False
            scene = engine.scenes[0]
            scene.scene.background = ( 1.0, 1.0, 1.0 )
            surface.actor.property.specular_color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.diffuse_color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.ambient_color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.line_width = 1.
            scene.scene.isometric_view()
        mlab.outline()
        outline = engine.scenes[0].children[0].children[0].children[0].children[0].children[1]
        outline.actor.property.specular_color = ( 0.0, 0.0, 0.0 )
        outline.actor.property.diffuse_color = ( 0.0, 0.0, 0.0 )
        outline.actor.property.ambient_color = ( 0.0, 0.0, 0.0 )
        outline.actor.property.color = ( 0.0, 0.0, 0.0 )

        mlab.show()
            
        
    def plot_3D( self ):
        x = self.compute3D_plot[0]
        y = self.compute3D_plot[1]
        z = self.compute3D_plot[2]
        # print x_axis, y_axis, z_axis
        if self.autowarp_bool:
            x = x / x[-1]
            y = y / y[-1]
            z = z / z[-1] * self.z_scale
        mlab.surf( x, y , z, representation = 'wireframe' )
        engine = Engine()
        engine.start()
        if len( engine.scenes ) == 75.5:
            engine.new_scene()
            surface = engine.scenes[0].children[0].children[0].children[0].children[0].children[0]
            surface.actor.mapper.scalar_range = np.array( [ 6.97602671, 8.8533387 ] )
            surface.actor.mapper.scalar_visibility = False
            scene = engine.scenes[0]
            scene.scene.background = ( 1.0, 1.0, 1.0 )
            surface.actor.property.specular_color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.diffuse_color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.ambient_color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.color = ( 0.0, 0.0, 0.0 )
            surface.actor.property.line_width = 1.
            scene.scene.isometric_view()
            mlab.xlabel( self.x_name3D )
            mlab.ylabel( self.y_name3D )
        mlab.outline()
        mlab.show()
    
    launch_profile = Button( 'launch' )
    click_counter = Int
    def _launch_profile_fired( self ):
        self.set_model()
        self.profile()
        
    material_ok = Button( 'OK' )
    def _material_ok_fired( self ):
        if self.load1:
            self.config_loader1()
        if self.load2:
            self.config_loader2()
        if self.load3:
            self.config_loader3()
        if self.load4:
            self.config_loader4()
        pass
    
    plot3D_launch = Button( 'launch' )
    def _plot3D_launch_fired( self ):
        self.load1 = False
        self.load2 = False
        self.load3 = False
        self.load4 = False
        self.compute3D_plot
        if self.save3D:
            self.save_plot3D()
        self.plot_3D()

        
    dropdownls3D = Enum( os.listdir( 'DATA' ) )
    '''
    def _dropdownls3D_default( self ):
        print os.listdir( os.curdir )
        os.chdir( 'DATA' )
        all_studies = os.listdir( os.curdir )
        os.chdir( os.pardir )
        return  enumerate( all_studies )
    def _get_dropdownls3D( self ):
        print os.listdir( os.curdir )
        os.chdir( 'DATA' )
        all_studies = os.listdir( os.curdir )
        os.chdir( os.pardir )
        return  enumerate( all_studies )
    '''       
    material_cont_ribbon = Group( HGroup( Item( name = 'reinf1_bool', label = 'reinf1' ), Item( 'standard_configurations_1', label = 'preconfiguration' ), Item( 'load1', label = 'load config' ) ),
                            Item( name = 'r1', label = 'r' ),
                             Item( name = 'V_f1', label = 'Vf' ),
                             Item( name = 'tau1' , label = 'tau' ),
                             Item( name = 'E_f1', label = 'Ef' ),
                             Item( name = 'xi_shape1', label = 'xi_shape' ),
                             Item( name = 'xi_sv01', label = 'xi_sv0' ),
                             Item( name = 'n_int1', label = 'integration_points' ),
                             Item( '_' ),
                             HGroup( Item( name = 'reinf2_bool', label = 'reinf2' ), Item( 'standard_configurations_2', label = 'preconfiguration' ), Item( 'load2', label = 'load config' ) ),
                             Item( name = 'r2', label = 'r' ),
                             Item( name = 'V_f2', label = 'Vf' ),
                             Item( name = 'tau2' , label = 'tau' ),
                             Item( name = 'E_f2', label = 'Ef' ),
                             Item( name = 'xi_shape2', label = 'xi_shape' ),
                             Item( name = 'xi_sv02', label = 'xi_sv0' ),
                             Item( name = 'n_int2', label = 'integration_points' ),
                               Item( '_' ),
                               'material_ok', label = 'Cont Fiber Reinf' )
    material_sf_ribbon = Group( HGroup( Item( name = 'reinf3_bool', label = 'reinf3' ), Item( 'standard_configurations_3', label = 'preconfiguration' ) , Item( 'load3', label = 'load config' ) ),
                            Item( name = 'r3', label = 'r' ),
                              Item( name = 'tau3', label = 'tau' ),
                              Item( name = 'phi3' , label = 'phi' ) ,
                               Item( name = 'V_f3', label = 'Vf' ),
                               Item( name = 'snub3' , label = 'snub' ),
                               Item( name = 'E_f3', label = 'Ef' ),
                               Item( name = 'xi3', label = 'xi' ),
                               Item( name = 'n_int3', label = 'integration_points' ),
                             HGroup( Item( name = 'reinf4_bool', label = 'reinf4' ), Item( 'standard_configurations_4', label = 'preconfiguration' ), Item( 'load4', label = 'load config' ) ),
                              Item( name = 'r4', label = 'r' ),
                              Item( name = 'tau4' , label = 'tau' ),
                              Item( name = 'phi4' , label = 'phi' ) ,
                               Item( name = 'V_f4', label = 'Vf' ),
                               Item( name = 'snub4' , label = 'snub' ),
                               Item( name = 'E_f4', label = 'Ef' ),
                               Item( name = 'xi4' , label = 'xi' ),
                               Item( name = 'n_int4', label = 'integration_points' ),
                               Item( '_' ),
                               'material_ok', label = 'Short Fiber Reinf' )
    profile_ribbon = Group( Item( name = 'w_profile' , label = 'w in [mm]' ), Item( name = 'launch_profile', label = 'launch' ), label = 'profile' )
    
    plot3D_ribbon = Group( Item( name = 'choice3D', label = 'Choice of : (0) sigma_max,(1) w' ),
                           Item( '_' ),
                           HGroup( Item( name = 'autowarp_bool', label = 'normed axes' ),
                                  Item( name = 'z_scale', label = 'scale factor of z' ),
                                  Item( name = 'y_scale', label = 'scale factor of y' ),
                                  Item( name = 'x_scale', label = 'scale factor of x' ) ),
                           Item( '_' ),
                           HGroup( Item( name = 'x_name3D', label = 'x-Axis variable' ),
                               Item( name = 'ls_xfrom_3D' , label = 'x-linspace begins at' ),
                               Item( name = 'ls_xtil_3D', label = 'goes until' ),
                                Item( name = 'ls_xdp_3D', label = 'DataPoints ' ) ),
                            Item( '_' ),
                            HGroup( Item( name = 'y_name3D' , label = ' y-axis variable' ),
                             Item( name = 'ls_yfrom_3D', label = 'y-linspace begins at' ),
                              Item( name = 'ls_ytil_3D' , label = 'goes until' ),
                               Item( name = 'ls_ydp_3D' , label = 'DataPoints ' ) ),
                          Item( '_' ),
                            Item( name = 'save3D', label = 'Save File ?' ),
                             Item( name = 'plot3d_saveto', label = 'foldername' ),
                             Item( name = 'plot3D_launch', label = 'start when ready' ) ,
                             Item( name = '_' ),
                             Item( 'dropdownls3D' ),
                            Item( name = 'Button_load_plot3D', label = 'load from file' ),
                                   
                             label = 'Plot3D' ) 
    # ,Item( 'y_name3D)',Item( 'z_name3D' ),Item( 'ls_yfrom_3D' ),Item( 'ls_ytil_3D' ),Item( 'ls_ydp_3D' )
    traits_view = View( material_cont_ribbon, material_sf_ribbon,
                   profile_ribbon,
                   plot3D_ribbon,
                   resizable = True,
                   width = 900, height = 800,
                   title = "Plot App"
                   # ChacoPlotItem( "xdata", "ydata", y_bounds = ( -10., 10. ), y_auto = False, resizable = True, show_label = False, x_label = "x", y_label = "y", title = "" ),
                   
                  ) 

if __name__ == '__main__':
    p = Plot_app()
    p.configure_traits()

