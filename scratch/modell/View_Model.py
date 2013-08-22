
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
from stats.spirrid import make_ogrid as orthogonalize
from mayavi import mlab
from mayavi.api import Engine
from numpy import array
import os

# from quaducom.meso.homogenized_crack_bridge.elastic_matrix.hom_CB_elastic_mtrx_py_loop import CompositeCrackBridgeLoop
class Plot_app( HasTraits ):
    # MATERIAL PARAMETERS
    # reinf1
    x_name3D = Enum( 'Ll', 'Lr', 'V_f1', 'V_f2', 'V_f3', 'V_f4', 'r1', 'r2', 'r3', 'r4', 'E_f1', 'E_f2', 'E_f3', 'E_f4' )
    y_name3D = Enum( 'V_f1', 'V_f2', 'V_f3', 'V_f4' )
    Ll = Float( 100. )
    Lr = Float( 100. )
    discr_amin = Int( 50 )
    
    ###############
    # reinf1
    r1 = Float( 0.00345 )
    tau1 = List
    def _tau1_default( self ):
        return ['weibull_min', 0.006, .23, .03 ],
    V_f1 = Float( 0.03 )
    E_f1 = Float( 240e3 )
    xi_shape1 = Float( 4. )
    xi_sv01 = Float( 0.0025 )
    n_int1 = Int( 50 )
    ###############
    # reinf2
    r2 = Float( 0.00345 )
    tau2 = List
    def _tau2_default( self ):
        return ['uniform', 0.01, 0.11]
    V_f2 = Float( 0.03 )
    E_f2 = Float( 180e3 )
    xi_shape2 = Float( 4. )
    xi_sv02 = Float( 0.0025 )
    n_int2 = Int( 50 )
###########################
    # reinf3
    r3 = Float( 0.1 )
    tau3 = Float( 1. )
    lf3 = Float( 30. ), 'E_f1'
    snub3 = Float( 2. )
    phi3 = List
    def _phi3_default( self ):
        return ['sin2x', 0., 1.]
    V_f3 = Float( 0.02 )
    E_f3 = Float( 180e3 )
    xi3 = Float( 100. )
    n_int3 = Int( 50 )
    # reinf4
    r4 = Float( 0.1 )
    tau4 = Float( 1. )
    lf4 = Float( 30. )
    snub4 = Float( 2. )
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
                                 Ll = 100.,
                                 Lr = 100.,
                                 discr_amin = 70 )
    ccb_view = Instance( CompositeCrackBridgeView )
    def _ccb_view_default( self ):
        return CompositeCrackBridgeView( model = self.model )
    
    '''
    paramstring_3D = Property( depends_on = '_plot3D_launch_fired' )
    def _get_paramstring_3D( self ):
        Ll = 'LlOO{}OO'.format( self.Ll )
        Lr = '{}'.format( self.Ll )
        discr_amin = '{}'.format( self.discr_amin )
        r1 = '{}'.format( self.r1 )
        tau1 = '{}'.format( self.tau1 )
        V_f1 = '{}'.format( self.V_f1 )
        E_f1 = '{}'.format( self.E_f1 )
        xi_shape1 = '{}'.format( self.xi_shape1 )
        xi_sv01 = '{}'.format( self.xi_sv01 )
        n_int1 = '{}'.format( self.n_int1 )
        r2 = '{}'.format( self.r2 )
        tau2 = '{}'.format( self.tau2 )
        lf2 = '{}'.format( self.lf2 )
        snub2 = '{}'.format( self.snub2 )
        phi2 = '{}'.format( self.phi2 )
        V_f2 = '{}'.format( self.V_f2 )
        E_f2 = '{}'.format( self.E_f2 )
        xi2 = '{}'.format( self.xi2 )
        n_int2 = '{}'.format( self.n_int2 )
        choice3D = '{}'.format( self.choice3D )
        x_name3D = '{}'.format( self.x_name3D )
        y_name3D = '{}'.format( self.y_name3D )
        ls_xfrom_3D = '{}'.format( self.ls_xfrom_3D )
        ls_xtil_3D = '{}'.format( self.ls_xtil_3D )
        ls_xdp_3D = '{}'.format( self.ls_xdp_3D )
        ls_yfrom_3D = '{}'.format( self.ls_yfrom_3D )
        ls_ytil_3D = '{}'.format( self.ls_ytil_3D )
        ls_ydp_3D = '{}'.format( self.ls_ydp_3D )
        clamps = '{}' * 28
        return clamps.format( Ll, Lr, discr_amin,
                                                                         r1, tau1, V_f1, E_f1,
                                                                         xi_shape1, xi_sv01, n_int1,
                                                                         r2, tau2, lf2, snub2, phi2,
                                                                         V_f2, E_f2, xi2 , n_int2,
                                                                         choice3D, x_name3D , y_name3D ,
                                                                         ls_xfrom_3D, ls_xtil_3D, ls_xdp_3D,
                                                                         ls_yfrom_3D, ls_ytil_3D, ls_ydp_3D 
                                                                           )
    
    def library_read( self, methodname ):
        os.chdir( 'Library' )
        os.chdir( methodname )
        if self.paramstring_3D in os.listdir( os.curdir ):   
            x = np.load( 'x.npy' )
            y = np.load( 'y.npy' )
            z = np.load( 'z.npy' )
            os.chdir( os.pardir )
            os.chdir( os.pardir )
            return 1, x, y, z
        else:
            os.chdir( os.pardir )
            os.chdir( os.pardir )
            return 0, 0, 0, 0
    
    def library_write( self, methodname ):
        os.chdir( 'Library' )
        os.chdir( methodname )
        os.chdir( os.pardir )
        os.chdir( os.pardir )
    '''
    
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
                                tau = RV( self.tau1[0], loc = np.float( self.tau1[1] ), scale = np.float( self.tau1[2] ) ),
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
   
    choice3D = Int( 0 )
    
    ls_xfrom_3D = Float( 0.00001 )
    ls_xtil_3D = Float( 0.1 )
    ls_xdp_3D = Int( 10 )
    ls_yfrom_3D = Float( 0.00001 )
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
            print'file is not saved. Foldername already exists'
            # let exception propagate if we just can't
            # cd into the specified directory

    compute3D_plot = Property()
    @cached_property
    def _get_compute3D_plot( self ):
        # lib_list = self.library_read( 'plot3D' )
        # if lib_list[0] == True:
        #    return lib_list[1], lib_list[2], lib_list[3]
        # else: 
            z_axis_y = []
            z_axis = []
            x_string_exec = '{}{}{}'.format( 'self.', self.x_name3D, '=x' )
            y_string_exec = '{}{}{}'.format( 'self.', self.y_name3D, '=y' )
            x_ls = np.linspace( self.ls_xfrom_3D, self.ls_xtil_3D, self.ls_xdp_3D )
            y_ls = np.linspace( self.ls_yfrom_3D, self.ls_ytil_3D, self.ls_ydp_3D )
            # method = [self.ccb_view.sigma_c_max[0], self.ccb_view.sigma_c_max[1]]
            for x in x_ls:
                exec x_string_exec
                for y in y_ls:
                    exec y_string_exec
                    self._material_ok_fired()
                    self.set_model()
                    z_i = self.ccb_view.sigma_c_max[self.choice3D]
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
        mlab.surf( x, y , z, representation = 'wireframe' )  # ,warp_scale = "auto"
        if 'config.py' in os.listdir( os.curdir ):
            execfile( 'config.py' )
        else:
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
        os.chdir( os.pardir )
        os.chdir( os.pardir )
        mlab.show()
            
        
    def plot_3D( self ):
        x_axis = self.compute3D_plot[0]
        y_axis = self.compute3D_plot[1]
        z_axis = self.compute3D_plot[2]
        # print x_axis, y_axis, z_axis
        engine = Engine()
        engine.start()
        if len( engine.scenes ) == 0:
            engine.new_scene()
        mlab.surf( x_axis, y_axis , z_axis , representation = 'wireframe' )  # , warp_scale = "auto"
        
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
        mlab.zlabel( 'z' )
        mlab.show()
    
    launch_profile = Button( 'launch' )
    click_counter = Int
    def _launch_profile_fired( self ):
        self.set_model()
        self.profile()
        
    material_ok = Button( 'OK' )
    def _material_ok_fired( self ):
        pass
    
    plot3D_launch = Button( 'launch' )
    def _plot3D_launch_fired( self ):
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
    material_cont_ribbon = Group( Item( name = 'reinf1_bool', label = 'reinf1' ),
                            Item( name = 'r1', label = 'r' ),
                             Item( name = 'V_f1', label = 'Vf' ),
                             Item( name = 'tau1' , label = 'tau' ),
                             Item( name = 'E_f1', label = 'Ef' ),
                             Item( name = 'xi_shape1', label = 'xi_shape' ),
                             Item( name = 'xi_sv01', label = 'xi_sv0' ),
                             Item( name = 'n_int1', label = 'integration_points' ),
                             Item( '_' ),
                             Item( name = 'reinf2_bool', label = 'reinf2' ),
                             Item( name = 'r2', label = 'r' ),
                             Item( name = 'V_f2', label = 'Vf' ),
                             Item( name = 'tau2' , label = 'tau' ),
                             Item( name = 'E_f2', label = 'Ef' ),
                             Item( name = 'xi_shape2', label = 'xi_shape' ),
                             Item( name = 'xi_sv02', label = 'xi_sv0' ),
                             Item( name = 'n_int2', label = 'integration_points' ),
                               Item( '_' ),
                               'material_ok', label = 'Cont Fiber Reinf' )
    material_sf_ribbon = Group( Item( name = 'reinf3_bool', label = 'reinf3' ),
                            Item( name = 'r3', label = 'r' ),
                              Item( name = 'tau3', label = 'tau' ),
                              Item( name = 'phi3' , label = 'phi' ) ,
                               Item( name = 'V_f3', label = 'Vf' ),
                               Item( name = 'snub3' , label = 'snub' ),
                               Item( name = 'E_f3', label = 'Ef' ),
                               Item( name = 'xi3', label = 'xi' ),
                               Item( name = 'n_int3', label = 'integration_points' ),
                             Item( name = 'reinf4_bool', label = 'reinf4' ),
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

