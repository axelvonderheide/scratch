from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

from enthought.chaco.chaco_plot_container_editor import \
     PlotContainerEditor
from enthought.chaco.tools.api import \
     PanTool, SimpleZoom
from enthought.chaco.api import \
     Plot, AbstractPlotData, ArrayPlotData

from numpy import \
     array, arange, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt, swapaxes, copy, loadtxt

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from os.path import join

from scipy.linalg import \
    eigh, inv

from ibvpy.core.tstepper import \
    TStepper as TS

from ibvpy.mats.mats_eval import \
    IMATSEval, MATSEval

from ibvpy.core.rtrace_eval import \
    RTraceEval

from ibvpy.api import \
    RTraceGraph,RTraceArraySnapshot

from mats2D_cmdm_polar_discr import \
    IPolarFn, IsotropicPolarFn, AnisotropicPolarFn

from ibvpy.mats.mats2D.mats2D_tensor import \
    map2d_eps_eng_to_mtx, map2d_sig_eng_to_mtx, map2d_eps_mtx_to_eng, map2d_sig_mtx_to_eng, \
    map2d_ijkl2mn, map2d_tns2_to_tns4, map2d_tns4_to_tns2, compliance_mapping2d, \
    get_D_plane_stress, get_D_plane_strain, get_C_plane_stress, get_C_plane_strain

from ibvpy.mats.mats3D.mats3D_tensor import \
    map3d_eps_eng_to_mtx, map3d_sig_eng_to_mtx, map3d_eps_mtx_to_eng, map3d_sig_mtx_to_eng, \
    map3d_ijkl2mn, map3d_tns2_to_tns4, map3d_tns4_to_tns2, compliance_mapping3d


#---------------------------------------------------------------------------
# Material time-step-evaluator for Microplane-Damage-Model
#---------------------------------------------------------------------------

class MA2DMicroplaneDamage( MATSEval ):
    '''
    Microplane Damage Model.
    '''
    implements( IMATSEval )

    #------------------------------------------------------------------------------------
    # Damage function specification
    #------------------------------------------------------------------------------------
    polar_fn_class = Trait( 'Isotropic polar function',#'Anisotropic polar function'
                            {'Isotropic polar function'   : IsotropicPolarFn,
                             'Anisotropic polar function' : AnisotropicPolarFn } )
    
    polar_fn = Property( Instance( IPolarFn ), depends_on = 'polar_fn_class' )
    @cached_property
    def _get_polar_fn(self):
        print 'MA2DMicroplaneDamage: get_polar_fn'
        return self.polar_fn_class_()
    
    def _set_polar_fn(self, value):
        print 'MA2DMicroplaneDamage: set_polar_fn'
        if not (self.polar_fn_class_ is value.__class__):
            raise ValueError, 'class mismatch %s != %s' % ( self.polar_fn_class_, value.__class__)
        self._polar_fn = value
            
    #---------------------------------------------------------------------------
    # Material parameters 
    #---------------------------------------------------------------------------

    E = Delegate( 'polar_fn' )
    nu = Delegate( 'polar_fn' )
    c_T = Delegate( 'polar_fn' )
    alpha_list = Delegate( 'polar_fn' )

    #---------------------------------------------------------------------------
    # Parameters of the numerical algorithm (integration)
    #---------------------------------------------------------------------------

    n_mp = Delegate('polar_fn')
    
    model_version   = Enum("compliance", "stiffness")
    stress_state    = Enum("plane_stress", "plane_strain")
    symmetrization  = Enum("product-type", "sum-type")

    dimensionality = Trait( '2D', {'2D' : 2, '3D' : 3 },
                 label = 'Dimensionality',
                 desc = 'Specify the dimensionality',
                 auto_set = False)
    
    @on_trait_change('n_mp, dimensionality')
    def update_arrays( self ):
        self._init_arrays()
        self.changed = True
        
    elastic_debug = Bool( False,
                          desc = 'Switch to elastic behavior - used for debugging',
                          auto_set = False)
    
    double_constraint = Bool( False,
                              desc = 'Use double constraint to evaluate microplane elastic and fracture energy (Option effects only the response tracers)',
                              auto_set = False)
   

    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event                     

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View(  Group( 
                                Group(
                                    Item('polar_fn_class', style='custom', show_label = False),
                                      Item('polar_fn', style='custom', show_label = False, resizable = True ),
                                      label='XXX Material parameters',
                                 ),
                                Group( Item('model_version', style = 'custom' ),
                                       Item('stress_state', style = 'custom' ),
                                       Item('symmetrization', style = 'custom' ),
                                       Item('dimensionality', style='custom'),
                                       Item('elastic_debug@'),
                                       Item('double_constraint@'),
                                       Spring(resizable = True),
                                       label='Configuration parameters', show_border=True,
                                       ),
                                layout = 'tabbed',
                                orientation = 'horizontal',
                                springy = True
                                ),
                        kind = 'modal',
                        resizable = True,
                        scrollable = True,
                        width = 0.6, height = 0.8,
                        buttons = ['OK','Cancel' ]
                        )

    #-----------------------------------------------------------------------------------------------
    # Private initialization methods
    #-----------------------------------------------------------------------------------------------

    def __init__( self, **kwtraits ):
        '''
        Subsidiary arrays required for the integration.
        they are set up only once prior to the computation
        '''
        super( MA2DMicroplaneDamage, self ).__init__( **kwtraits )
        self._init_arrays()
        print 'MA2DMicroplaneDamage: __init__ called'

    def _init_arrays( self ):
        '''
        Redefine the arrays:
        MPN = microplane normals
        MPW = microplane weights
        MPNN = dyadic product of microplane weights
        '''
        # check the dimensionality and assign the attribute 'n_dim' to self
        self._check_dimensionality()
        # 2D-case:
        if self.n_dim == 2:
            # get the normal vectors of the microplanes
            self._MPN = array([[ cos( alpha ), sin( alpha )] for alpha in self.alpha_list ])
            # get the weights of the microplanes
            self._MPW = ones(self.n_mp) / self.n_mp * 2
            # get the dyadic product of the microplane normals
            self._MPNN = array( [ outer(mpn,mpn) for mpn in self._MPN ] )
        # 3D-case:
        else:
            # Use the numerical integration formula for the n-sphere by STROUD (1979)
            # with a number of microplanes equal to 28:
            self.n_mp = 28
            # microplane normals:
            self._MPN = array([[.577350259, .577350259, .577350259], \
                         [.577350259, .577350259,-.577350259], \
                         [.577350259,-.577350259, .577350259], \
                         [.577350259,-.577350259,-.577350259], \
                         [.935113132, .250562787, .250562787], \
                         [.935113132, .250562787,-.250562787], \
                         [.935113132,-.250562787, .250562787], \
                         [.935113132,-.250562787,-.250562787], \
                         [.250562787, .935113132, .250562787], \
                         [.250562787, .935113132,-.250562787], \
                         [.250562787,-.935113132, .250562787], \
                         [.250562787,-.935113132,-.250562787], \
                         [.250562787, .250562787, .935113132], \
                         [.250562787, .250562787,-.935113132], \
                         [.250562787,-.250562787, .935113132], \
                         [.250562787,-.250562787,-.935113132], \
                         [.186156720, .694746614, .694746614], \
                         [.186156720, .694746614,-.694746614], \
                         [.186156720,-.694746614, .694746614], \
                         [.186156720,-.694746614,-.694746614], \
                         [.694746614, .186156720, .694746614], \
                         [.694746614, .186156720,-.694746614], \
                         [.694746614,-.186156720, .694746614], \
                         [.694746614,-.186156720,-.694746614], \
                         [.694746614, .694746614, .186156720], \
                         [.694746614, .694746614,-.186156720], \
                         [.694746614,-.694746614, .186156720], \
                         [.694746614,-.694746614,-.186156720]]) 
            # microplane weights:
            # (values in the array must be multiplied by 6 (cf. [Baz05]))
            self._MPW = array([.0160714276, .0160714276, .0160714276, .0160714276, .0204744730, \
                         .0204744730, .0204744730, .0204744730, .0204744730, .0204744730, \
                         .0204744730, .0204744730, .0204744730, .0204744730, .0204744730, \
                         .0204744730, .0158350505, .0158350505, .0158350505, .0158350505, \
                         .0158350505, .0158350505, .0158350505, .0158350505, .0158350505, \
                         .0158350505, .0158350505, .0158350505 ]) * 6.0
            # dyadic product of the microplane normals
            self._MPNN = array( [ outer(mpn,mpn) for mpn in self._MPN ] )
            
         
    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------

    def get_state_array_size( self ):
        # In the state array the largest equivalent microplane 
        # strains reached in the loading history are saved
        return self.n_mp


    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        state_arr_size = self.get_state_array_size()
        sctx.mats_state_array = zeros(state_arr_size, 'float_')
        #
        self.update_state_on = False
        

    def _check_dimensionality( self ):
        if self.dimensionality == '2D':
            print 'dimensionality: 2D-case'
            self.n_dim = 2
            self.n_eng = 3
        elif self.dimensionality == '3D':
            print 'dimensionality: 3D-case'
            self.n_dim = 3    
            self.n_eng = 6
        

    elasticity_tensors = Property( depends_on = 'E, nu, dimensionality, stress_state' )
    @cached_property        
    def _get_elasticity_tensors( self ):
        '''
        Intialize the fourth order elasticity tensor 
        for 3D or 2D plane strain or 2D plane stress
        '''
        # ----------------------------------------------------------------------------
        # Lame constants calculated from E and nu
        # ----------------------------------------------------------------------------
        E   = self.E
        nu  = self.nu
        
        # first Lame paramter
        la = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
        # second Lame parameter (shear modulus)
        mu = E / ( 2 + 2 * nu ) 

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 3D-case
        # -----------------------------------------------------------------------------------------------------
        D4_e_3D = zeros((3,3,3,3),dtype=float)
        C4_e_3D = zeros((3,3,3,3),dtype=float)
        delta = identity(3)
        for i in range(0,3):
            for j in range(0,3):
                for k in range(0,3):
                    for l in range(0,3):
                        D4_e_3D[i,j,k,l] = la * delta[i,j] * delta[k,l] + \
                                            mu * ( delta[i,k] * delta[j,l] + delta[i,l] * delta[j,k] )
                        # elastic compliance tensor:                 
                        C4_e_3D[i,j,k,l] = (1+nu)/(2*E) * \
                                            ( delta[i,k] * delta[j,l] + delta[i,l]* delta[j,k] ) - \
                                              nu / E * delta[i,j] * delta[k,l]                    
       
        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 2D-case
        # -----------------------------------------------------------------------------------------------------
        # 1. step: Get the (6x6)-elasticity and compliance matrices 
        #          for the 3D-case:
        D2_e_3D = map3d_tns4_to_tns2(D4_e_3D)
        C2_e_3D = map3d_tns4_to_tns2(C4_e_3D)

        # 2. step: Get the (3x3)-elasticity and compliance matrices 
        #          for the 2D-cases plane stress and plane strain:
        D2_e_2D_plane_stress = get_D_plane_stress( D2_e_3D )
        D2_e_2D_plane_strain = get_D_plane_strain( D2_e_3D )
        C2_e_2D_plane_stress = get_C_plane_stress( C2_e_3D )
        C2_e_2D_plane_strain = get_C_plane_strain( C2_e_3D )

        # (3D case):                                       
        if self.dimensionality == '3D':
            D2_e = D2_e_3D
        # (2D case):
        if self.dimensionality == '2D' and self.stress_state == 'plane_stress':
            D2_e = D2_e_2D_plane_stress  

        if self.dimensionality == '2D' and self.stress_state == 'plane_strain':
            D2_e = D2_e_2D_plane_strain 
            
        # 3. step: Get the fourth order elasticity and compliance tensors
        #          for the 2D-cases plane stress and plane strain (D4.shape = (2,2,2,2)) 
        D4_e_2D_plane_stress = map2d_tns2_to_tns4( D2_e_2D_plane_stress )
        D4_e_2D_plane_strain = map2d_tns2_to_tns4( D2_e_2D_plane_strain )
        C4_e_2D_plane_stress = map2d_tns2_to_tns4( C2_e_2D_plane_stress )
        C4_e_2D_plane_strain = map2d_tns2_to_tns4( C2_e_2D_plane_strain )

        # -----------------------------------------------------------------------------------------------------
        # assign the fourth order elasticity and compliance tensors as return values
        # -----------------------------------------------------------------------------------------------------
        # (3D case):                                       
        if self.dimensionality == '3D':
            D4_e = D4_e_3D
            C4_e = C4_e_3D

        # (2D case):
        if self.dimensionality == '2D' and self.stress_state == 'plane_stress':
            #print 'stress state:   plane-stress'
            D4_e = D4_e_2D_plane_stress 
            C4_e = C4_e_2D_plane_stress 

        if self.dimensionality == '2D' and self.stress_state == 'plane_strain':
            #print 'stress state:   plane-strain'
            D4_e = D4_e_2D_plane_strain
            C4_e = C4_e_2D_plane_strain

        return D4_e, C4_e, D2_e


    D4_e = Property
    def _get_D4_e(self):
        '''
        Return the elasticity tensor
        '''
        return self.elasticity_tensors[0]


    C4_e = Property
    def _get_C4_e(self):
        '''
        Return the elastic compliance tensor
        '''
        return self.elasticity_tensors[1]

    # -----------------------------------------------------------------------------------------------------
    # Get the 3x3-elasticity and compliance matrix (used only for elastic debug)
    # -----------------------------------------------------------------------------------------------------
    D2_e = Property( depends_on = 'stress_state, E, nu, dimensionality' )
    @cached_property
    def _get_D2_e( self ):
        print 'elastic_debug:  switched on!' 
        return self.elasticity_tensors[2]

        
    identity_tns = Property( depends_on = 'dimensionality' )
    @cached_property
    def _get_identity_tns(self):
        '''
        Get the identity matrix (used only in formula for sum-type symmetrization)
        '''     
        # (2D case):                                       
        if self.dimensionality == '2D':
            return identity(2)
        else:
            return identity(3)

        
    #-----------------------------------------------------------------------------------------------------
    # Prepare evaluation - get the 4th order damage tensor beta4 (or damage effect tensor M4)
    #-----------------------------------------------------------------------------------------------------
    
    def _get_e_vct_arr(self, eps_eng ):
        '''
        Projects the strain tensor onto the microplanes and returns a list of microplane strain
        vectors. Method is used both by stiffness and compliance version to derive the list 'phi_arr'
        or 'psi_arr'! In case of the compliance version the kinematic constraint is not assumed in 
        the derivation of the formula for the damage compliance tensor 'C_mdm', e.g. the construction
        of the damage effect tensor 'M4'.
        '''
        # Switch from engineering notation to tensor notation for the apparent strains
        if self.n_dim == 2:
            # 2D-case:
            eps_mtx = map2d_eps_eng_to_mtx(eps_eng)
        else:
            # 3D-case:
            eps_mtx = map3d_eps_eng_to_mtx(eps_eng)
        # Projection of apparent strain onto the individual microplanes
        # slower: e_vct_arr = array( [ dot( eps_mtx, mpn ) for mpn in self._MPN ] )
        # slower: e_vct_arr = transpose( dot( eps_mtx, transpose(self._MPN) ))
        # due to the symmetry of the strain tensor eps_mtx = transpose(eps_mtx) and so this is equal
        e_vct_arr = dot( self._MPN, eps_mtx ) 
        return e_vct_arr

    
    def _get_e_equiv_arr(self, e_vct_arr ):
        '''
        Returns a list of the microplane equivalent strains 
        based on the list of microplane strain vectors
        '''        
        # magnitude of the normal strain vector for each microplane
        #@todo: faster numpy functionality possible?
        e_N_arr = array([ dot( e_vct, mpn  ) for e_vct, mpn in zip(e_vct_arr, self._MPN) ])
        # positive part of the normal strain magnitude for each microplane
        e_N_pos_arr = ( fabs( e_N_arr ) + e_N_arr ) / 2
        # normal strain vector for each microplane
        #@todo: faster numpy functionality possible?
        e_N_vct_arr = array( [ self._MPN[i,:] * e_N_arr[i] for i in range(0,self.n_mp) ] )
        # tangent strain ratio
        c_T = self.c_T
        # tangential strain vector for each microplane
        e_T_vct_arr = e_vct_arr - e_N_vct_arr
        # squared tangential strain vector for each microplane
        e_TT_arr = array( [ inner( e_T_vct, e_T_vct ) for e_T_vct in e_T_vct_arr ] )
        # equivalent strain for each microplane
        e_equiv_arr = arr_sqrt( e_N_pos_arr * e_N_pos_arr + c_T * e_TT_arr ) 
        return e_equiv_arr


    def _get_e_max( self, e_equiv_arr, e_max_arr ):
        '''
        Compares the equivalent microplane strain of a single microplane with the 
        maximum strain reached in the loading history for the entire array
        '''
        bool_e_max = e_equiv_arr >= e_max_arr
        e_max_arr[bool_e_max] = e_equiv_arr[bool_e_max]
        return e_max_arr
    

    def _get_state_variables(self, sctx, eps_app_eng):
        '''
        Compares the list of current equivalent microplane strains with 
        the values in the state array and returns the higher values 
        '''
        e_vct_arr = self._get_e_vct_arr(eps_app_eng)
        e_equiv_arr = self._get_e_equiv_arr( e_vct_arr )
        e_max_arr_old = sctx.mats_state_array
        e_max_arr_new = self._get_e_max(e_equiv_arr, e_max_arr_old)
        # Get the highest value from the current state array.
#        self.e_max_value = max( e_max_arr_new )
        return e_max_arr_new

    e_max_value = Float(0.0)

    
    def _get_phi_arr(self, sctx, eps_app_eng ):
        '''
        Returns a list of the integrity factors for all microplanes.
        '''
        e_max_arr = self._get_state_variables( sctx, eps_app_eng )    
        return self.polar_fn.get_phi_arr( sctx, e_max_arr )
        

    def _get_phi_mtx(self, sctx, eps_app_eng):
        '''
        Returns the 2nd order damage tensor 'phi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_arr   = self._get_phi_arr( sctx, eps_app_eng )
        # integration terms for each microplanes
        #@todo: faster numpy functionality possible?
        phi_mtx_arr = array( [ phi_arr[i] * self._MPNN[i,:,:] * self._MPW[i]
                                for i in range(0,self.n_mp) ] )
        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        phi_mtx = phi_mtx_arr.sum(0) 
        return phi_mtx


    def _get_psi_mtx(self, sctx, eps_app_eng):
        '''
        Returns the 2nd order damage effect tensor 'psi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_arr   = self._get_phi_arr( sctx, eps_app_eng )
        psi_mtx_arr = array( [ 1. / phi_arr[i] * self._MPNN[i,:,:] * self._MPW[i]
                                for i in range(0,self.n_mp) ] )
        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        psi_mtx = psi_mtx_arr.sum(0)
        return psi_mtx


    def _get_beta_tns_product_type(self, phi_mtx):
        '''
        Returns the 4th order damage tensor 'beta4' using product-type symmetrization
        '''
        # (cf. [Baz97], Eq.(87))
        n_dim = self.n_dim
        # Get the direction of the principle damage coordinates (pdc):
        phi_eig_value, phi_eig_mtx = eigh( phi_mtx )
        phi_eig_value_real = array([ pe.real for pe in phi_eig_value] )                
        # transform phi_mtx to PDC:
        # (assure that besides the diagonal the entries are exactly zero)
        phi_pdc_mtx = zeros((n_dim,n_dim),dtype=float)
        for i in range(n_dim):
            phi_pdc_mtx[i,i] = phi_eig_value_real[i]
        # w_mtx = tensorial square root of the second order damage tensor:
        w_pdc_mtx = arr_sqrt( phi_pdc_mtx )
        # transform the matrix w back to x-y-coordinates:
        w_mtx = dot(dot(phi_eig_mtx, w_pdc_mtx),transpose(phi_eig_mtx))
        # beta_ijkl = w_ik * w_jl (cf. [Baz 97])
        # exploiting numpy-functionality (faster).
        # Method 'outer' returns beta_ijkl = w_ij * w_kl,
        # therefore the axis j and k need to be swapped
        beta4_ = outer(w_mtx, w_mtx).reshape(n_dim,n_dim,n_dim,n_dim)
        beta4 = beta4_.swapaxes(1,2)
        return beta4

    
    def _get_beta_tns_sum_type(self, phi_mtx):
        '''
        Returns the 4th order damage tensor 'beta4' using sum-type symmetrization
        '''
        # (cf. [Jir99], Eq.(21))
        n_dim = self.n_dim
        delta = self.identity_tns
        # use numpy functionality to evaluate [Jir99], Eq.(21) 
        beta_ijkl = outer(phi_mtx, delta).reshape(n_dim,n_dim,n_dim,n_dim)
        beta_ikjl = beta_ijkl.swapaxes(1,2)
        beta_iljk = beta_ikjl.swapaxes(2,3)
        beta_jlik = beta_iljk.swapaxes(0,1)
        beta_jkil = beta_jlik.swapaxes(2,3)
        beta4 = 0.25 * ( beta_ikjl + beta_iljk + beta_jkil + beta_jlik )
        return beta4


    def _get_M_tns_product_type(self, psi_mtx):
        '''
        Returns the 4th order damage effect tensor 'M4' using product-type symmetrization
        '''
        n_dim = self.n_dim
        # Get the direction orthogonal to the principle damage coordinates (pdc):
        # @todo: is this direction orthogonal? Which one do we want?
        psi_eig_value, psi_eig_mtx = eigh( psi_mtx )
        psi_eig_value_real = array([ pe.real for pe in psi_eig_value] )
        # transform phi_mtx to PDC:
        # (assure that besides the diagonal the entries are exactly zero)
        psi_pdc_mtx = zeros((n_dim, n_dim),dtype=float)
        for i in range(n_dim):
            psi_pdc_mtx[i,i] = psi_eig_value_real[i]
        # second order damage effect tensor:
        w_hat_pdc_mtx = arr_sqrt( psi_pdc_mtx )
        # transform the matrix w back to x-y-coordinates:
        w_hat_mtx = dot(dot(psi_eig_mtx, w_hat_pdc_mtx),transpose(psi_eig_mtx))
        # M_ijkl = w_hat_ik * w_hat_lj (cf. Eq.(5.62) Script Prag Jirasek (2007))
        #        = w_hat_ik * w_hat_jl (w is a symmetric tensor) 
        # Exploiting numpy-functionality using the
        # method 'outer' (returns M_ijkl = w_hat_ij * w_hat_kl),
        # therefore the axis j and k need to be swapped
        M4_ = outer(w_hat_mtx, w_hat_mtx).reshape(n_dim,n_dim,n_dim,n_dim)
        M4 = M4_.swapaxes(1,2)
        return M4

        
    def _get_M_tns_sum_type(self, psi_mtx):
        '''
        Returns the 4th order damage effect tensor 'M4' using sum-type symmetrization
        '''
        n_dim = self.n_dim
        delta = self.identity_tns
        # use numpy functionality to evaluate [Jir99], Eq.(21) 
        M_ijkl = outer(psi_mtx, delta).reshape(n_dim,n_dim,n_dim,n_dim)
        M_ikjl = M_ijkl.swapaxes(1,2)
        M_iljk = M_ikjl.swapaxes(2,3)
        M_jlik = M_iljk.swapaxes(0,1)
        M_jkil = M_jlik.swapaxes(2,3)
        M4 = 0.25 * ( M_ikjl + M_iljk + M_jkil + M_jlik )
        return M4


    #-----------------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------------
    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1, eps_avg = None ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        # -----------------------------------------------------------------------------------------------
        # for debugging purposes only: if elastic_debug is switched on, linear elastic material is used 
        # -----------------------------------------------------------------------------------------------
        if self.elastic_debug:
            # NOTE: This must be copied otherwise self.D2_e gets modified when
            # essential boundary conditions are inserted
            D2_e = copy( self.D2_e )
            sig_eng = tensordot( D2_e, eps_app_eng, [[1],[0]])
            return sig_eng, D2_e
       
        # -----------------------------------------------------------------------------------------------
        # update state variables
        # -----------------------------------------------------------------------------------------------
        if sctx.update_state_on:
#            print 'in cmdm: update_state_array!!!!!!!!!!!!!!'
            eps_n = eps_app_eng - d_eps
            e_max = self._get_state_variables( sctx, eps_n)
            sctx.mats_state_array[:] = e_max
            
        #------------------------------------------------------------------------------------------------
        # stiffness version:
        #------------------------------------------------------------------------------------------------
        if self.model_version == 'stiffness' :

            #------------------------------------------------------------------------------------
            # Damage tensor (2th order):  
            #------------------------------------------------------------------------------------
            phi_mtx = self._get_phi_mtx(sctx, eps_app_eng)
            #------------------------------------------------------------------------------------
            # Damage tensor (4th order) using product-type symmetrization: 
            #------------------------------------------------------------------------------------
            if self.symmetrization == 'product-type':
                # (cf. [Baz97])
                beta4 = self._get_beta_tns_product_type(phi_mtx)
                            
            #-------------------------------------------------------------------------------------
            # Damage tensor (4th order) using sum-type symmetrization: 
            #-------------------------------------------------------------------------------------
            elif self.symmetrization == 'sum-type':
                # (cf. [Jir99] Eq.(21))
                beta4 = self._get_beta_tns_sum_type(phi_mtx)
                
            #-------------------------------------------------------------------------------------
            # Damaged stiffness tensor calculated based on the damage tensor beta4:
            #-------------------------------------------------------------------------------------
            # (cf. [Jir99] Eq.(7): C = beta * D_e * beta^T), 
            # minor symmetry is tacitly assumed ! (i.e. beta_ijkl = beta_jilk)
            D4_mdm = tensordot( beta4, tensordot( self.D4_e, beta4, [[2,3],[2,3]] ), [[2,3],[0,1]] )

            #-------------------------------------------------------------------------------------
            # Reduction of the fourth order tensor to a matrix assuming minor and major symmetry:
            #-------------------------------------------------------------------------------------
            # 2D-case:
            if self.n_dim == 2:
                D2_mdm = map2d_tns4_to_tns2(D4_mdm)
            # 3D-case:
            else:
                D2_mdm = map3d_tns4_to_tns2(D4_mdm)
     
                    
        #------------------------------------------------------------------------------------------------
        # compliance version:
        #------------------------------------------------------------------------------------------------
        elif self.model_version == 'compliance' :

            #------------------------------------------------------------------------------------
            # Damage effect tensor (2th order):  
            #------------------------------------------------------------------------------------
            psi_mtx = self._get_psi_mtx(sctx, eps_app_eng)

            #------------------------------------------------------------------------------------
            # Damage effect tensor (4th order) using product-type-symmetrization:  
            #------------------------------------------------------------------------------------
            if self.symmetrization == 'product-type':
                M4 = self._get_M_tns_product_type(psi_mtx)

            #------------------------------------------------------------------------------------
            # Damage effect tensor (4th order) using sum-type-symmetrization:  
            #------------------------------------------------------------------------------------
            elif self.symmetrization == 'sum-type':
                # (cf. [Jir99], Eq.(30)
                M4 = self._get_M_tns_sum_type(psi_mtx)

            #-------------------------------------------------------------------------------------
            # Damaged compliance tensor calculated based on the damage effect tensor M4:
            #-------------------------------------------------------------------------------------
            # (cf. [Jir99] Eq.(8): C = M^T * C_e * M, 
            # minor symmetry is tacitly assumed ! (i.e. M_ijkl = M_jilk)
            C4_mdm = tensordot( M4, tensordot( self.C4_e, M4, [[2,3],[0,1]] ), [[0,1],[0,1]] )

            #-------------------------------------------------------------------------------------
            # Reduction of the fourth order tensor to a matrix assuming minor and major symmetry: 
            #-------------------------------------------------------------------------------------
            # 2D-case:
            if self.n_dim == 2:
                C2_mdm = map2d_tns4_to_tns2(C4_mdm)
                # multiply the resulting D matrix with the factors for compliance mapping
                D2_mdm = inv(compliance_mapping2d( C2_mdm ))
            # 3D-case:
            else:
                C2_mdm = map3d_tns4_to_tns2(C4_mdm)
                # multiply the resulting D matrix with the factors for compliance mapping
                D2_mdm = inv(compliance_mapping3d( C2_mdm ))

        #----------------------------------------------------------------------------------------
        # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
        #----------------------------------------------------------------------------------------
        sig_eng = tensordot( D2_mdm, eps_app_eng, [[1],[0]])
        return sig_eng, D2_mdm


    #---------------------------------------------------------------------------------------------
    # Control variables and update state method
    #---------------------------------------------------------------------------------------------

    def new_cntl_var(self):
        return zeros( self.n_eng, float_ )

    def new_resp_var(self):
        return zeros( self.n_eng, float_ )

    def update_state(self, sctx, eps_app_eng ):
        '''
        Update state method is called upon an accepted time-step. 
        Here just set the flag on to make the update afterwards in the method itself.
        '''
        print 'in update-state'
        self.update_state_on = True


    #---------------------------------------------------------------------------------------------
    # Response trace evaluators
    #---------------------------------------------------------------------------------------------

    def get_eps_app( self, sctx, eps_app_eng ):
        if self.dimensionality == '2D':
            return map2d_eps_eng_to_mtx( eps_app_eng )
        else: # 3D
            return map3d_eps_eng_to_mtx( eps_app_eng )
    
    def get_sig_app( self, sctx, eps_app_eng ):
        # @TODO
        # the stress calculation is performed twice - it might be
        # cached. But not in the spatial integration scheme.
        sig_app_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        #return sig_app_eng
        if self.dimensionality == '2D':
            return map2d_sig_eng_to_mtx( sig_app_eng )
        else: # 3D
            return map3d_sig_eng_to_mtx( sig_app_eng )
    
    def get_microplane_integrity(self, sctx, eps_app_eng ):
        phi_arr = self._get_phi_arr(sctx, eps_app_eng)
        return phi_arr

    def get_sig_norm( self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return array( [ scalar_sqrt( sig_eng[0]**2 + sig_eng[1]**2 ) ] )
    
    def get_phi_mtx( self, sctx, eps_app_eng ):
        return self._get_phi_mtx( sctx, eps_app_eng )
    
    def get_phi_pdc( self, sctx, eps_app_eng ):
        phi_mtx = self._get_phi_mtx( sctx, eps_app_eng )
        # Get the direction of the principle damage coordinates (pdc):
        phi_eig_value, phi_eig_mtx = eigh( phi_mtx )
        phi_eig_value_real = array([ pe.real for pe in phi_eig_value] )
        phi_pdc_mtx = zeros([self.n_dim,self.n_dim])
        for i in range(self.n_dim):
            phi_pdc_mtx[i,i] = phi_eig_value_real[i]
        return array( [ phi_pdc_mtx[0,0] ])

# ------------------------------------------
# SUBSIDARY METHODS used only for the response tracer:
# ------------------------------------------

    ### Projection: Methods used for projection: 

    #  '_get_e_vct_arr'  -  has been defined above (see 'get_corr_pred') 
    #   (also used as subsidary method for '_get_e_s_vct_arr'.) 
    
    def _get_s_vct_arr(self, sig_eng ):
        '''
        Projects the stress tensor onto the microplanes and
        returns a list of microplane strain vectors.
        (Subsidary method for '_get_e_s_vct_arr'.) 
        '''
        # Switch from engineering notation to tensor notation for the apparent strains
        if self.n_dim == 2:
            # 2D-case:
            sig_mtx = map2d_sig_eng_to_mtx(sig_eng)
        else:
            # 3D-case:
            sig_mtx = map3d_sig_mtx_to_eng(sig_eng)
        # Projection of apparent strain onto the individual microplanes
        # slower: s_vct_arr = array( [ dot( sig_mtx, mpn ) for mpn in self._MPN ] )
        s_vct_arr = dot( self._MPN, sig_mtx ) 
        return s_vct_arr
    
    
    ### Equiv: Equivalent microplane parts:
        
    #  '_get_e_equiv_arr'  -  has been defined above (see 'get_corr_pred') 

    def _get_s_equiv_arr(self, s_vct_arr ):
        '''
        Get microplane equivalent stress. 
        (Subsidary method for '_get_e_s_equiv_arr'.) 
        '''        
        # The same method is used to calculate the equivalent stresses
        # and the equivalent strains
        s_equiv_arr = self._get_e_equiv_arr(s_vct_arr)
        return s_equiv_arr
    
    
    ### N: normal microplane parts    

    def _get_e_N_arr(self, e_vct_arr):
        '''
        Returns a list of the microplane normal strains (scalar) 
        based on the list of microplane strain vectors
        (Subsidary method for '_get_e_s_N_arr'.) '''                
        # magnitude of the normal strain vector for each microplane
        e_N_arr = array([ dot( e_vct, mpn  ) for e_vct, mpn in zip(e_vct_arr, self._MPN) ])
        return e_N_arr
        # @todo: check if the direct calculation of e_N using MPNN makes sense. 
        #        Note that e_T needs the calculation of e_N_vct ! 
        #        s_N_arr = array( [ s_vct_arr[i] * self._MPNN[i,:,:] for i in range(0,self.n_mp) ] )
        #        return s_N_arr[0,:].flatten()

    def _get_s_N_arr(self, s_vct_arr):
        # the same method is used for the calculation of the 
        # normal parts of the microplane strain and stress vector
        s_N_arr = self._get_e_N_arr(s_vct_arr)
        return s_N_arr

    ### T: tangential microplane parts    

    def _get_e_T_arr(self, e_vct_arr):
        '''
        Returns a list of the microplane shear strains (scalar)
        based on the list of microplane strain vectors
        (Subsidary method for '_get_e_s_T_arr'.) 
        '''        
        # magnitude of the normal strain vector for each microplane
        e_N_arr = self._get_e_N_arr(e_vct_arr)
        # normal strain vector for each microplane
        e_N_vct_arr = array( [ self._MPN[i,:] * e_N_arr[i] for i in range(0,self.n_mp) ] )
        # tangential strain vector for each microplane
        e_T_vct_arr = e_vct_arr - e_N_vct_arr
        # squared tangential strain vector for each microplane
        e_TT_arr = array( [ inner( e_T_vct, e_T_vct ) for e_T_vct in e_T_vct_arr ] )
        # equivalent strain for each microplane
        e_T_arr = arr_sqrt( e_TT_arr )      
        return e_T_arr

    def _get_s_T_arr(self, s_vct_arr):
        # the same method is used for the calculation of the 
        # tangential parts of the microplane strain and stress vector
        s_T_arr = self._get_e_T_arr(s_vct_arr)
        return s_T_arr
    
    
# -------------------------------------------------------------
# Get the microplane strain and stresses either based on a 
# consistently derived model-version or the double constraint
# -------------------------------------------------------------
    
    def _get_e_s_vct_arr(self, sctx, eps_app_eng):
        '''
        Returns two arrays containing the microplane strain and stress vectors 
        either assuming a double constraint or consistently derived based on the specified model version
        '''
        #----------------------------------------------------------------------------
        # DOUBLE CONSTRAINT
        #----------------------------------------------------------------------------
        if self.double_constraint:
            '''
            Return a pair of microplane stress and strain vectors based on the 
            double constrain, i.e kinematic AND static constraint !
            The connection between the apparent strain and stress tensor is established
            with D_mdm based on the chosen model version (e.g. stiffness or compliance) 
            '''
            # microplane strain vectors obtained by projection (kinematic constraint):
            e_app_vct_arr   = self._get_e_vct_arr( eps_app_eng )
            # get the corresponding macroscopic stresses
            sig_app_eng, D_mtx  = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
            # microplane stress vectors obtained by projection (static constraint)
            s_app_vct_arr   = self._get_s_vct_arr(sig_app_eng)
            return e_app_vct_arr, s_app_vct_arr
        
        #----------------------------------------------------------------------------
        # CONSISTENTLY DERIVED pair of microplane strain and stress vectors
        #----------------------------------------------------------------------------
        '''
        Returns two arrays containing the microplane strain and stress vectors 
        consistently derived based on the specified model version, i.e. either 'stiffness' or 'compliance'
        '''
        #-------------------
        # stiffness version:
        #-------------------
        if self.model_version == 'stiffness' :
            ### microplane equivalent strains obtained by projection (kinematic constraint)
            e_app_vct_arr   = self._get_e_vct_arr( eps_app_eng )
            
            ### microplane equivalent stresses calculated based on corresponding 'beta' and 'phi_mtx'
            # 2nd order damage tensor:
            phi_mtx = self._get_phi_mtx(sctx, eps_app_eng)
            # 4th order damage tensor:
            if self.symmetrization == 'product-type':
                beta4 = self._get_beta_tns_product_type(phi_mtx)
            elif self.symmetrization == 'sum-type':
                beta4 = self._get_beta_tns_sum_type(phi_mtx)
                
            # apparent strain tensor:
            eps_app_mtx = map2d_eps_eng_to_mtx(eps_app_eng)
            # effective strain tensor:
            eps_eff_mtx = tensordot(beta4, eps_app_mtx, [[0,1],[0,1]])
            # effective stress tensor:
            sig_eff_mtx = tensordot( self.D4_e, eps_eff_mtx, [[2,3],[0,1]])           
            # effective microplane stresses obtained by projection (static constraint)
            s_eff_vct_arr = array( [ dot( sig_eff_mtx, mpn ) for mpn in self._MPN ] )
            # microplane scalar damage variable (integrity):
            phi_arr = self._get_phi_arr( sctx, eps_app_eng )
            # apparent microplane stresses 
            s_app_vct_arr = array( [ dot( phi, s_eff_vct ) for phi, s_eff_vct in zip(phi_arr,s_eff_vct_arr) ] )

        #--------------------
        # compliance version:
        #--------------------
        if self.model_version == 'compliance' :
            # get the corresponding macroscopic stresses
            sig_app_eng, D_mtx  = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )

            ### microplane equivalent stress obtained by projection (static constraint)
            s_app_vct_arr   = self._get_s_vct_arr(sig_app_eng)

            ### microplane equivalent strains calculated based on corresponding 'M' and 'psi_mtx'
            # 2nd order damage effect tensor:
            psi_mtx = self._get_psi_mtx(sctx, eps_app_eng)
            # 4th order damage effect tensor:
            if self.symmetrization == 'product-type':
                M4 = self._get_M_tns_product_type(psi_mtx)
            elif self.symmetrization == 'sum-type':
                M4 = self._get_M_tns_sum_type(psi_mtx)

            # apparent stress tensor:
            sig_app_mtx = map2d_sig_eng_to_mtx(sig_app_eng)
            # effective stress tensor:
            sig_eff_mtx = tensordot(M4, sig_app_mtx, [[2,3],[0,1]])
            # effective strain tensor:
            eps_eff_mtx = tensordot(self.C4_e, sig_eff_mtx, [[2,3],[0,1]])
            # effective microplane strains obtained by projection (kinematic constraint)
            e_eff_vct_arr = array( [ dot( eps_eff_mtx, mpn ) for mpn in self._MPN ] )
            # microplane scalar damage variable (integrity):
            phi_arr = self._get_phi_arr( sctx, eps_app_eng )
            psi_arr = 1./phi_arr
            # apparent microplane strains 
            e_app_vct_arr = array( [ dot( psi, e_eff_vct ) for psi, e_eff_vct in zip(psi_arr,e_eff_vct_arr) ] )

        return e_app_vct_arr, s_app_vct_arr


    # Equiv:
    def _get_e_s_equiv_arr(self, sctx, eps_app_eng):
        '''
        Returns two arrays containing the corresponding equivalent
        microplane strains and stresses.
        '''                                         
        e_app_vct_arr, s_app_vct_arr = self._get_e_s_vct_arr( sctx, eps_app_eng )
        
        # equivalent strain for each microplane
        e_equiv_arr = self._get_e_equiv_arr( e_app_vct_arr )
        # equivalent strain for each microplane
        s_equiv_arr = self._get_s_equiv_arr( s_app_vct_arr )
        return e_equiv_arr, s_equiv_arr

    def get_e_equiv_arr(self, sctx, eps_app_eng):
        '''
        Return a list of equivalent microplane strains.
        '''
        e_equiv_arr, s_equiv_arr = self._get_e_s_equiv_arr(sctx, eps_app_eng)
        return e_equiv_arr
    
    def get_s_equiv_arr(self, sctx, eps_app_eng):
        '''
        Return a list of equivalent microplane stresses.
        '''
        e_equiv_arr, s_equiv_arr = self._get_e_s_equiv_arr(sctx, eps_app_eng)
        return s_equiv_arr 
    

    # N: normal microplane strain and stresses:
    def _get_e_s_N_arr(self, sctx, eps_app_eng):
        '''
        Returns two arrays containing the corresponding normal 
        microplane strains and stress components.
        '''                                         
        e_app_vct_arr, s_app_vct_arr = self._get_e_s_vct_arr( sctx, eps_app_eng )
        # normal strain for each microplane
        e_N_arr = self._get_e_N_arr( e_app_vct_arr )
        # normal strain for each microplane
        s_N_arr = self._get_s_N_arr( s_app_vct_arr )
        return e_N_arr, s_N_arr

    def get_e_N_arr(self, sctx, eps_app_eng):
        '''
        Return a list of normal microplane strains components. 
        '''
        e_N_arr, s_N_arr = self._get_e_s_N_arr(sctx, eps_app_eng)
        return e_N_arr
    
    def get_s_N_arr(self, sctx, eps_app_eng):
        '''
        Return a list of normal microplane stresses components. 
        '''
        e_N_arr, s_N_arr = self._get_e_s_N_arr(sctx, eps_app_eng)
        return s_N_arr 


    # T: tangential microplane strain and stresses:
    def _get_e_s_T_arr(self, sctx, eps_app_eng):
        '''
        Returns two arrays containing the corresponding tangential 
        microplane strains and stress components. 
        '''                                         
        e_app_vct_arr, s_app_vct_arr = self._get_e_s_vct_arr( sctx, eps_app_eng )
        # shear strain for each microplane
        e_T_arr = self._get_e_T_arr( e_app_vct_arr )
        # shear strain for each microplane
        s_T_arr = self._get_s_T_arr( s_app_vct_arr )
        return e_T_arr, s_T_arr

    def get_e_T_arr(self, sctx, eps_app_eng):
        '''
        Return a list of tangential microplane strains components. 
        '''
        e_T_arr, s_T_arr = self._get_e_s_T_arr(sctx, eps_app_eng)
        return e_T_arr
    
    def get_s_T_arr(self, sctx, eps_app_eng):
        '''
        Return a list of tangential microplane stresses components. 
        '''
        e_T_arr, s_T_arr = self._get_e_s_T_arr(sctx, eps_app_eng)
        return s_T_arr 

### @todo: remove / old

# @todo: remove: this method and method 'get_integ' in 'polar_fn'
### old implementation: this assumes a decoupled reaction of all microplanes
    def get_fracture_energy_arr(self, sctx, eps_app_eng):
        '''
        Get the microplane contributions to the fracture energy
        '''
        e_max_arr = self._get_state_variables(sctx, eps_app_eng)
        fracture_energy_arr = self.polar_fn.get_fracture_energy_arr( sctx, e_max_arr )
        return fracture_energy_arr

    def get_fracture_energy(self, sctx, eps_app_eng):
        '''
        Get the macroscopic fracture energy as a weighted sum of all mircoplane contributions
        '''
        e_max_arr = self._get_state_variables(sctx, eps_app_eng)
        fracture_energy_arr = self.polar_fn.get_fracture_energy_arr( sctx, e_max_arr )
        return array( [dot( self._MPW, fracture_energy_arr )], float )
    

    def get_e_equiv_projection(self, sctx, eps_app_eng):
        '''
        Return a list of equivalent microplane strains as the projection of the strain tensor
        (kinematic constraint)
        '''
        e_vct_arr = self._get_e_vct_arr(eps_app_eng)
        e_equiv_arr = self._get_e_equiv_arr(e_vct_arr)
#        print 'e_equiv_arr:  ', e_equiv_arr
        return e_equiv_arr

###

    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return { 'eps_app'                  : self.get_eps_app,
                 'sig_app'                  : self.get_sig_app,
                 'sig_norm'                 : self.get_sig_norm,
                 'phi_mtx'                  : self.get_phi_mtx,
                 'phi_pdc'                  : self.get_phi_pdc,
                 'microplane_damage'        : RTraceEval( eval = self.get_microplane_integrity,
                                                          ts = self ),
                                                          
                 'e_equiv_arr' : RTraceEval( eval = self.get_e_equiv_arr,
                                             ts = self ),
                 's_equiv_arr' : RTraceEval( eval = self.get_s_equiv_arr,
                                             ts = self ),
                 'e_N_arr'     : RTraceEval( eval = self.get_e_N_arr,
                                             ts = self ),
                 's_N_arr'     : RTraceEval( eval = self.get_s_N_arr,
                                             ts = self ),
                 'e_T_arr'     : RTraceEval( eval = self.get_e_T_arr,
                                             ts = self ),
                 's_T_arr'     : RTraceEval( eval = self.get_s_T_arr,
                                             ts = self ),

                 'equiv_projection'         : RTraceEval( eval = self.get_e_equiv_projection,
                                                          ts = self ),
                 'fracture_energy_arr'      : self.get_fracture_energy_arr,
                 'fracture_energy'          : self.get_fracture_energy 
                 }

# @todo: - temporary alias rename the class and test it all

MA2DCompositeMicroplaneDamage = MA2DMicroplaneDamage


if __name__ == '__main__':
    #--------------------------------------------------------------------------------
    # Example using the mats2d_explore 
    #--------------------------------------------------------------------------------
    from ibvpy.mats.mats2D.mats2D_explore import MATS2DExplore
    from ibvpy.mats.mats2D.mats2D_rtrace_cylinder import MATS2DRTraceCylinder
    from mats2D_cmdm_rtrace_Gf_mic import MATS2DMicroplaneDamageTraceGfmic,\
        MATS2DMicroplaneDamageTraceEtmic, MATS2DMicroplaneDamageTraceUtmic
    from mats2D_cmdm_rtrace_Gf_mac import MATS2DMicroplaneDamageTraceGfmac,\
        MATS2DMicroplaneDamageTraceEtmac, MATS2DMicroplaneDamageTraceUtmac
    
    mats2D_explore = \
        MATS2DExplore( mats2D_eval = MA2DMicroplaneDamage( elastic_debug = False ),
                       rtrace_list = [ 
                                       MATS2DRTraceCylinder(name = 'Laterne',
                                                            var_axis    = 'time', idx_axis = 0,
                                                            var_surface = 'microplane_damage',
                                                            update_on = 'update' ),
                                                            
                                       RTraceGraph(name = 'strain 0 - stress 0',
                                                   var_x = 'eps_app', idx_x = 0,
                                                   var_y = 'sig_app', idx_y = 0,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'strain 0 - strain 1',
                                                   var_x = 'eps_app', idx_x = 0,
                                                   var_y = 'eps_app', idx_y = 1,
                                                   record_on = 'update' ),
                                       RTraceGraph(name = 'stress 0 - stress 1',
                                                   var_x = 'sig_app', idx_x = 0,
                                                   var_y = 'sig_app', idx_y = 1,
                                                   record_on = 'update' ),

                                       RTraceGraph(name = 'time - microplane damage',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'microplane_damage', idx_y = 0,
                                                   update_on = 'update' ),

                                       # e_equiv, s_equiv
                                       RTraceGraph(name = 'e_equiv - s_equiv',
                                                   var_x = 'e_equiv_arr', idx_x = 0,
                                                   var_y = 's_equiv_arr', idx_y = 0,
                                                   update_on = 'update' ),
                                    
                                       # e_N, s_N:            
                                       RTraceGraph(name = 'e_N - s_N',
                                                   var_x = 'e_N_arr', idx_x = 0,
                                                   var_y = 's_N_arr', idx_y = 0,
                                                   update_on = 'update' ),

                                       # e_T, s_T:            
                                       RTraceGraph(name = 'e_T - s_T',
                                                   var_x = 'e_T_arr', idx_x = 0,
                                                   var_y = 's_T_arr', idx_y = 0,
                                                   update_on = 'update' ),

                                       RTraceArraySnapshot(name = 'equiv_projection',
                                                           var = 'equiv_projection',
                                                           record_on = 'update' ),
                                                   
                                       RTraceArraySnapshot(name = 'microplane damage',
                                                           var = 'microplane_damage',
                                                           record_on = 'update' ),

                                       RTraceArraySnapshot(name = 'e_equiv',
                                                           var = 'e_equiv_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 's_equiv',
                                                           var = 's_equiv_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 'e_N',
                                                           var = 'e_N_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 's_N',
                                                           var = 's_N_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 'e_T',
                                                           var = 'e_T_arr',
                                                           record_on = 'update' ),
                                       RTraceArraySnapshot(name = 's_T',
                                                           var = 's_T_arr',
                                                           record_on = 'update' ),
                                                           
#                                       # G_f_mic: microplane fracture energy:
#                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_equiv',
#                                                                        var_x = 'e_equiv_arr', idx_x = 0,
#                                                                        var_y = 's_equiv_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_N',
#                                                                        var_x = 'e_N_arr', idx_x = 0,
#                                                                        var_y = 's_N_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceGfmic(name = 'G_f_mic_T',
#                                                                        var_x = 'e_T_arr', idx_x = 0,
#                                                                        var_y = 's_T_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       # E_t_mic: microplane total energy                                                                       
#                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_equiv',
#                                                                        var_x = 'e_equiv_arr', idx_x = 0,
#                                                                        var_y = 's_equiv_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_N',
#                                                                        var_x = 'e_N_arr', idx_x = 0,
#                                                                        var_y = 's_N_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceEtmic(name = 'E_t_mic_T',
#                                                                        var_x = 'e_T_arr', idx_x = 0,
#                                                                        var_y = 's_T_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       # U_t_mic: microplane elastic energy                                                                        
#                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_equiv',
#                                                                        var_x = 'e_equiv_arr', idx_x = 0,
#                                                                        var_y = 's_equiv_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_N',
#                                                                        var_x = 'e_N_arr', idx_x = 0,
#                                                                        var_y = 's_N_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
#                                       MATS2DMicroplaneDamageTraceUtmic(name = 'U_t_mic_T',
#                                                                        var_x = 'e_T_arr', idx_x = 0,
#                                                                        var_y = 's_T_arr', idx_y = 0,
#                                                                        record_on = 'update' ),
                                
                                       # direction 11:                  
                                       # G_f_mac: macroscopic fracture energy:
                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_11',

                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        update_on = 'update' ),
                                       # E_t_mac: macroscopic total energy:
                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        update_on = 'update' ),
                                       # U_t_mac: macroscopic elastic energy:
                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_11',
                                                                        var_x = 'eps_app', idx_x = 0,
                                                                        var_y = 'sig_app', idx_y = 0,
                                                                        update_on = 'update' ),

#                                       # direction 22:
#                                       # G_f_mac: macroscopic fracture energy:
#                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_22',
#                                                                        var_x = 'eps_app', idx_x = 1,
#                                                                        var_y = 'sig_app', idx_y = 1,
#                                                                        record_on = 'update' ),
#                                       # E_t_mac: macroscopic total energy:
#                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_22',
#                                                                        var_x = 'eps_app', idx_x = 1,
#                                                                        var_y = 'sig_app', idx_y = 1,
#                                                                        record_on = 'update' ),
#                                       # U_t_mac: macroscopic elastic energy:
#                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_22',
#                                                                        var_x = 'eps_app', idx_x = 1,
#                                                                        var_y = 'sig_app', idx_y = 1,
#                                                                        record_on = 'update' ),
#                                       
#                                       # direction 12:                                     
#                                       # G_f_mac: macroscopic fracture energy:
#                                       MATS2DMicroplaneDamageTraceGfmac(name = 'G_f_mac_12',
#                                                                        var_x = 'eps_app', idx_x = 2,
#                                                                        var_y = 'sig_app', idx_y = 2,
#                                                                        record_on = 'update' ),
#                                       # E_t_mac: macroscopic total energy:
#                                       MATS2DMicroplaneDamageTraceEtmac(name = 'E_t_mac_12',
#                                                                        var_x = 'eps_app', idx_x = 2,
#                                                                        var_y = 'sig_app', idx_y = 2,
#                                                                        record_on = 'update' ),
#                                       # U_t_mac: macroscopic elastic energy:
#                                       MATS2DMicroplaneDamageTraceUtmac(name = 'U_t_mac_12',
#                                                                        var_x = 'eps_app', idx_x = 2,
#                                                                        var_y = 'sig_app', idx_y = 2,
#                                                                        update_on = 'update' ),

#                                       RTraceArraySnapshot(name = 'fracture energy contributions',
#                                                           var = 'fracture_energy_arr',

#                                                           update_on = 'update' ),

##                                     ### decoupled energy contributions for G_f
#                                     RTraceGraph(name = 'time - G_f',
#                                                  var_x = 'time', idx_x = 0,
#                                                  var_y = 'fracture_energy', idx_y = 0,
#                                                  record_on = 'update' ),
#                                     ###
                                       RTraceGraph(name = 'time - sig_norm',
                                                   var_x = 'time', idx_x = 0,
                                                   var_y = 'sig_norm', idx_y = 0,
                                                   record_on = 'update' ),
#                                       RTraceGraph(name = 'time - phi_pdc',
#                                                   var_x = 'time', idx_x = 0,
#                                                   var_y = 'phi_pdc', idx_y = 0,
#                                                   record_on = 'update' ),
#                                       # e_equiv_projection:            
#                                       RTraceGraph(name = 'e_equiv_projection - s_equiv',
#                                                   var_x = 'equiv_projection', idx_x = 0,
#                                                   var_y = 's_equiv', idx_y = 0,
#                                                   record_on = 'update' ),

                                     ]
                       )

        

    # Elastic debug:    
    mats2D_explore.mats2D_eval.elastic_debug = False
    print 'elastic_debug: ', mats2D_explore.mats2D_eval.elastic_debug
    
    
    
    #---------------------------
    # Default settings:
    #---------------------------

    mats2D_explore.mats2D_eval.polar_fn.n_mp = 6
    print 'n_mp:          ', mats2D_explore.mats2D_eval.polar_fn.n_mp

    mats2D_explore.alpha_degree = 0.
    print 'alpha_degree:  ', mats2D_explore.alpha_degree

#    mats2D_explore.mats2D_eval.stress_state = 'plane_strain'
    mats2D_explore.mats2D_eval.stress_state = 'plane_stress'

#    mats2D_explore.mats2D_eval.symmetrization = 'product-type'
    mats2D_explore.mats2D_eval.symmetrization = 'sum-type'
    print 'symmetrization:', mats2D_explore.mats2D_eval.symmetrization

    mats2D_explore.mats2D_eval.model_version = 'compliance'
#   mats2D_explore.mats2D_eval.model_version = 'stiffness'
    print 'model_version: ', mats2D_explore.mats2D_eval.model_version

    # specify if reinforcement is used or not, i.e. isotropic, quasi-brittel case or
    # anisotropic, tension-stiffening case is considered.


    test_fitted_phi_fn = True
    #---------------------------
    # test_fitted_phi_fn
    #---------------------------
    if test_fitted_phi_fn:
        print 'test_fitted_phi_fn = True'
        # Default settings for 'polar_fn_class':
        mats2D_explore.mats2D_eval.polar_fn_class = 'Isotropic polar function'
        # Default settings of 'PolarFnBase' for 'phi_fn_class':
        mats2D_explore.mats2D_eval.polar_fn.phi_fn_class = 'General'
#        fitted_phi_fn = loadtxt(join( 'apps', 'cmdm2d_lab', 'mats_lab', 'target_data', 'experimental_data', 'fitted_phi_fn.out'))
        fitted_phi_fn = loadtxt( 'fitted_phi_fn.out' )
        print 'fitted_phi_fn', fitted_phi_fn
        x = fitted_phi_fn[0]
        y = fitted_phi_fn[1]
        print 'x', x
        print 'y', y
        mats2D_explore.mats2D_eval.polar_fn.phi_fn.mfn.set( xdata = x, ydata = y )
        mats2D_explore.mats2D_eval.polar_fn.phi_fn.mfn.data_changed = True



    reinforced = False    

    #---------------------------
    # Isotropic polar function
    #---------------------------
    if reinforced == False and test_fitted_phi_fn == False:
        # Default settings for 'polar_fn_class':
        mats2D_explore.mats2D_eval.polar_fn_class = 'Isotropic polar function'
        # Default settings of 'PolarFnBase' for 'phi_fn_class':
        mats2D_explore.mats2D_eval.polar_fn.phi_fn_class = 'QuasiBrittle'
        print'aaa'
    #---------------------------
    # Anisotropic polar function
    #---------------------------
    elif reinforced == True and test_fitted_phi_fn == False:
        mats2D_explore.mats2D_eval.polar_fn_class = 'Anisotropic polar function'
    
        # Only available settings of 'AnisotropicPolarFn' for 'phi_fn_class':
        mats2D_explore.mats2D_eval.polar_fn.phi_fn_class = 'TensionStiffening'
        
        # default parameters of 'phi_fn':
        mats2D_explore.mats2D_eval.polar_fn.phi_fn.Epp =  59e-6 # 0.5e-3
        mats2D_explore.mats2D_eval.polar_fn.phi_fn.Efp = 191e-6 # 4.0e-3
        mats2D_explore.mats2D_eval.polar_fn.phi_fn.Dfp = 0.4
        mats2D_explore.mats2D_eval.polar_fn.phi_fn.Elimit = 1.e-3 
        
        # Polar function: parameter definition: 
        mats2D_explore.mats2D_eval.polar_fn.varied_params = ['Dfp']
        mats2D_explore.mats2D_eval.polar_fn.varpars['Dfp'].polar_fn.set( \
                                                             phi_residual = 1.,
                                                             phi_quasibrittle = 0.0,
                                                             alpha = 0.,
                                                             delta_trans = Pi/4. ,
                                                             delta_alpha = Pi/8. ) 


#    mats2D_explore.mats2D_eval.configure_traits()
    mats2D_explore.tloop.eval()
    mats2D_explore.tloop.setup()
    mats2D_explore.tloop.rtrace_mngr.rtrace_bound_list[1].redraw()

    print 'xdata: ', mats2D_explore.tloop.rtrace_mngr.rtrace_bound_list[1].trace.xdata
    print 'ydata: ', mats2D_explore.tloop.rtrace_mngr.rtrace_bound_list[1].trace.ydata
    
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = mats2D_explore )
    ibvpy_app.main()
    
#-------------------------------
# @todo's: 
#-------------------------------
#
#1)Renaming:
#  - 'MA2DMicroplaneDamage' to 'MATS2DCompositeMicroplaneDamage'
#  Check if further renaming would be helpful (propositions):
#  - 'PolarFn' to something like 'PolarDiscr' to avoid ambiguity with 'polar_fn = MfnPolar'
#  - filename 'cmdm_phi_fn' to 'cmdm_phi_fn' and 'cmdm_polar_discr' to 'cmdm_polar_discr'
#  - all vectorized versions with suffix '_vectorized'  
#
#2) Check if caching of n_mp is not working properly. It seams to work now. 
#   What was the problem?
#
#3) Task done.
#
#4) Note: parameters specified above correspond to the parameters for paper CST2008.
#   There is a deviation in the value for the stresses under tension stiffening for eps = 2E-3
#   it's about 7 MPa instead of 12 as is the case in the paper. Choices for 'stress_state' 
#   and 'symmetrization' are not relevant in that context. What is the reason? 

    
    