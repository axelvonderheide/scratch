from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

# Chaco imports
from enthought.chaco.chaco_plot_container_editor import \
     PlotContainerEditor
from enthought.chaco.tools.api import \
     PanTool, SimpleZoom
from enthought.chaco.api import \
     Plot, AbstractPlotData, ArrayPlotData

from numpy import \
     array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt, swapaxes

from math import \
    pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from scipy.linalg import \
    eig, inv


from ibvpy.core.tstepper import \
    TStepper as TS

from ibvpy.mats.mats_eval import \
    IMATSEval, MATSEval

from ibvpy.core.rtrace_eval import \
    RTraceEval

from ibvpy.api import \
    RTraceGraph, RTraceArraySnapshot

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

class MATS2DMicroplaneDamage( MATSEval ):
    '''
    Microplane Damage Model.
    '''
    implements( IMATSEval )

    #------------------------------------------------------------------------------------
    # Damage function specification
    #------------------------------------------------------------------------------------
    polar_fn_class = Trait( 'Isotropic polar function',
                            {'Isotropic polar function'   : IsotropicPolarFn,
                             'Anisotropic polar function' : AnisotropicPolarFn } )

    polar_fn = Property( Instance( IPolarFn ), depends_on = 'polar_fn_class' )
    @cached_property
    def _get_polar_fn( self ):
        return self.polar_fn_class_()

    def _set_polar_fn( self, value ):
        if not ( self.polar_fn_class_ is value.__class__ ):
            raise ValueError, 'class mismatch %s != %s' % ( self.polar_fn_class_, value.__class__ )
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

    n_mp = Delegate( 'polar_fn' )

    model_version = Enum( "compliance", "stiffness" )
    stress_state = Enum( "plane_strain", "plane_stress" )
    symmetrization = Enum( "product-type", "sum-type" )

    dimensionality = Trait( '2D', {'2D' : 2, '3D' : 3 },
                 label = 'Dimensionality',
                 desc = 'This value is inactive yet',
                 auto_set = False )

    @on_trait_change( 'n_mp, dimensionality' )
    def update_arrays( self ):
        self._init_arrays()
        self.changed = True

    elastic_debug = Bool( False,
                           desc = 'switch to elastic behavior - used for debugging',
                           auto_set = False )

    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View( VSplit( Group( Item( 'polar_fn_class', style = 'custom', show_label = False ),
                                      Item( 'polar_fn', style = 'custom', show_label = False ),
                                      label = 'Material parameters',
                                      show_border = True ),
                                Group( Item( 'model_version', style = 'custom' ),
                                       Item( 'stress_state', style = 'custom' ),
                                       Item( 'symmetrization', style = 'custom' ),
                                       Item( 'dimensionality', style = 'custom' ),
                                       Item( 'elastic_debug@' ),
                                       Spring( resizable = True ),
                                       label = 'Configuration parameters', show_border = True,
                                       ),
                                ),
                        resizable = True
                        )

    #-----------------------------------------------------------------------------------------------
    # Private initialization methods
    #-----------------------------------------------------------------------------------------------

    def __init__( self, **kwtraits ):
        '''
        Subsidiary arrays required for the integration.
        they are set up only once prior to the computation
        '''
        super( MATS2DMicroplaneDamage, self ).__init__( **kwtraits )
        self._init_arrays()

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
            self._MPN = array( [[ cos( alpha ), sin( alpha )] for alpha in self.alpha_list ] )
            # get the weights of the microplanes
            self._MPW = ones( self.n_mp ) / self.n_mp * 2
            # get the dyadic product of the microplane normals
            self._MPNN = array( [ outer( mpn, mpn ) for mpn in self._MPN ] )
        # 3D-case:
        else:
            # Use the numerical integration formula for the n-sphere by STROUD (1979)
            # with a number of microplanes equal to 28:
            self.n_mp = 28
            # microplane normals:
            self._MPN = array( [[.577350259, .577350259, .577350259], \
                         [.577350259, .577350259, -.577350259], \
                         [.577350259, -.577350259, .577350259], \
                         [.577350259, -.577350259, -.577350259], \
                         [.935113132, .250562787, .250562787], \
                         [.935113132, .250562787, -.250562787], \
                         [.935113132, -.250562787, .250562787], \
                         [.935113132, -.250562787, -.250562787], \
                         [.250562787, .935113132, .250562787], \
                         [.250562787, .935113132, -.250562787], \
                         [.250562787, -.935113132, .250562787], \
                         [.250562787, -.935113132, -.250562787], \
                         [.250562787, .250562787, .935113132], \
                         [.250562787, .250562787, -.935113132], \
                         [.250562787, -.250562787, .935113132], \
                         [.250562787, -.250562787, -.935113132], \
                         [.186156720, .694746614, .694746614], \
                         [.186156720, .694746614, -.694746614], \
                         [.186156720, -.694746614, .694746614], \
                         [.186156720, -.694746614, -.694746614], \
                         [.694746614, .186156720, .694746614], \
                         [.694746614, .186156720, -.694746614], \
                         [.694746614, -.186156720, .694746614], \
                         [.694746614, -.186156720, -.694746614], \
                         [.694746614, .694746614, .186156720], \
                         [.694746614, .694746614, -.186156720], \
                         [.694746614, -.694746614, .186156720], \
                         [.694746614, -.694746614, -.186156720]] )
            # microplane weights:
            # (values in the array must be multiplied by 6 (cf. [Baz05]))
            self._MPW = array( [.0160714276, .0160714276, .0160714276, .0160714276, .0204744730, \
                         .0204744730, .0204744730, .0204744730, .0204744730, .0204744730, \
                         .0204744730, .0204744730, .0204744730, .0204744730, .0204744730, \
                         .0204744730, .0158350505, .0158350505, .0158350505, .0158350505, \
                         .0158350505, .0158350505, .0158350505, .0158350505, .0158350505, \
                         .0158350505, .0158350505, .0158350505 ] ) * 6.0
            # dyadic product of the microplane normals
            self._MPNN = array( [ outer( mpn, mpn ) for mpn in self._MPN ] )


    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------

    def get_state_array_size():
        # In the state array the largest equivalent microplane 
        # strains reached in the loading history are saved
        return self.n_mp

    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        state_arr_size = self.get_state_array_size()
        sctx.mats_state_array = zeros( state_arr_size, 'float_' )
        #
        self.update_state_on = False
        self._setup_elasticity_tensors()


    def _check_dimensionality( self ):
        if self.dimensionality == '2D':
            print 'dimensionality: 2D-case'
            self.n_dim = 2
            self.n_eng = 3
        elif self.dimensionality == '3D':
            print 'dimensionality: 3D-case'
            self.n_dim = 3
            self.n_eng = 6


    def _setup_elasticity_tensors( self ):
        '''
        Intialize the fourth order elasticity tensor for 3D or 2D plane strain or 2D plane stress
        '''
        # ----------------------------------------------------------------------------
        # Lame constants calculated from E and nu
        # ----------------------------------------------------------------------------
        E = self.E
        nu = self.nu
        # first Lame paramter
        la = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
        # second Lame parameter (shear modulus)
        mu = E / ( 2 + 2 * nu )

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 3D-case
        # -----------------------------------------------------------------------------------------------------
        D4_e_3D = zeros( [3, 3, 3, 3] )
        C4_e_3D = zeros( [3, 3, 3, 3] )
        delta = identity( 3 )
        for i in range( 0, 3 ):
            for j in range( 0, 3 ):
                for k in range( 0, 3 ):
                    for l in range( 0, 3 ):
                        D4_e_3D[i, j, k, l] = la * delta[i, j] * delta[k, l] + \
                                            mu * ( delta[i, k] * delta[j, l] + delta[i, l] * delta[j, k] )
                        # elastic compliance tensor:                 
                        C4_e_3D[i, j, k, l] = ( 1 + nu ) / ( 2 * E ) * \
                                            ( delta[i, k] * delta[j, l] + delta[i, l] * delta[j, k] ) - \
                                              nu / E * delta[i, j] * delta[k, l]

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 2D-case
        # -----------------------------------------------------------------------------------------------------
        # 1. step: Get the (6x6)-elasticity and compliance matrices 
        #          for the 3D-case:
        D2_e_3D = map3d_tns4_to_tns2( D4_e_3D )
        C2_e_3D = map3d_tns4_to_tns2( C4_e_3D )

        # 2. step: Get the (3x3)-elasticity and compliance matrices 
        #          for the 2D-cases plane stress and plane strain:
        D2_e_2D_plane_stress = get_D_plane_stress( D2_e_3D )
        D2_e_2D_plane_strain = get_D_plane_strain( D2_e_3D )
        C2_e_2D_plane_stress = get_C_plane_stress( C2_e_3D )
        C2_e_2D_plane_strain = get_C_plane_strain( C2_e_3D )

        # 3. step: Get the fourth order elasticity and compliance tensors
        #          for the 2D-cases plane stress and plane strain (D4.shape = (2,2,2,2)) 
        D4_e_2D_plane_stress = map2d_tns2_to_tns4( D2_e_2D_plane_stress )
        D4_e_2D_plane_strain = map2d_tns2_to_tns4( D2_e_2D_plane_strain )
        C4_e_2D_plane_stress = map2d_tns2_to_tns4( C2_e_2D_plane_stress )
        C4_e_2D_plane_strain = map2d_tns2_to_tns4( C2_e_2D_plane_strain )

        # -----------------------------------------------------------------------------------------------------
        # assign the fourth order elasticity and compliance tensors to self
        # -----------------------------------------------------------------------------------------------------
        # (3D case):                                       
        if self.dimensionality == '3D':
            self.D4_e = D4_e_3D
            self.C4_e = C4_e_3D

        # (2D case):
        if self.dimensionality == '2D' and self.stress_state == 'plane_stress':
            print 'stress state:   plane-stress'
            self.D4_e = D4_e_2D_plane_stress
            self.C4_e = C4_e_2D_plane_stress

        if self.dimensionality == '2D' and self.stress_state == 'plane_strain':
            print 'stress state:   plane-strain'
            self.D4_e = D4_e_2D_plane_strain
            self.C4_e = C4_e_2D_plane_strain

        # -----------------------------------------------------------------------------------------------------
        # assign the 3x3 - elasticity and compliance matrix to self (used only for elastic debug)
        # -----------------------------------------------------------------------------------------------------
        if self.elastic_debug == True:
            print 'elastic_debug:  switched on!'
            # (3D case):                                       
            if self.dimensionality == '3D':
                self.D2_e = D2_e_3D

            # (2D case):
            if self.dimensionality == '2D' and self.stress_state == 'plane_stress':
                self.D2_e = D2_e_2D_plane_stress

            if self.dimensionality == '2D' and self.stress_state == 'plane_strain':
                self.D2_e = D2_e_2D_plane_strain

        # -----------------------------------------------------------------------------------------------------
        # assign the identity matrix to self (used only for sum-type symmetrization)
        # -----------------------------------------------------------------------------------------------------
        if self.symmetrization == 'sum-type':
            # (2D case):                                       
            if self.dimensionality == '2D':
                self.identity_tns = identity( 2 )
            else:
                self.identity_tns = identity( 3 )


    #-----------------------------------------------------------------------------------------------------
    # Prepare evaluation - get the 4th order damage tensor beta4 (or damage effect tensor M4)
    #-----------------------------------------------------------------------------------------------------

    def _get_e_vct_list( self, eps_eng, MPN = None ):
        '''
        Projects the strain tensor on the microplanes and
        returns a list of microplane strain vectors
         
        '''
        # Switch from engineering notation to tensor notation for the apparent strains
        if self.n_dim == 2:
            # 2D-case:
            eps_mtx = map2d_eps_eng_to_mtx( eps_eng )
        else:
            # 3D-case:
            eps_mtx = map3d_eps_mtx_to_eng( eps_eng )
        if MPN == None:
            MPN = self._MPN
        # Projection of apparent strain onto the individual microplanes
        e_vct_list = array( [ dot( eps_mtx, mpn ) for mpn in MPN ] )
        return e_vct_list


    def _get_e_equiv_list( self, e_vct_list ):
        '''
        Returns a list of the microplane equivalent strains
        '''
        # magnitude of the normal strain vector for each microplane
        e_N_list = array( [ dot( e_vct, mpn ) for e_vct, mpn in zip( e_vct_list, self._MPN ) ] )
        # positive part of the normal strain magnitude for each microplane
        e_N_pos_list = ( fabs( e_N_list ) + e_N_list ) / 2
        # normal strain vector for each microplane
        e_N_vct_list = array( [ self._MPN[i, :] * e_N_list[i] for i in range( 0, self.n_mp ) ] )
        # tangent strain ratio
        c_T = self.c_T
        # tangential strain vector for each microplane
        e_T_vct_list = e_vct_list - e_N_vct_list
        # squared tangential strain vector for each microplane
        e_TT_list = array( [ inner( e_T_vct, e_T_vct ) for e_T_vct in e_T_vct_list ] )
        # equivalent strain for each microplane
        e_equiv_list = sqrt( e_N_pos_list * e_N_pos_list + c_T * e_TT_list )
        return e_equiv_list


    def _get_e_max( self, e_equiv, e_max ):
        '''
        Compares the equivalent microplane strain of a single microplane with the 
        maximum strain reached in the loading history
        '''
        if e_equiv >= e_max:
            e_max = e_equiv
        return e_max


    def _get_state_variables( self, sctx, eps_app_eng ):
        '''
        Compares the list of current equivalent microplane strains with 
        the values in the state array and returns the higher values 
        '''
        get_e_max_vectorized = frompyfunc( self._get_e_max, 2, 1 )
        e_vct_list = self._get_e_vct_list( eps_app_eng )
        e_equiv_list = self._get_e_equiv_list( e_vct_list )
        e_max_list_old = sctx.mats_state_array
        e_max_list_new = get_e_max_vectorized( e_equiv_list, e_max_list_old )
        return e_max_list_new


    def _get_phi_list( self, sctx, eps_app_eng ):
        '''
        Returns a list of the integrity factors for all microplanes.
        '''
        e_max_list = self._get_state_variables( sctx, eps_app_eng )
        return self.polar_fn.get_phi_list( e_max_list )


    def _get_phi_mtx( self, sctx, eps_app_eng ):
        '''
        Returns the 2nd order damage tensor 'phi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_list = self._get_phi_list( sctx, eps_app_eng )
        # integration terms for each microplanes
        phi_vct_list = array( [ phi_list[i] * self._MPNN[i, :, :] * self._MPW[i]
                                for i in range( 0, self.n_mp ) ] )
        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        phi_mtx = phi_vct_list.sum( 0 )
        return phi_mtx


    def _get_psi_mtx( self, sctx, eps_app_eng ):
        '''
        Returns the 2nd order damage effect tensor 'psi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_list = self._get_phi_list( sctx, eps_app_eng )
        # integration terms for each microplanes
        psi_vct_list = array( [ 1. / phi_list[i] * self._MPNN[i, :, :] * self._MPW[i]
                                for i in range( 0, self.n_mp ) ] )
        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        psi_mtx = psi_vct_list.sum( 0 )
        return psi_mtx


    def _get_beta_tns_product_type( self, phi_mtx ):
        '''
        Returns the 4th order damage tensor 'beta4' using product-type symmetrization
        '''
        # (cf. [Baz97], Eq.(87))
        n_dim = self.n_dim
        # Get the direction of the principle damage coordinates (pdc):
        phi_eig_value, phi_eig_mtx = eig( phi_mtx )
        phi_eig_value_real = array( [ pe.real for pe in phi_eig_value] )
        # transform phi_mtx to PDC:
        # (assure that besides the diagonal the entries are exactly zero)
        phi_pdc_mtx = zeros( [n_dim, n_dim] )
        for i in range( n_dim ):
            phi_pdc_mtx[i, i] = phi_eig_value_real[i]
        # w_mtx = tensorial square root of the second order damage tensor:
        w_pdc_mtx = arr_sqrt( phi_pdc_mtx )
        # transform the matrix w back to x-y-coordinates:
        w_mtx = dot( dot( phi_eig_mtx, w_pdc_mtx ), transpose( phi_eig_mtx ) )
        # beta_ijkl = w_ik * w_jl (cf. [Baz 97])
        # exploiting numpy-functionality (faster).
        # Method 'outer' returns beta_ijkl = w_ij * w_kl,
        # therefore the axis 2 and 3 need to be swapped
        beta4_ = outer( w_mtx, w_mtx ).reshape( n_dim, n_dim, n_dim, n_dim )
        beta4 = beta4_.swapaxes( 1, 2 )
        return beta4


    def _get_beta_tns_sum_type( self, phi_mtx ):
        '''
        Returns the 4th order damage tensor 'beta4' using sum-type symmetrization
        '''
        # (cf. [Jir99], Eq.(21))
        n_dim = self.n_dim
        delta = self.identity_tns
        beta4 = zeros( [n_dim, n_dim, n_dim, n_dim] )
        for i in range( 0, n_dim ):
            for j in range( 0, n_dim ):
                for k in range( 0, n_dim ):
                    for l in range( 0, n_dim ):
                        beta4[i, j, k, l] = 0.25 * ( phi_mtx[i, k] * delta[j, l] + phi_mtx[i, l] * delta[j, k] + \
                                                  phi_mtx[j, k] * delta[i, l] + phi_mtx[j, l] * delta[i, k] )
        return beta4


    def _get_M_tns_product_type( self, psi_mtx ):
        '''
        Returns the 4th order damage effect tensor 'M4' using product-type symmetrization
        '''
        n_dim = self.n_dim
        # Get the direction orthogonal to the principle damage coordinates (pdc):
        # @todo: is this direction orthogonal? Which one do we want?
        psi_eig_value, psi_eig_mtx = eig( psi_mtx )
        psi_eig_value_real = array( [ pe.real for pe in psi_eig_value] )
        # transform phi_mtx to PDC:
        # (assure that besides the diagonal the entries are exactly zero)
        psi_pdc_mtx = zeros( [n_dim, n_dim] )
        for i in range( n_dim ):
            psi_pdc_mtx[i, i] = psi_eig_value_real[i]
        # second order damage effect tensor:
        w_hat_pdc_mtx = arr_sqrt( psi_pdc_mtx )
        # transform the matrix w back to x-y-coordinates:
        w_hat_mtx = dot( dot( psi_eig_mtx, w_hat_pdc_mtx ), transpose( psi_eig_mtx ) )
        # M_ijkl = w_hat_ik * w_hat_lj (cf. Eq.(5.62) Script Prag Jirasek (2007))
        #        = w_hat_ik * w_hat_jl (w is a symmetric tensor) 
        # Exploiting numpy-functionality using the
        # method 'outer' (returns M_ijkl = w_hat_ij * w_hat_kl),
        # therefore the axis 2 and 3 need to be swapped
        M4_ = outer( w_hat_mtx, w_hat_mtx ).reshape( n_dim, n_dim, n_dim, n_dim )
        M4 = M4_.swapaxes( 1, 2 )
        return M4


    def _get_M_tns_sum_type( self, psi_mtx ):
        '''
        Returns the 4th order damage effect tensor 'M4' using sum-type symmetrization
        '''
        n_dim = self.n_dim
        delta = self.identity_tns
        # (cf. [Jir99], Eq.(30))
        M4 = zeros( [n_dim, n_dim, n_dim, n_dim] )
        for i in range( 0, n_dim ):
            for j in range( 0, n_dim ):
                for k in range( 0, n_dim ):
                    for l in range( 0, n_dim ):
                        M4[i, j, k, l] = 0.25 * ( psi_mtx[i, k] * delta[j, l] + psi_mtx[i, l] * delta[j, k] + \
                                               psi_mtx[j, k] * delta[i, l] + psi_mtx[j, l] * delta[i, k] )
        return M4


    #-----------------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------------
    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1, eps_avg = None ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
#        print 'eps_mtx000: get_corr_pred =', map2d_eps_eng_to_mtx(eps_app_eng)
        # -----------------------------------------------------------------------------------------------
        # for debugging purposes only: if elastic_debug is switched on, linear elastic material is used 
        # -----------------------------------------------------------------------------------------------
        if self.elastic_debug:
            D2_e = self.D2_e
            sig_eng = tensordot( D2_e, eps_app_eng, [[1], [0]] )
            return sig_eng, D2_e

        # -----------------------------------------------------------------------------------------------
        # update state variables
        # -----------------------------------------------------------------------------------------------
        if sctx.update_state_on:
            eps_n = eps_app_eng - d_eps
            e_max = self._get_state_variables( sctx, eps_n )
            sctx.mats_state_array[:] = e_max

        #------------------------------------------------------------------------------------------------
        # stiffness version:
        #------------------------------------------------------------------------------------------------
        if self.model_version == 'stiffness' :

            #------------------------------------------------------------------------------------
            # Damage tensor (2th order):  
            #------------------------------------------------------------------------------------
            phi_mtx = self._get_phi_mtx( sctx, eps_app_eng )

            #------------------------------------------------------------------------------------
            # Damage tensor (4th order) using product-type symmetrization: 
            #------------------------------------------------------------------------------------
            if self.symmetrization == 'product-type':
                # (cf. [Baz97])
                beta4 = self._get_beta_tns_product_type( phi_mtx )

            #-------------------------------------------------------------------------------------
            # Damage tensor (4th order) using sum-type symmetrization: 
            #-------------------------------------------------------------------------------------
            elif self.symmetrization == 'sum-type':
                # (cf. [Jir99] Eq.(21))
                beta4 = self._get_beta_tns_sum_type( phi_mtx )

            #-------------------------------------------------------------------------------------
            # Damaged stiffness tensor calculated based on the damage tensor beta4:
            #-------------------------------------------------------------------------------------
            # (cf. [Jir99] Eq.(7): C = beta * D_e * beta^T)
            D4_mdm = tensordot( beta4, tensordot( self.D4_e, beta4, [[2, 3], [2, 3]] ), [[2, 3], [0, 1]] )

            #-------------------------------------------------------------------------------------
            # Reduction of the fourth order tensor to a matrix assuming minor and major symmetry:
            #-------------------------------------------------------------------------------------
            # 2D-case:
            if self.n_dim == 2:
                D2_mdm = map2d_tns4_to_tns2( D4_mdm )
            # 3D-case:
            else:
                D4_mdm = map3d_tns4_to_tns2( D4_mdm )


        #------------------------------------------------------------------------------------------------
        # compliance version:
        #------------------------------------------------------------------------------------------------
        elif self.model_version == 'compliance' :

            #------------------------------------------------------------------------------------
            # Damage effect tensor (2th order):  
            #------------------------------------------------------------------------------------
            psi_mtx = self._get_psi_mtx( sctx, eps_app_eng )

            #------------------------------------------------------------------------------------
            # Damage effect tensor (4th order) using product-type-symmetrization:  
            #------------------------------------------------------------------------------------
            if self.symmetrization == 'product-type':
                M4 = self._get_M_tns_product_type( psi_mtx )

            #------------------------------------------------------------------------------------
            # Damage effect tensor (4th order) using sum-type-symmetrization:  
            #------------------------------------------------------------------------------------
            elif self.symmetrization == 'sum-type':
                # (cf. [Jir99], Eq.(30)
                M4 = self._get_M_tns_sum_type( psi_mtx )

            #-------------------------------------------------------------------------------------
            # Damaged compliance tensor calculated based on the damage effect tensor M4:
            #-------------------------------------------------------------------------------------
            # (cf. [Jir99] Eq.(8): C = M^T * C_e * M)
            C4_mdm = tensordot( M4, tensordot( self.C4_e, M4, [[2, 3], [0, 1]] ), [[0, 1], [0, 1]] )

            #-------------------------------------------------------------------------------------
            # Reduction of the fourth order tensor to a matrix assuming minor and major symmetry: 
            #-------------------------------------------------------------------------------------
            # 2D-case:
            if self.n_dim == 2:
                C2_mdm = map2d_tns4_to_tns2( C4_mdm )
                # multiply the resulting D matrix with the factors for compliance mapping
                D2_mdm = inv( compliance_mapping2d( C2_mdm ) )
            # 3D-case:
            else:
                C2_mdm = map3d_tns4_to_tns2( C4_mdm )
                # multiply the resulting D matrix with the factors for compliance mapping
                D2_mdm = inv( compliance_mapping3d( C2_mdm ) )

        #----------------------------------------------------------------------------------------
        # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
        #----------------------------------------------------------------------------------------
        sig_eng = tensordot( D2_mdm, eps_app_eng, [[1], [0]] )
        return sig_eng, D2_mdm


    #---------------------------------------------------------------------------------------------
    # Control variables and update state method
    #---------------------------------------------------------------------------------------------

    def new_cntl_var( self ):
        return zeros( self.n_eng, float_ )

    def new_resp_var( self ):
        return zeros( self.n_eng, float_ )

    def update_state( self, sctx, eps_app_eng ):
        '''
        Update state method is called upon an accepted time-step. 
        Here just set the flag on to make the update afterwards in the method itself.
        '''
        self.update_state_on = True


    #---------------------------------------------------------------------------------------------
    # Response trace evaluators
    #---------------------------------------------------------------------------------------------

    def get_eps_app( self, sctx, eps_app_eng ):
        return eps_app_eng

    def get_sig_app( self, sctx, eps_app_eng ):
        # @TODO
        # the stress calculation is performed twice - it might be
        # cached. But not in the spatial integration scheme.
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return sig_eng

    def get_sig_norm( self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return array( [ scalar_sqrt( sig_eng[0] ** 2 + sig_eng[1] ** 2 ) ] )

    def get_phi_pdc( self, sctx, eps_app_eng ):
        phi_mtx = self._get_phi_mtx( sctx, eps_app_eng )
        # Get the direction of the principle damage coordinates (pdc):
        phi_eig_value, phi_eig_mtx = eig( phi_mtx )
        phi_eig_value_real = array( [ pe.real for pe in phi_eig_value] )
        phi_pdc_mtx = zeros( [self.n_dim, self.n_dim] )
        for i in range( self.n_dim ):
            phi_pdc_mtx[i, i] = phi_eig_value_real[i]
        return array( [ phi_pdc_mtx[0, 0] ] )

    def get_microplane_damage( self, sctx, eps_app_eng ):
        phi_list = self._get_phi_list( sctx, eps_app_eng )
        return phi_list

    def _get_s_vct_list( self, sig_eng ):
        '''
        Projects the stress tensor onto the microplanes and
        returns a list of microplane strain vectors 
        '''
        # Switch from engineering notation to tensor notation for the apparent strains
        if self.n_dim == 2:
            # 2D-case:
            sig_mtx = map2d_sig_eng_to_mtx( sig_eng )
        else:
            # 3D-case:
            sig_mtx = map3d_sig_mtx_to_eng( sig_eng )
        # Projection of apparent strain onto the individual microplanes
        s_vct_list = array( [ dot( sig_mtx, mpn ) for mpn in self._MPN ] )
        return s_vct_list

    def _get_s_equiv_list( self, s_vct_list ):
        '''
        Get microplane equivalent stress
        '''
        # magnitude of the normal strain vector for each microplane
        s_N_list = array( [ dot( s_vct, mpn ) for s_vct, mpn in zip( s_vct_list, self._MPN ) ] )
        # positive part of the normal strain magnitude for each microplane
        s_N_pos_list = ( fabs( s_N_list ) + s_N_list ) / 2
        # normal strain vector for each microplane
        s_N_vct_list = array( [ self._MPN[i, :] * s_N_list[i] for i in range( 0, self.n_mp ) ] )
        # tangent strain ratio
        c_T = self.c_T
        # tangential strain vector for each microplane
        s_T_vct_list = s_vct_list - s_N_vct_list
        # squared tangential strain vector for each microplane
        s_TT_list = array( [ inner( s_T_vct, s_T_vct ) for s_T_vct in s_T_vct_list ] )
        # equivalent strain for each microplane
        s_equiv_list = sqrt( s_N_pos_list * s_N_pos_list + c_T * s_TT_list )
        return s_equiv_list

    def get_e_vct_list( self, sctx, eps_app_eng, MPN = None ):
        '''
        Return a list with current strain vectors on each microplane
        '''
        return self._get_e_vct_list( eps_app_eng, MPN )



    def get_e_s_equiv_list( self, sctx, eps_app_eng ):
        '''
        Return a list of equivalent microplane strains consistently derived based on the 
        specified model version, i.e. stiffness or compliance.
        '''
        #-----------------------
        # stiffness version:
        #-----------------------
        if self.model_version == 'stiffness' :
            ### microplane equivalent strains obtained by projection (kinematic constraint)
            e_vct_list = self._get_e_vct_list( eps_app_eng )
            e_equiv_list = self._get_e_equiv_list( e_vct_list )

            ### microplane equivalent stresses calculated based on corresponding 'beta' and 'phi_mtx'
            # 2nd order damage tensor:
            phi_mtx = self._get_phi_mtx( sctx, eps_app_eng )
            # 4th order damage tensor:
            if self.symmetrization == 'product-type':
                print 'stiffness + prod_type'
                beta4 = self._get_beta_tns_product_type( phi_mtx )
            elif self.symmetrization == 'sum-type':
                print 'stiffness + sum_type'
                beta4 = self._get_beta_tns_sum_type( phi_mtx )
            # apparent strain tensor:
            eps_app_mtx = map2d_eps_eng_to_mtx( eps_app_eng )

            # effective strain tensor:
            eps_eff_mtx = tensordot( beta4, eps_app_mtx, [[0, 1], [0, 1]] )
            # effective stress tensor:
            sig_eff_mtx = tensordot( self.D4_e, eps_eff_mtx, [[2, 3], [0, 1]] )

            # effective microplane stresses obtained by projection (static constraint)
            s_eff_vct_list = array( [ dot( sig_eff_mtx, mpn ) for mpn in self._MPN ] )
            # apparent microplane stresses 
            s_app_vct_list = array( [ dot( phi_mtx, s_eff_vct ) for s_eff_vct in s_eff_vct_list ] )
            # equivalent strain for each microplane
            s_equiv_list = self._get_s_equiv_list( s_app_vct_list )

#            # NOTE (1): the summation of the dyadic product of s_vct_list and n_vct yields the same stress tensor then get_corr_pred!
#            sig_app_eng, D = self.get_corr_pred(sctx, eps_app_eng, 0, 0, 0)
#            sig_app_mtx = map2d_sig_eng_to_mtx(sig_app_eng)
#            print 'sig_app_mtx0', sig_app_mtx
#            sig_app_mtx_list = array([outer(s_app_vct_list_i, MPN)*MPW for s_app_vct_list_i, MPN, MPW in zip(s_app_vct_list, self._MPN, self._MPW)])
#            sig_app_mtx = sig_app_mtx_list.sum(0)
#            if self.symmetrization == 'sum-type':
#                sig_app_mtx_sym = 0.5*(sig_app_mtx + transpose(sig_app_mtx))
#            elif self.symmetrization == 'product-type':
#                # @todo: verify if this in genneral is the correct way to symmetrize a matrix using product-type symmetrization
#                # Note that a problem arises here for negative values in eps_mtx.
#                # In phi_mtx only positiv values arise and therefore this way of 
#                # symmetrization is possible (cf. example [Jir99])  
#                sig_app_mtx_sym = sqrt(sig_app_mtx * transpose(sig_app_mtx))
#            print 'sig_app_mtx1', sig_app_mtx_sym

        #-----------------------
        # compliance version:
        #-----------------------
        elif self.model_version == 'compliance' :
            # get the corresponding macroscopic stresses
            sig_app_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )

            ### microplane equivalent stress obtained by projection (static constraint)
            s_vct_list = self._get_s_vct_list( sig_app_eng )
            s_equiv_list = self._get_s_equiv_list( s_vct_list )

            ### microplane equivalent strains calculated based on corresponding 'M' and 'psi_mtx'
            # 2nd order damage effect tensor:
            psi_mtx = self._get_psi_mtx( sctx, eps_app_eng )
            # 4th order damage effect tensor:
            if self.symmetrization == 'product-type':
                M4 = self._get_M_tns_product_type( psi_mtx )
            elif self.symmetrization == 'sum-type':
                M4 = self._get_M_tns_sum_type( psi_mtx )

            # apparent stress tensor:
            sig_app_mtx = map2d_sig_eng_to_mtx( sig_app_eng )
            # effective stress tensor:
            sig_eff_mtx = tensordot( M4, sig_app_mtx, [[2, 3], [0, 1]] )
            # effective strain tensor:
            eps_eff_mtx = tensordot( self.C4_e, sig_eff_mtx, [[2, 3], [0, 1]] )

            # effective microplane strains obtained by projection (kinematic constraint)
            e_eff_vct_list = array( [ dot( eps_eff_mtx, mpn ) for mpn in self._MPN ] )
            # apparent microplane strains 
            e_app_vct_list = array( [ dot( psi_mtx, e_eff_vct ) for e_eff_vct in e_eff_vct_list ] )
            # equivalent strain for each microplane
            e_equiv_list = self._get_e_equiv_list( e_app_vct_list )

            eps_app_mtx = map2d_eps_eng_to_mtx( eps_app_eng )
            print 'eps_app_mtx0', eps_app_mtx

#            # NOTE (1): compare with direct kinematic constraint
#            # this does NOT yield the same result!
#            e_app_vct_list2 = self._get_e_vct_list(eps_app_eng)
#            e_equiv_list2 = self._get_e_equiv_list(e_app_vct_list2)
#            print 'e_equiv_list:   ', e_equiv_list
#            print 'e_equiv_list2:  ', e_equiv_list2

#            # NOTE (2): the summation of the dyadic product of e_vct_list and n_vct yields the original strain tensor!
#            eps_app_mtx_list = array([outer(e_app_vct_list_i, MPN)*MPW for e_app_vct_list_i, MPN, MPW in zip(e_app_vct_list, self._MPN, self._MPW)])
#            eps_app_mtx = eps_app_mtx_list.sum(0)
#            if self.symmetrization == 'sum-type':
#                eps_app_mtx_sym = 0.5*(eps_app_mtx + transpose(eps_app_mtx))
#            elif self.symmetrization == 'product-type':
#                # @todo: verify if this in genneral is the correct way to symmetrize a matrix using product-type symmetrization
#                # Note that a problem arises here for negative values in eps_mtx.
#                # In phi_mtx only positiv values arise and therefore this way of 
#                # symmetrization is possible (cf. example [Jir99])  
#                eps_app_mtx_sym = sqrt(eps_app_mtx * transpose(eps_app_mtx))
#            print 'eps_app_mtx1', eps_app_mtx_sym

        return array( e_equiv_list ), array( s_equiv_list )


    def get_e_equiv_list( self, sctx, eps_app_eng ):
        '''
        Return a list of equivalent microplane strains consistently derived based on the 
        specified model version, i.e. stiffness or compliance.
        '''
        e_equiv_list, s_equiv_list = get_e_s_equiv_list( sctx, eps_app_eng )
        return e_equiv_list


    def get_s_equiv_list( self, sctx, eps_app_eng ):
        '''
        Return a list of equivalent microplane stresses consistently derived based on the 
        specified model version, i.e. stiffness or compliance.
        '''
        e_equiv_list, s_equiv_list = get_e_s_equiv_list( sctx, eps_app_eng )
        return s_equiv_list


    def get_fracture_energy( self, sctx, eps_app_eng ):
        '''
        Get the macroscopic fracture energy as a weighted sum of all mircoplane contributions
        '''
        e_max_list = self._get_state_variables( sctx, eps_app_eng )
        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
        return array( [dot( self._MPW, fracture_energy_list )], float )



# old implementation: this assumes a decoupled reaction of all microplanes
    def get_fracture_energy_list( self, sctx, eps_app_eng ):
        '''
        Get the microplane contributions to the fracture energy
        '''
        e_max_list = self._get_state_variables( sctx, eps_app_eng )
        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
        return fracture_energy_list

#    def get_fracture_energy(self, sctx, eps_app_eng):
#        '''
#        Get the macroscopic fracture energy as a weighted sum of all mircoplane contributions
#        '''
#        e_max_list = self._get_state_variables(sctx, eps_app_eng)
#        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
#        return array( [dot( self._MPW, fracture_energy_list )], float )

    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default( self ):
        return { 'eps_app'                  : self.get_eps_app,
                 'sig_app'                  : self.get_sig_app,
                 'sig_norm'                 : self.get_sig_norm,
                 'phi_pdc'                  : self.get_phi_pdc,
                 'microplane_damage'        : RTraceEval( eval = self.get_microplane_damage,
                                                          ts = self ),
                 'e_vct_list'               : self.get_e_vct_list,
                 'e_equiv'                  : self.get_e_equiv,
                 's_equiv'                  : self.get_s_equiv,
                 'fracture_energy_list'     : self.get_fracture_energy_list,
                 'fracture_energy'          : self.get_fracture_energy }






#--------------------------------------------------------------------------------
# Example 
#--------------------------------------------------------------------------------

# @todo - temporary alias rename the class and test it all

MA2DCompositeMicroplaneDamage = MATS2DMicroplaneDamage



#
# Supply the boundary conditions and construct the response graphs.
#
def get_value_and_coeff( a, alpha ):
    '''
    @TODO comment
    '''
    ca = cos( alpha )
    sa = sin( alpha )
    coeff = -( sa / ca )
    value = a / ca
    return value, coeff


from ibvpy.core.tloop import TLoop, TLine
from ibvpy.api import BCDof
from ibvpy.core.ibv_model import IBVModel

class MATS2DMicroplaneDamageExplore( IBVModel ):

    alpha_degree = Range( 0., 360., 0.,
                          label = 'Loading angle',
                          auto_set = False )

    ### 
    bcond_alpha = Instance( BCDof )
    def _bcond_alpha_default( self ):
        alpha = Pi * self.alpha_degree / 180
        value, coeff = get_value_and_coeff( 1., alpha )
        return  BCDof( var = 'u', dof = 0, value = value,
                      link_dofs = [1],
                      link_coeffs = [coeff],
                      time_function = lambda t: t )

    @on_trait_change( 'alpha_degree' )
    def update_bcond_alpha( self ):
        alpha = Pi * self.alpha_degree / 180
        value, coeff = get_value_and_coeff( 1., alpha )
        self.bcond_alpha.value = value
        self.bcond_alpha.link_coeffs[0] = coeff

    def _tloop_default( self ):

        from mats2D_cmdm_phi_fn import PhiFnStrainHardening
        elastic_debug = False
        # tseval for a material model
        #
#        tseval  = MATS2DMicroplaneDamage( elastic_debug = elastic_debug,
#                                        polar_fn_class = 'Anisotropic damage function' )
#        tseval.polar_fn.varied_params = ['Dfp']
        from mats2D_rtrace_cylinder import MATS2DRTraceCylinder

        tseval = MATS2DMicroplaneDamage( elastic_debug = elastic_debug )

        ts = TS( tse = tseval,
                 bcond_list = [ self.bcond_alpha
                             ],
                 rtrace_list = [ RTraceGraph( name = 'strain 0 - stress 0',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'strain 0 - stress 1',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'strain 0 - strain 1',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'eps_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'stress 0 - stress 1',
                                      var_x = 'sig_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'time - sig_norm',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'sig_norm', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'time - phi_pdc',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'phi_pdc', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'time - G_f',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'fracture_energy', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'e_equiv - s_equiv',
                                      var_x = 'e_equiv', idx_x = 0,
                                      var_y = 's_equiv', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph( name = 'time - microplane damage',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'microplane_damage', idx_y = 0,
                                      record_on = 'update' ),
                             MATS2DRTraceCylinder( name = 'Laterne',
                                      var_axis = 'time', idx_axis = 0,
                                      var_surface = 'microplane_damage',
                                      record_on = 'update' ),
                             RTraceArraySnapshot( name = 'fracture energy contributions',
                                      var = 'fracture_energy_list',
                                      record_on = 'update' ),
                             RTraceArraySnapshot( name = 'microplane damage',
                                      var = 'microplane_damage',
                                      record_on = 'update' ),
                             ]
                             )

        # Put the time-stepper into the time-loop
        #
        if elastic_debug:
            tmax = 1.
            n_steps = 1
        else:
            # Put the time-stepper into the time-loop
            #

            tmax = 0.001
            # tmax = 0.0006
            n_steps = 60
            #n_steps = 3

        tloop = TLoop( tstepper = ts,
                    KMAX = 100, RESETMAX = 0,
                    tline = TLine( min = 0.0, step = tmax / n_steps, max = tmax ) )

        return tloop

    traits_view = View( Item( 'alpha_degree@', label = 'angle' ),
                              resizable = True,
                              width = 1.0,
                              height = 1.0
                              )



def construct_fail_envelope():

    elastic_debug = False
    # Tseval for a material model
    #
    tseval = MATS2DMicroplaneDamage( elastic_debug = elastic_debug )


    value, coeff = get_value_and_coeff( 1., 0.0 )

    bcond_alpha = BCDof( var = 'u', dof = 0, value = value,
                     link_dofs = [1],
                     link_coeffs = [coeff],
                     time_function = lambda t: t )

    ts = TS( tse = tseval,
             bcond_list = [ self.bcond_alpha
                         ],
             rtrace_list = [ RTraceGraph( name = 'strain 0 - stress 0',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 0,
                                  record_on = 'update' ),
                         RTraceGraph( name = 'strain 1 - stress 1',
                                  var_x = 'eps_app', idx_x = 1,
                                  var_y = 'sig_app', idx_y = 1,
                                  record_on = 'update' ),
                         RTraceGraph( name = 'strain 0 - stress 1',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 1,
                                  record_on = 'update' ),
                         RTraceGraph( name = 'strain 1 - stress 0',
                                  var_x = 'eps_app', idx_x = 1,
                                  var_y = 'sig_app', idx_y = 0,
                                  record_on = 'update' ),
                         RTraceGraph( name = 'strain 0 - strain 1',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'eps_app', idx_y = 1,
                                  record_on = 'update' ),
                         ]
                         )

    # Put the time-stepper into the time-loop
    #
    if elastic_debug:
        tmax = 1.
        n_steps = 1
    else:
        tmax = 0.001
        # tmax = 0.0006
        n_steps = 100

    tl = TL( ts = ts,
             DT = tmax / n_steps, KMAX = 100, RESETMAX = 0,
             T = TRange( min = 0.0, max = tmax ) )

    from numpy import argmax

    alpha_arr = linspace( -Pi / 2 * 1.05, 2 * ( Pi / 2. ) + Pi / 2. * 0.05, 20 )

    sig0_m_list = []
    sig1_m_list = []
    eps0_m_list = []
    eps1_m_list = []

    for alpha in alpha_arr:

        value, coeff = get_value_and_coeff( 1., alpha )
        bcond_alpha.value = value
        bcond_alpha.link_coeffs[0] = coeff

        tl.eval()

        eps0_sig0 = tl.rv_mngr.rv_list[0]
        eps1_sig1 = tl.rv_mngr.rv_list[1]

        sig0_midx = argmax( fabs( eps0_sig0.trace.ydata ) )
        sig1_midx = argmax( fabs( eps1_sig1.trace.ydata ) )

        sig0_m = eps0_sig0.trace.ydata[ sig0_midx ]
        sig1_m = eps1_sig1.trace.ydata[ sig1_midx ]

        eps0_m = eps0_sig0.trace.xdata[ sig0_midx ]
        eps1_m = eps1_sig1.trace.xdata[ sig1_midx ]

        sig0_m_list.append( sig0_m )
        sig1_m_list.append( sig1_m )
        eps0_m_list.append( eps0_m )
        eps1_m_list.append( eps1_m )

    from math_func import MFnLineArray

    sig_plot = MFnLineArray( xdata = sig0_m_list,
                              ydata = sig1_m_list )
    eps_plot = MFnLineArray( xdata = eps0_m_list,
                              ydata = eps1_m_list )
    sig_plot.configure_traits()

    # Put the time-loop into the simulation-framework and map the
    # object to the user interface.
    #

if __name__ == '__main__':
    # consturct_fail_envelope()
    # Tseval for a material model
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp( ibv_resource = MATS2DMicroplaneDamageExplore() )
    ibvpy_app.main()
