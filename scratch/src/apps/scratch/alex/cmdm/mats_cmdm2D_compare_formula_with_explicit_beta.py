from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

# Chaco imports
from enthought.chaco2.chaco2_plot_editor import \
     Chaco2PlotEditor, \
     Chaco2PlotItem
from enthought.chaco2.chaco2_plot_container_editor import \
     PlotContainerEditor
from enthought.chaco2.tools.api import \
     PanTool, SimpleZoom
from enthought.chaco2.api import \
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

from ibvpy.api import \
    RTraceGraph,RTraceArraySnapshot

from mats2D_cmdm_polar_discr import \
    IPolarFn, IsotropicPolarFn, AnisotropicPolarFn


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
    def _get_polar_fn(self):
        return self.polar_fn_class_()
    
    def _set_polar_fn(self, value):
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
    
    model_version   = Enum("stiffness","compliance")
    stress_state    = Enum("plane_strain","plane_stress")
    symmetrization  = Enum("product-type","sum-type")

    dimensionality = Trait( '2D', {'2D' : 2, '3D' : 3 },
                 label = 'Dimensionality',
                 desc = 'This value is inactive yet',
                 auto_set = False)
        
    @on_trait_change('n_mp, dimensionality')
    def update_arrays( self ):
        self._init_arrays()
        self.changed = True
        
    elastic_debug = Bool( False,
                           desc = 'switch to elastic behavior - used for debugging',
                           auto_set = False)

    # This event can be used by the clients to trigger an action upon
    # the completed reconfiguration of the material model
    #
    changed = Event                     

    #---------------------------------------------------------------------------------------------
    # View specification
    #---------------------------------------------------------------------------------------------

    view_traits = View( VSplit( Group(Item('polar_fn_class', style='custom', show_label = False),
                                      Item('polar_fn', style='custom', show_label = False),
                                      label='Material parameters',
                                      show_border=True),
                                Group( Item('model_version', style = 'custom' ),
                                       Item('stress_state', style = 'custom' ),
                                       Item('symmetrization', style = 'custom' ),
                                       Item('dimensionality', style='custom'),
                                       Item('elastic_debug@'),
                                       Spring(resizable = True),
                                       label='Configuration parameters', show_border=True,
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
        Redefine the arrays
        '''
        # get the normal vectors of the microplanes
        self._MPN = array([[ cos( alpha ), sin( alpha )] for alpha in self.alpha_list ])
        # get the weights of the microplanes
        self._MPW = ones(self.n_mp) / self.n_mp * 2
        # get the dyadic product of the microplane normals
        self._MPNN = array( [ outer(mpn,mpn) for mpn in self._MPN ] )
         
    #-----------------------------------------------------------------------------------------------
    # Setup for computation within a supplied spatial context
    #-----------------------------------------------------------------------------------------------

    def get_state_array_size( self ):
        return self.n_mp

    def setup( self, sctx ):
        '''
        Intialize state variables.
        '''
        print '\n'
        print '--- SETTING UP ---'
        state_arr_size = self.get_state_array_size()
        sctx.mats_state_array = zeros(state_arr_size, 'float_')
        ##
        self.update_state_on = False
        self._check_dimensionality()
        self._setup_elasticity_tensors()
        print '\n'
        

    def _check_dimensionality( self ):
        # check dimensionality:
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
        E   = self.E
        nu  = self.nu
        # first Lame paramter
        la = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
        # second Lame parameter (shear modulus)
        mu = E / ( 2 + 2 * nu ) 

        # -----------------------------------------------------------------------------------------------------
        # Get the fourth order elasticity and compliance tensors for the 3D-case
        # -----------------------------------------------------------------------------------------------------
        D4_e_3D = zeros([3,3,3,3])
        C4_e_3D = zeros([3,3,3,3])
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
        D2_e_2D_plane_stress = self._get_D_plane_stress( D2_e_3D )
        D2_e_2D_plane_strain = self._get_D_plane_strain( D2_e_3D )
        C2_e_2D_plane_stress = self._get_C_plane_stress( C2_e_3D )
        C2_e_2D_plane_strain = self._get_C_plane_strain( C2_e_3D )

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
                
        
    def new_cntl_var(self):
        return zeros( 3, float_ )

    def new_resp_var(self):
        return zeros( 3, float_ )

    def _get_e_equiv_list(self, eps_app_eng ):
        '''
        Get microplane equivalent strains.
        '''        
                
        # tangent strain ratio
        c_T = self.c_T

        # apparent strain tensor
        eps_app_mtx = array( [[eps_app_eng[0]     , eps_app_eng[2] / 2.],
                              [eps_app_eng[2] / 2., eps_app_eng[1]]] )

        # apparent strain projection onto the individual microplanes
        e_app_vct_list = array( [ inner( eps_app_mtx, transpose( mpn ) ) for mpn in self._MPN ] )

        # magnitude of the normal strain vector for each microplane
        e_N_list = array( [ vdot( mpnn, eps_app_mtx ) for mpnn in self._MPNN ] )

        # positive part of the normal strain magnitude for each microplane
        e_N_pos_list = ( fabs( e_N_list ) + e_N_list ) / 2

        # normal strain vector for each microplane
        e_N_vct_list = array( [ self._MPN[i,:] * e_N_list[i] for i in range(0,self.n_mp) ] )

        # tangential strain vector for each microplane
        e_T_vct_list = e_app_vct_list - e_N_vct_list

        # squared tangential strain vector for each microplane
        e_TT_list = array( [ inner( e_T_vct, e_T_vct ) for e_T_vct in e_T_vct_list ] )

        # equivalent strain for each microplane
        e_equiv_list = sqrt( e_N_pos_list * e_N_pos_list + c_T * e_TT_list )
        
        return e_equiv_list
    

    def _get_e_max(self, e_equiv, e_max):
        '''
        compare the equivalent microplane strain with the 
        maximum strain reached in the loading history
        '''
        if e_equiv >= e_max:
            e_max = e_equiv
        return e_max
    
    
    def _get_state_variables(self, sctx, eps_app_eng):
        '''
        Compares the list of current equivalent microplane strains with 
        the values in the state array and returns the higher values 
        '''
        get_e_max_vectorized = frompyfunc( self._get_e_max, 2, 1 )
        e_equiv_list = self._get_e_equiv_list( eps_app_eng )
        e_max_list_old = sctx.mats_state_array
        e_max_list_new = get_e_max_vectorized(e_equiv_list, e_max_list_old)
        return e_max_list_new

    
    def _get_phi_list(self, sctx, eps_app_eng ):
        '''
        Get integrity factors for all microplanes.
        '''
        e_max_list = self._get_state_variables( sctx, eps_app_eng )    
        return self.polar_fn.get_phi_list( e_max_list )
    

    def _get_phi_mtx(self, sctx, eps_app_eng):
        '''
        Get the 2nd order damage tensor 'phi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_list   = self._get_phi_list( sctx, eps_app_eng )
        # integration terms for each microplanes
        phi_vct_list = array( [ phi_list[i] * self._MPNN[i,:,:] * self._MPW[i]
                                for i in range(0,self.n_mp) ] )
        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        phi_mtx = phi_vct_list.sum(0) 
        return phi_mtx


    def _get_psi_mtx(self, sctx, eps_app_eng):
        '''
        Get the 2nd order damage effect tensor 'psi_mtx'
        '''
        # scalar integrity factor for each microplane
        phi_list   = self._get_phi_list( sctx, eps_app_eng )
        
        # integration terms for each microplanes
        psi_vct_list = array( [ 1. / phi_list[i] * self._MPNN[i,:,:] * self._MPW[i]
                                for i in range(0,self.n_mp) ] )

        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        psi_mtx = psi_vct_list.sum(0) 
        return psi_mtx


    def _get_pdc_mtx(self, mtx):
        '''
        Transform the matrix phi_mtx (or psi_mtx) into principle damage coordinates and 
        assure that besides the diagonal the entries are exactly zero
        '''
        eig_value, eig_mtx = eig( mtx )
        eig_value_real = array([ pe.real for pe in eig_value] )
        pdc_mtx = zeros([self.n_dim,self.n_dim])
        for i in range(self.n_dim):
            pdc_mtx[i,i] = eig_value_real[i]
        return pdc_mtx

        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, d_eps, tn, tn1, eps_avg = None ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        
        #############################
        # old implementation using a stacked 3d phi_mtx or psi_mtx and a loop over an
        # formula for the components of the fourth order damage tensor or compliance tensor.
        # The 3D-matrixes are then reduced to 2D matrices. The results of this old implementation
        # are compared with the results of the new implementation using an explicitly define beta
        # tensor or damage effect tensor M.
        # Both implementations yield the same results !!! 
        
        n_d = 3
        
        n_mp = self.n_mp
        E   = self.E
        nu  = self.nu

        MPN  = self._MPN
        MPW  = self._MPW
        MPNN = self._MPNN
        # rank-four tensor including damage effect
        self.D4s_mdm = zeros([n_d,n_d,n_d,n_d])
        self.C4c_mdm = zeros([n_d,n_d,n_d,n_d])
        self.D4s_e   = zeros([n_d,n_d,n_d,n_d])
        self.C4c_e   = zeros([n_d,n_d,n_d,n_d])
        # rank-two tensor (matrix) including damage effect
        self.D2s_mdm = zeros([6,6])
        self.C2c_mdm = zeros([6,6])
        self.D2s_e   = zeros([6,6])
        self.C2c_e   = zeros([6,6])
        
        E   = self.E
        nu  = self.nu
        # first Lame paramter
        la = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
        # second Lame parameter (shear modulus)
        mu = E / ( 2 + 2 * nu ) 

        phi_list   = self._get_phi_list( sctx, eps_app_eng )
        
                
        # integration terms for each microplanes
        phi_vct_list = array( [ phi_list[i] * MPNN[i,:,:] * MPW[i]
                                for i in range(0,n_mp) ] )

        # integration terms for each microplanes
        psi_vct_list = array( [ 1. / phi_list[i] * MPNN[i,:,:] * MPW[i]
                                for i in range(0,n_mp) ] )

        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        phi_mtx2d = phi_vct_list.sum(0) 
        psi_mtx2d = psi_vct_list.sum(0) 

        phi_mtx = vstack( [hstack( [phi_mtx2d, [[0],[0]]] ), [[0,0,1]]] )
        psi_mtx = vstack( [hstack( [psi_mtx2d, [[0],[0]]] ), [[0,0,1]]] )
        
        
        # Lame constants calculated from E and nu
        # first Lame paramter calculated from E, nu
        la = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
        # second Lame parameter (shear modulus) calculated from E, nu
        mu = E / ( 2 + 2 * nu ) 
        
        n_d = 3
             
        # rank-four tensor including damage effect
        D4s_mdm = self.D4s_mdm
        C4c_mdm = self.C4c_mdm
        D4s_e = self.D4s_e
        C4c_e = self.C4c_e
        # rank-two tensor (matrix) including damage effect
        D2s_mdm = self.D2s_mdm
        C2c_mdm = self.C2c_mdm
        D2s_e = self.D2s_e
        C2c_e = self.C2c_e

        delta = identity(3)
        # Fill the rank-four matrices
        for i in range(0,n_d):
            for j in range(0,n_d):
                for k in range(0,n_d):
                    for l in range(0,n_d):
                        #--------------------------------------------------------------------------------
                        # Damaged material (formula by Bazant 1997)
                        #--------------------------------------------------------------------------------
                        D4s_mdm[i,j,k,l] = la * phi_mtx[i,j] * phi_mtx[k,l] + \
                                           mu * ( phi_mtx[i,k] * phi_mtx[j,l] + phi_mtx[i,l] * phi_mtx[j,k] )
                        # analogous derivation for the damaged compliance tensor (substitute delta by phi_mtx)
                        C4c_mdm[i,j,k,l] = (1+nu)/(2*E) * \
                                           ( psi_mtx[i,k] * psi_mtx[j,l] + psi_mtx[i,l]* psi_mtx[j,k] ) - \
                                           nu / E * psi_mtx[i,j] * psi_mtx[k,l]
                        #--------------------------------------------------------------------------------
                        # Elastic material 
                        #--------------------------------------------------------------------------------
                        # elastic stiffness tensor: 
                        # Note: in the 2D case (n_d = 2) this yields the matrix for plane strain)
                        D4s_e[i,j,k,l] = la * delta[i,j] * delta[k,l] + \
                                        mu * ( delta[i,k] * delta[j,l] + delta[i,l] * delta[j,k] )
                        # elastic compliance tensor:                 
                        C4c_e[i,j,k,l] = (1+nu)/(2*E) * \
                                         ( delta[i,k] * delta[j,l] + delta[i,l]* delta[j,k] ) - \
                                         nu / E * delta[i,j] * delta[k,l]                     
                        #--------------------------------------------------------------------------------
                        # map the components of the fourth-order tensor to the corresponding engineering notation
                        # using minor and major symmetry (=condensation to second-order tensor, i.e. 6x6-matrix)
                        #--------------------------------------------------------------------------------
                        D2s_mdm[map3d_ijkl2mn(i,j,k,l)] =  D4s_mdm[i,j,k,l]
                        C2c_mdm[map3d_ijkl2mn(i,j,k,l)] =  C4c_mdm[i,j,k,l]
                        #
                        D2s_e[map3d_ijkl2mn(i,j,k,l)] =  D4s_e[i,j,k,l]
                        C2c_e[map3d_ijkl2mn(i,j,k,l)] =  C4c_e[i,j,k,l]
                        
                        
        # compliance version:
        C2cc_mdm = self._compliance_mapping3d( C2c_mdm )
        D2c_mdm = inv( C2cc_mdm )
     

        # reduce 3d stiffness matrix to 2d
        if self.stress_state == 'plane_stress':
            D2r_mdm_s = self._get_D_plane_stress( D2s_mdm )
            D2r_mdm_c = self._get_D_plane_stress( D2c_mdm )
        else: # plane strain
            D2r_mdm_s = self._get_D_plane_strain( D2s_mdm )
            D2r_mdm_c = self._get_D_plane_strain( D2c_mdm )
        
        print 'D2r_mdm_s ', D2r_mdm_s, '\n'
        print 'D2r_mdm_c ', D2r_mdm_c, '\n'
        
        
        
        #############################
        
        
        
        
        
        
        
        # dimensionality:
        n_dim = self.n_dim

        
        # ------------------------------------------------------------------------------------------------
        # for debugging purposes only: if elastic_debug is switched on, linear elastic material is used 
        # ------------------------------------------------------------------------------------------------
        if self.elastic_debug:
            D2_e = self.D2_e
            sig_eng = tensordot( D2_e, eps_app_eng, [[1],[0]])
            return sig_eng, D2_e

        
        # ------------------------------------------------------------------------------------------------
        # update state variables
        # ------------------------------------------------------------------------------------------------
        if sctx.update_state_on:
            #print "in us"
            eps_n = eps_app_eng - d_eps
            e_max = self._get_state_variables( sctx, eps_n)
            sctx.mats_state_array[:] = e_max
        
        
        #-----------------------------------------------------------------------------------------------
        # stiffness version:
        #-----------------------------------------------------------------------------------------------
        
        if self.model_version == 'stiffness' :

            # get the damage tensor phi_mtx in form of a matrix
            phi_mtx = self._get_phi_mtx(sctx, eps_app_eng)
            
            # Get the direction of the principle damage coordinates (pdc):
            phi_eig_value, phi_eig_mtx = eig( phi_mtx )
            phi_eig_value_real = array([ pe.real for pe in phi_eig_value] )
            
                
            #---------------------------------------------------------------------------------------------
            # Product symmetrization using explicitly expressed beta tensor (=damage tensor)
            #---------------------------------------------------------------------------------------------
            # The case for product symmetrization yields exactly the same result as 
            # the formula for the damaged stiffness tensor (cf. [Baz97])
            if self.symmetrization == 'product-type':
                # transform phi_mtx to PDC:
                # (assure that besides the diagonal the entries are exactly zero)
                phi_pdc_mtx = zeros([n_dim,n_dim])
                for i in range(n_dim):
                    phi_pdc_mtx[i,i] = phi_eig_value_real[i]
                # w_mtx = tensorial square root of the second order damage tensor:
                w_pdc_mtx = arr_sqrt( phi_pdc_mtx )
                # transform the matrix w back to x-y-coordinates:
                w_mtx = dot(dot(phi_eig_mtx, w_pdc_mtx),transpose(phi_eig_mtx))
                # beta_ijkl = w_ik * w_jl (cf. [Baz 97])
                # exploiting numpy-functionality (faster).
                # Method 'outer' returns beta_ijkl = w_ij * w_kl,
                # therefore the axis 2 and 3 need to be swapped
                beta4_ = outer(w_mtx, w_mtx).reshape(n_dim,n_dim,n_dim,n_dim)
                beta4 = beta4_.swapaxes(1,2)
            
            
            #---------------------------------------------------------------------------------------------
            # Sum-type symmetrization using explicitly expressed beta tensor
            #---------------------------------------------------------------------------------------------
            elif self.symmetrization == 'sum-type':
                # (cf. [Jir99])
                beta4 = zeros([n_dim,n_dim,n_dim,n_dim])
                for i in range(0,n_dim):
                    for j in range(0,n_dim):
                        for k in range(0,n_dim):
                            for l in range(0,n_dim):
                                beta4[i,j,k,l] = 0.25 * ( phi_mtx[i,k] * delta[j,l] + phi_mtx[i,l] * delta[j,k] +\
                                                          phi_mtx[j,k] * delta[i,l] + phi_mtx[j,l] * delta[i,k] ) 
    
            #-----------------------------------------------------------------------------------
            # Calculate the damaged stiffness tensor based on the damage tensor beta4:
            #-----------------------------------------------------------------------------------
            D4_mdm = tensordot( beta4, tensordot( self.D4_e, beta4, [[2,3],[2,3]] ), [[2,3],[0,1]] )

            #-----------------------------------------------------------------------------------
            # Reduce the fourth order tensor to a matrix assuming minor and major symmetry
            #-----------------------------------------------------------------------------------

            # 2D-case:
            if n_dim == 2:
                D2_mdm = map2d_tns4_to_tns2(D4_mdm)
            # 3D-case:
            else:
                D4_mdm = map3d_tns4_to_tns2(D4_mdm)
     
                    
        #-----------------------------------------------------------------------------------------------
        # compliance version:
        #-----------------------------------------------------------------------------------------------
        
        elif self.model_version == 'compliance' :

            # get the damage tensor phi_mtx in form of a matrix
            psi_mtx = self._get_psi_mtx(sctx, eps_app_eng)

            # Get the direction of the principle damage coordinates (pdc):
            psi_eig_value, psi_eig_mtx = eig( psi_mtx )
            psi_eig_value_real = array([ pe.real for pe in psi_eig_value] )

                    
            #---------------------------------------------------------------------------------------------
            # Product symmetrization using explicitly expressed M tensor (=damage effect tensor)
            #---------------------------------------------------------------------------------------------
            if self.symmetrization == 'product-type':
                # transform phi_mtx to PDC:
                # (assure that besides the diagonal the entries are exactly zero)
                psi_pdc_mtx = zeros([n_dim,n_dim])
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
                # therefore the axis 2 and 3 need to be swapped
                M4_ = outer(w_hat_mtx, w_hat_mtx).reshape(n_dim,n_dim,n_dim,n_dim)
                M4 = M4_.swapaxes(1,2)
            
            #---------------------------------------------------------------------------------------------
            # Sum-type symmetrization using explicitly expressed beta tensor
            #---------------------------------------------------------------------------------------------
            elif self.symmetrization == 'sum-type':
                # (cf. [Jir99])
                M4 = zeros([n_dim,n_dim,n_dim,n_dim])
                for i in range(0,n_dim):
                    for j in range(0,n_dim):
                        for k in range(0,n_dim):
                            for l in range(0,n_dim):
                                M4[i,j,k,l] = 0.25 * ( psi_mtx[i,k] * delta[j,l] + psi_mtx[i,l] * delta[j,k] +\
                                                       psi_mtx[j,k] * delta[i,l] + psi_mtx[j,l] * delta[i,k] ) 
    
            #-----------------------------------------------------------------------------------
            # Calculate the damaged compliance tensor based on the damage effect tensor M4:
            #-----------------------------------------------------------------------------------
            C4_mdm = tensordot( M4, tensordot( self.C4_e, M4, [[2,3],[2,3]] ), [[2,3],[0,1]] )
            
            #-----------------------------------------------------------------------------------
            # Reduce the fourth order tensor to a matrix assuming minor and major symmetry and 
            # multiply the resulting D matrix with the factors for compliance mapping
            #-----------------------------------------------------------------------------------

            # 2D-case:
            if n_dim == 2:
                C2_mdm = map2d_tns4_to_tns2(C4_mdm)
                D2_mdm = inv(self._compliance_mapping2d( C2_mdm ))
            # 3D-case:
            else:
                C2_mdm = map3d_tns4_to_tns2(C4_mdm)
                D2_mdm = inv(self._compliance_mapping3d( C2_mdm ))
            
        #-----------------------------------------------------------------------------------
        # Return stresses (corrector) and damaged secant stiffness matrix (predictor)
        #-----------------------------------------------------------------------------------
        sig_eng = tensordot( D2_mdm, eps_app_eng, [[1],[0]])
        
        print 'D2_mdm ', D2_mdm
        
        return sig_eng, D2_mdm



    #---------------------------------------------------------------------------------------------
    # Subsidiary methods realizing configurable features
    #---------------------------------------------------------------------------------------------
  
    def _compliance_mapping2d( self, C_mtx_2d ):
        '''
        Compliance matrix must reflect the symmetry assumption 
        The components of the compliance matrix are multiplied with factor 1,2 or 4 depending on their 
        position in the matrix (due to symmetry and switching of tensorial to engineering notation
        (Note: gamma_xy = 2*epsilon_xy, etc.). Necessary when evaluating D=inv(C). 
        '''
        idx1 = [0,1]
        idx2 = [2]
        C11 = C_mtx_2d[ix_(idx1,idx1)]
        C12 = C_mtx_2d[ix_(idx1,idx2)]
        C21 = C_mtx_2d[ix_(idx2,idx1)]
        C22 = C_mtx_2d[ix_(idx2,idx2)]
        return vstack( [ hstack( [C11,   2*C12] ),
                         hstack( [2*C21, 4*C22] ) ] )


    def _compliance_mapping3d( self, C_mtx_3d ):
        '''
        Compliance matrix must reflect the symmetry assumption 
        The components of the compliance matrix are multiplied with factor 1,2 or 4 depending on their 
        position in the matrix (due to symmetry and switching of tensorial to engineering notation
        (Note: gamma_xy = 2*epsilon_xy, etc.). Necessary when evaluating D=inv(C). 
        '''
        idx1 = [0,1,2]
        idx2 = [3,4,5]
        C11 = C_mtx_3d[ix_(idx1,idx1)]
        C12 = C_mtx_3d[ix_(idx1,idx2)]
        C21 = C_mtx_3d[ix_(idx2,idx1)]
        C22 = C_mtx_3d[ix_(idx2,idx2)]
        return vstack( [ hstack( [C11,   2*C12] ),
                         hstack( [2*C21, 4*C22] ) ] )


    def _get_D_plane_stress( self, D_mtx_3d ):
        '''reduce the 6x6-elasticity matrix for the 3D-case to a 
        3x3 matrix for the 2D case assuming plane stress (sig_yz=0, sig_zz=0, sig_xz=0)
        '''
        idx2 = [0,1,5]
        idx3 = [2,3,4]
        D22 = D_mtx_3d[ix_(idx2,idx2)]
        D23 = D_mtx_3d[ix_(idx2,idx3)]
        D32 = D_mtx_3d[ix_(idx3,idx2)]
        D33 = D_mtx_3d[ix_(idx3,idx3)]
        D_pstress_term = dot( dot( D23, inv( D33 ) ), D32 )
        D_mtx_2d = D22 - D_pstress_term
        return D_mtx_2d


    def _get_D_plane_strain( self, D_mtx_3d ):
        '''reduce the 6x6-elasticity matrix for the 3D-case to a 
        3x3 matrix for the 2D case assuming plane strain (eps_yz=0, eps_zz=0, eps_xz=0)
        '''
        idx2 = [0,1,5]
        D_mtx_2d = D_mtx_3d[ix_(idx2,idx2)]
        return D_mtx_2d


    def _get_C_plane_stress( self, C_mtx_3d ):
        '''reduce the 6x6-compliance matrix for the 3D-case to a 
        3x3 matrix for the 2D case assuming plane stress (sig_yz=0, sig_zz=0, sig_xz=0)
        '''
        idx2 = [0,1,5]
        C_mtx_2d = C_mtx_3d[ix_(idx2,idx2)]
        return C_mtx_2d


    def _get_C_plane_strain( self, C_mtx_3d ):
        '''reduce the 6x6-compliance matrix for the 3D-case to a 
        3x3 matrix for the 2D case assuming plane strain (eps_yz=0, eps_zz=0, eps_xz=0)
        '''
        idx2 = [0,1,5]
        idx3 = [2,3,4]
        C22 = C_mtx_3d[ix_(idx2,idx2)]
        C23 = C_mtx_3d[ix_(idx2,idx3)]
        C32 = C_mtx_3d[ix_(idx3,idx2)]
        C33 = C_mtx_3d[ix_(idx3,idx3)]
        C_pstrain_term = dot( dot( C23, inv( C33 ) ), C32 )
        C_mtx_2d = C22 - C_pstrain_term
        return C_mtx_2d
    

    #---------------------------------------------------------------------------------------------
    # Update state method called upon an accepted time-step
    #---------------------------------------------------------------------------------------------

    def update_state(self, sctx, eps_app_eng ):
        '''
        Here just set the flag on to make the update afterwards in the method itself.
        '''
        self.update_state_on = True

    #---------------------------------------------------------------------------------------------
    # Response trace evaluators
    #---------------------------------------------------------------------------------------------

    def get_eps_app( self, sctx, eps_app_eng ):
        return eps_app_eng
    
    def get_eps_integ( self, sctx, eps_app_eng ):
        eps_app_eng
        
        return 

    def get_sig_app( self, sctx, eps_app_eng ):
        # @TODO
        # the stress calculation is performed twice - it might be
        # cached. But not in the spatial integration scheme.
        # 
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return sig_eng

    def get_sig_norm( self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0, 0 )
        return array( [ scalar_sqrt( sig_eng[0]**2 + sig_eng[1]**2 ) ] )
    
    def get_phi_pdc( self, sctx, eps_app_eng ):
        phi_mtx = self._get_phi_mtx( sctx, eps_app_eng )
        # Get the direction of the principle damage coordinates (pdc):
        phi_eig_value, phi_eig_mtx = eig( phi_mtx )
        phi_eig_value_real = array([ pe.real for pe in phi_eig_value] )
        phi_pdc_mtx = zeros([self.n_dim,self.n_dim])
        for i in range(self.n_dim):
            phi_pdc_mtx[i,i] = phi_eig_value_real[i]
        return array( [ phi_pdc_mtx[0,0] ])
    
    def get_microplane_damage(self, sctx, eps_app_eng ):
        phi_list = self._get_phi_list(sctx, eps_app_eng)
        return phi_list

    def get_fracture_energy_list(self, sctx, eps_app_eng):
        '''
        Get the microplane contributions to the fracture energy
        '''
        e_max_list = self._get_state_variables(sctx, eps_app_eng)
        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
        return fracture_energy_list

    def get_fracture_energy(self, sctx, eps_app_eng):
        '''
        Get the macroscopic fracture energy as a weighted sum of all mircoplane contributions
        '''
        e_max_list = self._get_state_variables(sctx, eps_app_eng)
        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
        return array( [dot( self._MPW, fracture_energy_list )], float )
    
    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return { 'eps_app'                  : self.get_eps_app,
                 'sig_app'                  : self.get_sig_app,
                 'sig_norm'                 : self.get_sig_norm,
                 'phi_pdc'                  : self.get_phi_pdc,
                 'microplane_damage'        : self.get_microplane_damage,
                 'fracture_energy_list'     : self.get_fracture_energy_list,
                 'fracture_energy'          : self.get_fracture_energy }

#---------------------------------------------------------------------------------------------
# Subsidiary index mapping functions for rank-four to rank-two tensors
#---------------------------------------------------------------------------------------------

def map2d_ijkl2mn(i,j,k,l):
    '''
    Map the four-rank indexes to the two-rank matrix using the major
    and minor symmetry.
    '''
    # 2D-case:
    # first two indices (ij)
    if i==0 and j==0:
        m = 0
    elif i==1 and j==1:
        m = 1
    elif (i==0 and j==1) or (i==1 and j==0):
        m = 2
        
    # second two indices (kl)
    if k==0 and l==0:
        n = 0
    elif k==1 and l==1:
        n = 1
    elif (k==0 and l==1) or (k==1 and l==0):
        n = 2
        
    return m,n

def map3d_ijkl2mn(i,j,k,l):
    '''
    Map the four-rank indexes to the two-rank matrix using the major
    and minor symmetry.
    '''
    # 3D-case:
    # first two indices (ij)
    if i==0 and j==0:
        m = 0
    elif i==1 and j==1:
        m = 1
    elif i==2 and j==2:
        m = 2
    elif (i==1 and j==2) or (i==2 and j==1):
        m = 3
    elif (i==0 and j==2) or (i==2 and j==0):
        m = 4
    elif (i==0 and j==1) or (i==1 and j==0):
        m = 5
    else:
        raise IndexError, 'error in the tensor index mapping'
        
    # second two indices (kl)
    if k==0 and l==0:
        n = 0
    elif k==1 and l==1:
        n = 1
    elif k==2 and l==2:
        n = 2
    elif (k==1 and l==2) or (k==2 and l==1):
        n = 3
    elif (k==0 and l==2) or (k==2 and l==0):
        n = 4
    elif (k==0 and l==1) or (k==1 and l==0):
        n = 5
    else:
        raise IndexError, 'error in the tensor index mapping'

    return m,n

#---------------------------------------------------------------------------------------------
# Subsidiary mapping functions for rank-two to rank-four tensor 
#---------------------------------------------------------------------------------------------

def map2d_tns2_to_tns4(tns2):
    '''map a matrix to a fourth order tensor assuming minor and major symmetry,
    e.g. D_mtx (3x3) in engineering notation to D_tns(2,2,2,2)).
    '''
    n_dim = 2
    tns4 = zeros([n_dim,n_dim,n_dim,n_dim])
    for i in range(0,n_dim):
        for j in range(0,n_dim):
            for k in range(0,n_dim):
                for l in range(0,n_dim):
                    tns4[i,j,k,l] = tns2[map2d_ijkl2mn(i,j,k,l)]
    return tns4
               
def map3d_tns2_to_tns4(tns2):
    '''map a matrix to a fourth order tensor assuming minor and major symmetry,
    e.g. D_mtx (6x6) in engineering notation to D_tns(3,3,3,3)).
    '''
    n_dim = 3
    tns4 = zeros([n_dim,n_dim,n_dim,n_dim])
    for i in range(0,n_dim):
        for j in range(0,n_dim):
            for k in range(0,n_dim):
                for l in range(0,n_dim):
                    tns4[i,j,k,l] = tns2[map3d_ijkl2mn(i,j,k,l)]
    return tns4
               

#---------------------------------------------------------------------------------------------
# Subsidiary mapping functions for rank-four to rank-two tensor 
#---------------------------------------------------------------------------------------------

def map3d_tns4_to_tns2(tns4):
    '''map a fourth order tensor to a matrix assuming minor and major symmetry,
    e.g. D_tns(3,3,3,3) to D_mtx (6x6) in engineering notation.
    (Note: Explicit assignment of components used for speedup.)
    '''
    n_eng = 6
    tns2 = zeros([n_eng,n_eng])
    
    tns2[0,0] =             tns4[0,0,0,0]
    tns2[0,1] = tns2[1,0] = tns4[0,0,1,1]
    tns2[0,2] = tns2[2,0] = tns4[0,0,2,2]
    tns2[0,3] = tns2[3,0] = tns4[0,0,1,2]
    tns2[0,4] = tns2[4,0] = tns4[0,0,0,2]
    tns2[0,5] = tns2[5,0] = tns4[0,0,0,1]
    
    tns2[1,1] =             tns4[1,1,1,1]
    tns2[1,2] = tns2[2,1] = tns4[1,1,2,2]
    tns2[1,3] = tns2[3,1] = tns4[1,1,1,2]
    tns2[1,4] = tns2[4,1] = tns4[1,1,0,2]
    tns2[1,5] = tns2[5,1] = tns4[1,1,0,1]
    
    tns2[2,2] =             tns4[2,2,2,2]
    tns2[2,3] = tns2[3,2] = tns4[2,2,1,2]
    tns2[2,4] = tns2[4,2] = tns4[2,2,0,2]
    tns2[2,5] = tns2[5,2] = tns4[2,2,0,1]
    
    tns2[3,3] =             tns4[1,2,1,2]
    tns2[3,4] = tns2[4,3] = tns4[1,2,0,2]
    tns2[3,5] = tns2[5,3] = tns4[1,2,0,1]
    
    tns2[4,4] =             tns4[0,2,0,2]
    tns2[4,5] = tns2[5,4] = tns4[0,2,0,1]
    
    tns2[5,5] =             tns4[0,1,0,1]
    
    return tns2

def map2d_tns4_to_tns2(tns4):
    '''map a fourth order tensor to a matrix assuming minor and major symmetry,
    e.g. D_tns(2,2,2,2) to D_mtx (3x3) in engineering notation.
    (Note: Explicit assignment of components used for speedup.)
    '''
    n_eng = 3
    tns2 = zeros([n_eng,n_eng])
    
    tns2[0,0] =             tns4[0,0,0,0]
    tns2[0,1] = tns2[1,0] = tns4[0,0,1,1]
    tns2[0,2] = tns2[2,0] = tns4[0,0,0,1]
    
    tns2[1,1] =             tns4[1,1,1,1]
    tns2[1,2] = tns2[2,1] = tns4[1,1,0,1]
    
    tns2[2,2] =             tns4[0,1,0,1]
    
    return tns2





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
    ca = cos(alpha)
    sa = sin(alpha)
    coeff = -(sa/ca)
    value = a / ca
    return value, coeff


from ibvpy.core.tloop import TLoop, TLine
from ibvpy.api import BCDof
from ibvpy.core.ibvp_solve import IBVPSolve as IS

class MATS2DMicroplaneDamageExplore( HasTraits ):

    alpha_degree = Range( 0., 360., 0., 
                          label = 'Loading angle',
                          auto_set = False)

    ### 
    bcond_alpha = Instance( BCDof )
    def _bcond_alpha_default( self ):
        alpha =  Pi * self.alpha_degree / 180
        value, coeff = get_value_and_coeff( 1., alpha )
        return  BCDof(var='u', dof = 0, value = value,
                      link_dofs = [1],
                      link_coeffs = [coeff],
                      time_function = lambda t: t )

    @on_trait_change('alpha_degree')
    def update_bcond_alpha( self ):
        alpha =  Pi * self.alpha_degree / 180
        value, coeff = get_value_and_coeff( 1., alpha )
        self.bcond_alpha.value = value
        self.bcond_alpha.link_coeffs[0] = coeff

    sim = Instance( IS )
    def _sim_default( self ):
        
        from mats2D_cmdm_phi_fn import PhiFnStrainHardening
        elastic_debug = False
        # tseval for a material model
        #
#        tseval  = MATS2DMicroplaneDamage( elastic_debug = elastic_debug,
#                                        polar_fn_class = 'Anisotropic damage function' )
#        tseval.polar_fn.varied_params = ['Dfp']

        tseval  = MATS2DMicroplaneDamage( elastic_debug = elastic_debug )
        
        ts = TS( tse = tseval,
                 bcond_list = [ self.bcond_alpha 
                             ],
                 rtrace_list = [ RTraceGraph(name = 'strain 0 - stress 0',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph(name = 'strain 1 - stress 1',
                                      var_x = 'eps_app', idx_x = 1,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTraceGraph(name = 'strain 0 - stress 1',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTraceGraph(name = 'strain 1 - stress 0',
                                      var_x = 'eps_app', idx_x = 1,
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
                             RTraceGraph(name = 'time - sig_norm',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'sig_norm', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph(name = 'time - phi_pdc',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'phi_pdc', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceGraph(name = 'time - G_f',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'fracture_energy', idx_y = 0,
                                      record_on = 'update' ),
                             RTraceArraySnapshot(name = 'fracture energy contributions',
                                      var = 'fracture_energy_list',
                                      record_on = 'update' ),
                             RTraceArraySnapshot(name = 'microplane damage',
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
    
            tl = TLoop( tstepper = ts,
                     DT=tmax/n_steps, KMAX = 100, RESETMAX = 0,
                     tline = TLine( min = 0.0,  max = tmax ) )
    
            return IS( tloop = tl )

    traits_view = View( Item('alpha_degree@', label = 'angle'),
                        Item('bcond_alpha@'),
                        Item('sim@'),
                        resizable = True,
                        width = 1.0,
                        height = 1.0
                        )

    traits_view = View( Item('alpha_degree@', label = 'angle'),
                        Item('bcond_alpha@'),
                        Item('sim@'),
                        resizable = True,
                        width = 1.0,
                        height = 1.0
                        )

    traits_view_alpha = View( Item('alpha_degree@', label = 'angle'),
                              resizable = True,
                              width = 1.0,
                              height = 1.0
                              )


    
def construct_fail_envelope():

    elastic_debug = False
    # Tseval for a material model
    #
    tseval  = MATS2DMicroplaneDamage( elastic_debug = elastic_debug )


    value, coeff = get_value_and_coeff( 1., 0.0 )

    bcond_alpha = BCDof(var='u', dof = 0, value = value,
                     link_dofs = [1],
                     link_coeffs = [coeff],
                     time_function = lambda t: t )
    
    ts = TS( tse = tseval,
             bcond_list = [ self.bcond_alpha 
                         ],
             rtrace_list = [ RTraceGraph(name = 'strain 0 - stress 0',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 0,
                                  record_on = 'update' ),
                         RTraceGraph(name = 'strain 1 - stress 1',
                                  var_x = 'eps_app', idx_x = 1,
                                  var_y = 'sig_app', idx_y = 1,
                                  record_on = 'update' ),
                         RTraceGraph(name = 'strain 0 - stress 1',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 1,
                                  record_on = 'update' ),
                         RTraceGraph(name = 'strain 1 - stress 0',
                                  var_x = 'eps_app', idx_x = 1,
                                  var_y = 'sig_app', idx_y = 0,
                                  record_on = 'update' ),
                         RTraceGraph(name = 'strain 0 - strain 1',
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
             DT=tmax/n_steps, KMAX = 100, RESETMAX = 0,
             T = TRange( min = 0.0,  max = tmax ) )

    from numpy import argmax

    alpha_arr = linspace( - Pi/2 * 1.05,  2*(Pi/2.) + Pi/2.*0.05, 20 )

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
    mme = MATS2DMicroplaneDamageExplore()
    #mme.configure_traits( view = 'traits_view_begehung' )
    mme.configure_traits()
