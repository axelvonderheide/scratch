
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

#from dacwt import DAC

from numpy import \
     array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt

from math import pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from scipy.linalg import eig, inv

from core.time_stepper import \
     TimeStepper as TS

from core.mats_eval import IMATSEval, MATSEval

from core.rv import RTO, RTOGraph, RTOArraySnapshot
from mats2D_cmdm_polar_discr import \
    IPolarFn, IsotropicPolarFn, AnisotropicPolarFn

from time import *

#---------------------------------------------------------------------------
# Material time-step-evaluator for Microplane-Damage-Model
#---------------------------------------------------------------------------

class MATS2DMicroplaneDamage( MATSEval ):
    '''
    Microplane Damage Model.
    '''
    implements( IMATSEval )

    #---------------------------------------------------------------------------
    # Sphere discretization in form of an array of microplanes
    #---------------------------------------------------------------------------
    
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


    n_d = Trait( '2D', {'2D' : 2, '3D' : 3 },
                 label = 'Dimensionality',
                 desc = 'This value is inactive yet',
                 auto_set = False)

#    if n_d == '2D':
#        n_eng = Int(3)
#        print 'n_eng = 3'
#    else:
#        n_eng = Int(6)
#        print 'n_eng = 3'
        
        
    @on_trait_change('n_mp')
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
                                       Item('n_d', style='custom'),
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
        print 'SETTING UP'
        state_arr_size = self.get_state_array_size()
        sctx.mats_state_array = zeros(state_arr_size, 'float_')
        ##
        self.update_state_on = False
        self._setup_interim_arrays()
        self._setup_elasticity_tensors()
        
        
    def _setup_interim_arrays( self ):
        '''
        Intialize intermediate arrays
        '''
        t1 = time() 
#        n_d = self.n_d
#        n_eng = self.n_eng
        #
        # 2D case:
        n_d = 2
        n_eng = 3
        # rank-four tensor including damage effect
        self.D4s_mdm_2D = zeros([n_d,n_d,n_d,n_d])
        self.C4c_mdm_2D = zeros([n_d,n_d,n_d,n_d])
        self.D4s_e_2D   = zeros([n_d,n_d,n_d,n_d])
        self.C4c_e_2D   = zeros([n_d,n_d,n_d,n_d])
        # rank-two tensor (matrix) including damage effect
        self.D2s_mdm_2D = zeros([n_eng,n_eng])
        self.C2c_mdm_2D = zeros([n_eng,n_eng])
        self.D2s_e_2D   = zeros([n_eng,n_eng])
        self.C2c_e_2D   = zeros([n_eng,n_eng])
        # 3D case:
        n_d = 3
        n_eng = 6
        # rank-four tensor including damage effect
        self.D4s_mdm_3D = zeros([n_d,n_d,n_d,n_d])
        self.C4c_mdm_3D = zeros([n_d,n_d,n_d,n_d])
        self.D4s_e_3D   = zeros([n_d,n_d,n_d,n_d])
        self.C4c_e_3D   = zeros([n_d,n_d,n_d,n_d])
        # rank-two tensor (matrix) including damage effect
        self.D2s_mdm_3D = zeros([n_eng,n_eng])
        self.C2c_mdm_3D = zeros([n_eng,n_eng])
        self.D2s_e_3D   = zeros([n_eng,n_eng])
        self.C2c_e_3D   = zeros([n_eng,n_eng])
        t2 = time()
        self.t_interim = t2 - t1
        
    def _setup_elasticity_tensors( self ):
        '''
        Intialize the fourth order elasticity tensor for plane strain or plane stress
        '''
        t1 = time()        
        
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
        # Get the fourth order eleasticity and compliance tensors for the 3d-case
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
        # Get the fourth order elasticity and compliance tensors for the 2D-cases
        # -----------------------------------------------------------------------------------------------------
        # Get the (6x6)-elasticity and compliance matrices (for the 3D-case)
        D2_e_3D = map3d_tns4_to_tns2(D4_e_3D)
        C2_e_3D = map3d_tns4_to_tns2(C4_e_3D)

        # Get the (3x3)-elasticity and compliance matrices for the 2D-cases plane stress and plane strain
        D2_e_2D_plane_stress = self._get_D_plane_stress( D2_e_3D )
        D2_e_2D_plane_strain = self._get_D_plane_strain( D2_e_3D )
        C2_e_2D_plane_stress = self._get_C_plane_stress( C2_e_3D )
        C2_e_2D_plane_strain = self._get_C_plane_strain( C2_e_3D )

        # Get the fourth order elasticity and compliance tensors for the 2D-cases plane stress and plane strain
        D4_e_2D_plane_stress = map2d_tns2_to_tns4( D2_e_2D_plane_stress )
        D4_e_2D_plane_strain = map2d_tns2_to_tns4( D2_e_2D_plane_strain )
        C4_e_2D_plane_stress = map2d_tns2_to_tns4( C2_e_2D_plane_stress )
        C4_e_2D_plane_strain = map2d_tns2_to_tns4( C2_e_2D_plane_strain )

        # -----------------------------------------------------------------------------------------------------
        # assign the fourth order elasticity and compliance tensors to self
        # -----------------------------------------------------------------------------------------------------
        # (3D case):                                       
#        if self.n_d == '3D':
        self.D4_e_3D = D4_e_3D
        self.C4_e = C4_e_3D

        # (2D case):
#        if self.n_d == '2D' and self.stress_state == 'plane_stress':
        self.D4_e_2D_plane_stress = D4_e_2D_plane_stress 
        self.C4_e = C4_e_2D_plane_stress 

#        if self.n_d == '2D' and self.stress_state == 'plane_strain':
        self.D4_e_2D_plane_strain = D4_e_2D_plane_strain
        self.C4_e = C4_e_2D_plane_strain

        t2 = time()        
        self.t_setup_De = t2 - t1
        
        self.number_of_steps = 0
        
        self.t_formula2D = 0
        self.t_formula3D = 0
        self.t_bb2D = 0
        self.t_bb3D = 0
        self.t_reduction = 0

        #-----------------------------------------------------------------------------------
        # Get the fourth order elastic stiffness tensor      
        #-----------------------------------------------------------------------------------
        # for the plane stress state (2D-case)
        D2_e_2Dstress = zeros([3,3])
        D2_e_2Dstress[0, 0] = 4.0 * (la + mu) * mu / (la + 2.0 * mu)
        D2_e_2Dstress[0, 1] = 2.0 * mu * la / (la + 2.0 * mu)
        D2_e_2Dstress[1, 0] = 2.0 * mu * la / (la + 2.0 * mu)
        D2_e_2Dstress[1, 1] = 4.0 * (la + mu) * mu / (la + 2.0 * mu)
        D2_e_2Dstress[2, 2] = mu
        # map the components of the (3,3)-stiffness matrix into the corresponding tensor notation
        D4_e_2Dstress = zeros([2,2,2,2])
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    for l in range(0,2):
                        # map the components of the stiffness matrix (3,3) into the corresponding tensor notation
                        D4_e_2Dstress[i,j,k,l] = D2_e_2Dstress[map2d_ijkl2mn(i,j,k,l)]           
        self.D4_e_2Dstress = D4_e_2Dstress
        
        # for the plane strain state (2D-case)
        D2_e_2Dstrain = zeros([3,3])
        D2_e_2Dstrain[0, 0] = la + 2 * mu
        D2_e_2Dstrain[0, 1] = la
        D2_e_2Dstrain[1, 0] = la
        D2_e_2Dstrain[1, 1] = la + 2 * mu
        D2_e_2Dstrain[2, 2] = mu
        # map the components of the (3,3)-stiffness matrix into the corresponding tensor notation
        D4_e_2Dstrain = zeros([2,2,2,2])
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    for l in range(0,2):
                        # map the components of the stiffness matrix (3,3) into the corresponding tensor notation
                        D4_e_2Dstrain[i,j,k,l] = D2_e_2Dstrain[map2d_ijkl2mn(i,j,k,l)]           
        self.D4_e_2Dstrain = D4_e_2Dstrain
#### 
        

        
    def new_cntl_var(self):
        return zeros( 3, float_ )

    def new_resp_var(self):
        return zeros( 3, float_ )

    def get_e_equiv_list(self, eps_app_eng ):
        '''
        Get microplane equivalent strains.
        '''
        self.number_of_steps += 1

        n_mp = self.n_mp
        
        # Prepare the microplane coefficients needed for the 
        # projection and integration
        MPN  = self._MPN
        MPNN = self._MPNN
        
        # tangent strain ratio
        c_T = self.c_T

        # apparent strain tensor
        eps_app_mtx = array( [[eps_app_eng[0]     , eps_app_eng[2] / 2.],
                              [eps_app_eng[2] / 2., eps_app_eng[1]]] )

        # apparent strain projection onto the individual microplanes
        e_app_vct_list = array( [ inner( eps_app_mtx, transpose( mpn ) ) for mpn in MPN ] )

        # magnitude of the normal strain vector for each microplane
        e_N_list = array( [ vdot( mpnn, eps_app_mtx ) for mpnn in MPNN ] )

        # positive part of the normal strain magnitude for each microplane
        e_N_pos_list = ( fabs( e_N_list ) + e_N_list ) / 2

        # normal strain vector for each microplane
        e_N_vct_list = array( [ MPN[i,:] * e_N_list[i] for i in range(0,n_mp) ] )

        # tangential strain vector for each microplane
        e_T_vct_list = e_app_vct_list - e_N_vct_list

        # squared tangential strain vector for each microplane
        e_TT_list = array( [ inner( e_T_vct, e_T_vct ) for e_T_vct in e_T_vct_list ] )

        # equivalent strain for each microplane
        e_equiv_list = sqrt( e_N_pos_list * e_N_pos_list + c_T * e_TT_list )
        return e_equiv_list
    

    def _get_e_max(self, e_equiv, e_max):
        if e_equiv >= e_max:
            e_max = e_equiv
        return e_max
    
    def get_e_max_list(self, sctx, eps_app_eng):
        '''
        Compare the current equivalent strain with the values in the state array
        and return the higher values 
        '''
        get_e_max_vect = frompyfunc( self._get_e_max, 2, 1 )
        e_equiv_list = self.get_e_equiv_list( eps_app_eng )
        e_max_list_old = sctx.mats_state_array
        e_max_list = get_e_max_vect(e_equiv_list, e_max_list_old)
        return e_max_list

    
    def _get_phi_list(self, sctx, eps_app_eng ):
        '''
        Get integrity factors for all microplanes.
        '''
        e_max_list = self.get_e_max_list( sctx, eps_app_eng )    
        return self.polar_fn.get_phi_list( e_max_list )
        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        n_mp = self.n_mp

        MPN  = self._MPN
        MPW  = self._MPW
        MPNN = self._MPNN


        # ----------------------------------------------------------------------------
        # Lame constants calculated from E and nu
        # ----------------------------------------------------------------------------
        E   = self.E
        nu  = self.nu
        # first Lame paramter
        la = E * nu / ( ( 1 + nu ) * ( 1 - 2 * nu ) )
        # second Lame parameter (shear modulus)
        mu = E / ( 2 + 2 * nu ) 
####

        # TODO Put it into the parameters - this is a temporary hack
        #
        if self.update_state_on:
            e_max_list = self.get_e_max_list( sctx, eps_app_eng )    
            sctx.mats_state_array[:] = e_max_list

        phi_list   = self._get_phi_list( sctx, eps_app_eng )
        
                
        # integration terms for each microplanes
        phi_vct_list = array( [ phi_list[i] * MPNN[i,:,:] * MPW[i]
                                for i in range(0,n_mp) ] )

        # integration terms for each microplanes
        psi_vct_list = array( [ 1. / phi_list[i] * MPNN[i,:,:] * MPW[i]
                                for i in range(0,n_mp) ] )

        # sum of contributions from all microplanes
        # sum over the first dimension (over the microplanes)
        phi_mtx_2D = phi_vct_list.sum(0) 
        psi_mtx_2D = psi_vct_list.sum(0) 

        phi_mtx_3D = vstack( [hstack( [phi_mtx_2D, [[0],[0]]] ), [[0,0,1]]] )
        psi_mtx_3D = vstack( [hstack( [psi_mtx_2D, [[0],[0]]] ), [[0,0,1]]] )

        phi_eig_val_2D, phi_eig_mtx_2D = eig( phi_mtx_2D )
        phi_eig_2D = array([ pe.real for pe in phi_eig_val_2D] )

        phi_eig_val_2D, phi_eig_mtx_2D = eig( phi_mtx_2D )
        phi_eig_2D = array([ pe.real for pe in phi_eig_val_2D] )
        phi_eig_3D = hstack([phi_eig_2D, array([1]) ])
        phi_eig_mtx_3D = vstack( [hstack( [phi_eig_mtx_2D, [[0],[0]]] ), [[0,0,1]]] )

        
        #### 2D case: ###-------------------------------------------------------------
        #### 2D case: ###-------------------------------------------------------------
        #### 2D case: ###-------------------------------------------------------------
        # rank-four tensor including damage effect
        D4s_mdm_2D = self.D4s_mdm_2D
        C4c_mdm_2D = self.C4c_mdm_2D
        D4s_e_2D = self.D4s_e_2D
        C4c_e_2D = self.C4c_e_2D
        # rank-two tensor (matrix) including damage effect
        D2s_mdm_2D = self.D2s_mdm_2D
        C2c_mdm_2D = self.C2c_mdm_2D
        D2s_e_2D = self.D2s_e_2D
        C2c_e_2D = self.C2c_e_2D
        

        t1 = time()
        n_d = 2
        delta = identity(n_d)
        # Fill the rank-four matrices
        for i in range(0,n_d):
            for j in range(0,n_d):
                for k in range(0,n_d):
                    for l in range(0,n_d):
                        #--------------------------------------------------------------------------------
                        # Damaged material (formula by Bazant 1997)
                        #--------------------------------------------------------------------------------
                        D4s_mdm_2D[i,j,k,l] = la * phi_mtx_2D[i,j] * phi_mtx_2D[k,l] + \
                                           mu * ( phi_mtx_2D[i,k] * phi_mtx_2D[j,l] + phi_mtx_2D[i,l] * phi_mtx_2D[j,k] )
                        #--------------------------------------------------------------------------------
                        # map the components of the fourth-order tensor to the corresponding engineering notation
                        # using minor and major symmetry (=condensation to second-order tensor, i.e. 6x6-matrix)
                        #--------------------------------------------------------------------------------
                        D2s_mdm_2D[map2d_ijkl2mn(i,j,k,l)] =  D4s_mdm_2D[i,j,k,l]
        t2 = time()
        self.t_formula2D += t2-t1
        
        print '### 2D-case: ### '
        print 'D2s_mdm_2D', D2s_mdm_2D
        print '\n'
        
        #---------------------------------------------------------------------------------------------
        # Product symmetrization using explicitly expressed beta tensor
        #---------------------------------------------------------------------------------------------
        if self.symmetrization == 'product-type':
            t1 = time()
            # this yields the same result as the applied formula for the damaged stiffness tensor D4s
            # transform the matrix phi_mtx_2D into principle damage coordinates
            # phi_pdc_mtx_2D = tensordot( tensordot( phi_eig_mtx_2D, phi_mtx_2D, [[0],[0]] ), phi_eig_mtx_2D, [[1],[0]] )
            # assure that besides the diagonal the entries are exactly zero (avoid NAN-error)

            phi_pdc_mtx_2D = zeros([2,2])
            for i in range(2):
                phi_pdc_mtx_2D[i,i] = phi_eig_2D[i]
    
            # second order damage tensor:
            w_pdc_mtx_2D = arr_sqrt( phi_pdc_mtx_2D )
            # transform the matrix w back to x-y-coordinates:
            w_mtx_2D = dot(dot(phi_eig_mtx_2D, w_pdc_mtx_2D),transpose(phi_eig_mtx_2D))
            #
            
            beta4_2D_ = outer(w_mtx_2D, w_mtx_2D).reshape(n_d,n_d,n_d,n_d)
            beta4_2D = beta4_2D_.swapaxes(1,2)

#            beta4_2D = zeros([n_d,n_d,n_d,n_d])
#            for i in range(0,n_d):
#                for j in range(0,n_d):
#                    for k in range(0,n_d):
#                        for l in range(0,n_d):
#                            beta4_2D[i,j,k,l] = w_mtx_2D[i,k] * w_mtx_2D[j,l]
            
            
        #---------------------------------------------------------------------------------------------
        # Sum-type symmetrization using explicitly expressed beta tensor
        #---------------------------------------------------------------------------------------------
        elif self.symmetrization == 'sum-type':
            # cf. Jirasek "Comments to MDM"
            beta4_2D = zeros([n_d,n_d,n_d,n_d])
            for i in range(0,n_d):
                for j in range(0,n_d):
                    for k in range(0,n_d):
                        for l in range(0,n_d):
                            beta4_2D[i,j,k,l] = 0.25 * ( phi_mtx_2D[i,k] * delta[j,l] + phi_mtx_2D[i,l] * delta[j,k] +\
                                                               phi_mtx_2D[j,k] * delta[i,l] + phi_mtx_2D[j,l] * delta[i,k] ) 

        #-----------------------------------------------------------------------------------
        # Calculate the damaged stiffness tensor based on the damage tensor beta4:
        #-----------------------------------------------------------------------------------
        #
        if self.stress_state == 'plane_stress':
            D4_e_2D = self.D4_e_2Dstress
        else:
            D4_e_2D = D4s_e_2D
        
        D4_bb_mdm_2D = tensordot( beta4_2D, tensordot( self.D4_e_2Dstrain, beta4_2D, [[2,3],[2,3]] ), [[2,3],[0,1]] )

        D2_bb_mdm_2D = zeros([3,3])
        for i in range(0,n_d):
            for j in range(0,n_d):
                for k in range(0,n_d):
                    for l in range(0,n_d):
                        D2_bb_mdm_2D[map2d_ijkl2mn(i,j,k,l)] = D4_bb_mdm_2D[i,j,k,l]

        t2 = time()
        self.t_bb2D += t2-t1
        

        print 'D2_bb_mdm_2D', D2_bb_mdm_2D
        print '\n'



        #### 3D case: ###-------------------------------------------------------------
        #### 3D case: ###-------------------------------------------------------------
        #### 3D case: ###-------------------------------------------------------------
        # rank-four tensor including damage effect
        D4s_mdm_3D = self.D4s_mdm_3D
        C4c_mdm_3D = self.C4c_mdm_3D
        D4s_e_3D = self.D4s_e_3D
        C4c_e_3D = self.C4c_e_3D
        # rank-two tensor (matrix) including damage effect
        D2s_mdm_3D = self.D2s_mdm_3D
        C2c_mdm_3D = self.C2c_mdm_3D
        D2s_e_3D = self.D2s_e_3D
        C2c_e_3D = self.C2c_e_3D

        t1 = time()
        n_d = 3
        delta = identity(n_d)
        # Fill the rank-four matrices
        for i in range(0,n_d):
            for j in range(0,n_d):
                for k in range(0,n_d):
                    for l in range(0,n_d):
                        #--------------------------------------------------------------------------------
                        # Damaged material (formula by Bazant 1997)
                        #--------------------------------------------------------------------------------
                        D4s_mdm_3D[i,j,k,l] = la * phi_mtx_3D[i,j] * phi_mtx_3D[k,l] + \
                                           mu * ( phi_mtx_3D[i,k] * phi_mtx_3D[j,l] + phi_mtx_3D[i,l] * phi_mtx_3D[j,k] )
                        #--------------------------------------------------------------------------------
                        # map the components of the fourth-order tensor to the corresponding engineering notation
                        # using minor and major symmetry (=condensation to second-order tensor, i.e. 6x6-matrix)
                        #--------------------------------------------------------------------------------
                        D2s_mdm_3D[map3d_ijkl2mn(i,j,k,l)] =  D4s_mdm_3D[i,j,k,l]
        t2=time()
        self.t_formula3D += t2-t1
        
        #---------------------------------------------------------------------------------------------
        # Product symmetrization using explicitly expressed beta tensor
        #---------------------------------------------------------------------------------------------
        if self.symmetrization == 'product-type':
            t1 = time()
            # transform the matrix phi_mtx_3D into principle damage coordinates
#            # alternative 1:
#            phi_eig_val_3D, phi_eig_mtx_3D = eig( phi_mtx_3D )
#            phi_eig_3D = array([ pe.real for pe in phi_eig_val] )
#            phi_pdc_mtx_3D = zeros([n_d,n_d])
#            for i in range(n_d):
#                phi_pdc_mtx_3D[i,i] = phi_eig_3D[i]
#            # @todo: is the order of the eigenvectors such that (0,0,1) is always the 3rd column vector?
#            #        are the vectors normalized in order to form an orthonormal basis? (Only with such an basis the
#            #        operation D = At*B*A is valid            
#            # phi_pdc_mtx_3D = tensordot( tensordot( phi_eig_mtx_3D, phi_mtx_3D, [[0],[0]] ), phi_eig_mtx_3D, [[1],[0]] )
#            # assure that besides the diagonal the entries are exactly zero (avoid NAN-error)
            # alternative 2:
            phi_pdc_mtx_3D = zeros([3,3])
            for i in range(3):
                phi_pdc_mtx_3D[i,i] = phi_eig_3D[i]
            # this assures that the third eigenvector is always (0,0,1)

            
            # second order damage tensor
            w_pdc_mtx_3D = arr_sqrt( phi_pdc_mtx_3D )
            # transform the matrix w back to x-y-coordinates:
            w_mtx_3D = dot(dot(phi_eig_mtx_3D,w_pdc_mtx_3D),transpose(phi_eig_mtx_3D))


            beta4_3D_ = outer(w_mtx_3D, w_mtx_3D).reshape(n_d,n_d,n_d,n_d)
            beta4_3D = beta4_3D_.swapaxes(1,2)
            
#            beta4_3D = zeros([n_d,n_d,n_d,n_d])
#            for i in range(0,n_d):
#                for j in range(0,n_d):
#                    for k in range(0,n_d):
#                        for l in range(0,n_d):
#                            beta4_3D[i,j,k,l] = w_mtx_3D[i,k] * w_mtx_3D[j,l]
        #---------------------------------------------------------------------------------------------
        # Sum-type symmetrization using explicitly expressed beta tensor
        #---------------------------------------------------------------------------------------------
        elif self.symmetrization == 'sum-type':
            # cf. Jirasek "Comments to MDM"
            beta4_3D = zeros([n_d,n_d,n_d,n_d])
            for i in range(0,n_d):
                for j in range(0,n_d):
                    for k in range(0,n_d):
                        for l in range(0,n_d):
                            beta4_3D[i,j,k,l] = 0.25 * ( phi_mtx_3D[i,k] * delta[j,l] + phi_mtx_3D[i,l] * delta[j,k] +\
                                                               phi_mtx_3D[j,k] * delta[i,l] + phi_mtx_3D[j,l] * delta[i,k] ) 

        #-----------------------------------------------------------------------------------
        # Calculate the damaged stiffness tensor based on the damage tensor beta4:
        #-----------------------------------------------------------------------------------
        D4_bb_mdm_3D = tensordot( beta4_3D, tensordot( self.D4_e_3D, beta4_3D, [[2,3],[2,3]] ), [[2,3],[0,1]] )

          
        
        D2_bb_mdm_3D = zeros([6,6])
        for i in range(0,n_d):
            for j in range(0,n_d):
                for k in range(0,n_d):
                    for l in range(0,n_d):
                        D2_bb_mdm_3D[map3d_ijkl2mn(i,j,k,l)] = D4_bb_mdm_3D[i,j,k,l]
                        # this yields the same result as the applied formula for the damaged stiffness tensor D4s
        

        t2 = time()
        self.t_bb3D += t2-t1


        # ------------------------------------------------------------------------------------------------
        # for debugging purposes only: if elastic_debug is switched on, linear elastic material is used 
        # ------------------------------------------------------------------------------------------------
        if self.elastic_debug:
#            C2cc_e = self._compliance_mapping3d( C2c_e )
#            D2c_e = inv( C2cc_e )
            # reduce 3d-matrix (6x6) to 2d-matrix (3x3) 
            # using assumptions for plane stress or plane strain  
            if self.stress_state == 'plane_stress':
                D2r_e = self._get_D_plane_stress( D2s_e )
            else: # plane strain
                D2r_e = self._get_D_plane_strain( D2s_e )
            sig_eng = tensordot( D2r_e, eps_app_eng, [[1],[0]])
            return sig_eng, D2r_e
        



        #-----------------------------------------------------------------------------------------------
        # Postprocess the material stiffness according to the specified configuration
        # ----------------------------------------------------------------------------------------------

        # compliance version:
#        C2cc_mdm_3D = self._compliance_mapping3d( C2c_mdm_3D )
#        D2c_mdm_3D = inv( C2cc_mdm_3D )
     
        # stiffness version:
        if self.model_version == 'stiffness' :
            D2_mdm_3D = D2s_mdm_3D
        else: # compliance
            D2_mdm_3D = D2c_mdm_3D


        # reduce 3d stiffness matrix to 2d
        if self.stress_state == 'plane_stress':
            t1 = time()
            D2r_mdm_3D = self._get_D_plane_stress( D2_mdm_3D )
            t2 = time()
            self.t_reduction += t2-t1
            D2r_bb_mdm_3D = self._get_D_plane_stress( D2_bb_mdm_3D )
        else: # plane strain
            t1 = time()
            D2r_mdm_3D = self._get_D_plane_strain( D2_mdm_3D )
            t2 = time()
            self.t_reduction += t2-t1
            D2r_bb_mdm_3D = self._get_D_plane_strain( D2_bb_mdm_3D )

        print '### 3D-case: ### '
        print 'D2r_mdm_3D', D2r_mdm_3D
        print '\n'

        print 'D2r_bb_mdm_3D', D2r_bb_mdm_3D
        print '\n'

        print '### time-evaluation : ###############'
        print 't_interim', self.t_interim
        print 't_setup_De', self.t_setup_De
        print 't_bb2D', self.t_bb2D
        print 't_bb3D', self.t_bb3D
        print 't_formula2D', self.t_formula2D
        print 't_formula3D', self.t_formula3D
        print 't_reduction', self.t_reduction
        print '### number of steps: ###' ,self.number_of_steps, '  ###  '
        print 't_setup    /t_formula2D     ', self.t_setup_De  / self.t_formula2D
        print 't_bb2D     /t_formula2D     ', self.t_bb2D      / self.t_formula2D 
        print 't_bb3D     /t_formula2D     ', self.t_bb3D      / self.t_formula2D 
        print 't_formula2D/t_formula2D     ', self.t_formula2D / self.t_formula2D 
        print 't_formula3D/t_formula2D     ', self.t_formula3D / self.t_formula2D 
        print 't_reduction/t_formula2D     ', self.t_reduction / self.t_formula2D 
        print '#####################################'
                   
        
        sig_eng = tensordot( D2r_mdm_3D, eps_app_eng, [[1],[0]])
        return sig_eng, D2r_mdm_3D

    #---------------------------------------------------------------------------------------------
    # Subsidiary methods realizing configurable features
    #---------------------------------------------------------------------------------------------
  
    def _compliance_mapping2d( self, C2_mtx ):
        '''
        Compliance matrix must reflect the symmetry assumption 
        The components of the tensor are multiplied with factor 1,2 or 4 depending on their 
        position in the matrix (due to symmetry and switching of tensorial to engineering notation
        (Note: gamma_xy = 2*epsilon_xy, etc.) 
        '''
        idx1 = [0,1]
        idx2 = [2]
        
        C11 = C2_mtx[ix_(idx1,idx1)]
        C12 = C2_mtx[ix_(idx1,idx2)]
        C21 = C2_mtx[ix_(idx2,idx1)]
        C22 = C2_mtx[ix_(idx2,idx2)]

        return vstack( [ hstack( [C11,   2*C12] ),
                         hstack( [2*C21, 4*C22] ) ] )


    def _compliance_mapping3d( self, C2_mtx ):
        '''
        Compliance matrix must reflect the symmetry assumption 
        The components of the tensor are multiplied with factor 1,2 or 4 depending on their 
        position in the matrix (due to symmetry and switching of tensorial to engineering notation
        (Note: gamma_xy = 2*epsilon_xy, etc.) 
        '''
        idx1 = [0,1,2]
        idx2 = [3,4,5]
        
        C11 = C2_mtx[ix_(idx1,idx1)]
        C12 = C2_mtx[ix_(idx1,idx2)]
        C21 = C2_mtx[ix_(idx2,idx1)]
        C22 = C2_mtx[ix_(idx2,idx2)]

        return vstack( [ hstack( [C11,   2*C12] ),
                         hstack( [2*C21, 4*C22] ) ] )

    def _get_D_plane_stress( self, D2_mtx ):
        idx2 = [0,1,5]
        idx3 = [2,3,4]
        D22 = D2_mtx[ix_(idx2,idx2)]
        D23 = D2_mtx[ix_(idx2,idx3)]
        D32 = D2_mtx[ix_(idx3,idx2)]
        D33 = D2_mtx[ix_(idx3,idx3)]
        D2_pstress_term = dot( dot( D23, inv( D33 ) ), D32 )
        D2r_mdm = D22 - D2_pstress_term
        return D2r_mdm

    def _get_D_plane_strain( self, D2_mtx ):
        idx2 = [0,1,5]
        D22 = D2_mtx[ix_(idx2,idx2)]
        return D22

    def _get_C_plane_strain( self, C2_mtx ):
        idx2 = [0,1,5]
        idx3 = [2,3,4]
        C22 = C2_mtx[ix_(idx2,idx2)]
        C23 = C2_mtx[ix_(idx2,idx3)]
        C32 = C2_mtx[ix_(idx3,idx2)]
        C33 = C2_mtx[ix_(idx3,idx3)]
        C2_pstrain_term = dot( dot( C23, inv( C33 ) ), C32 )
        C2r_mdm = C22 - C2_pstrain_term
        return C2r_mdm

    def _get_C_plane_stress( self, C2_mtx ):
        idx2 = [0,1,5]
        C22 = C2_mtx[ix_(idx2,idx2)]
        return C22


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
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0 )
        return sig_eng

    def get_sig_norm( self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0 )
        return array( [ scalar_sqrt( sig_eng[0]**2 + sig_eng[1]**2 ) ] )
    
    def get_microplane_damage(self, sctx, eps_app_eng ):
        phi_list = self._get_phi_list(sctx, eps_app_eng)
        return phi_list

    def get_fracture_energy_list(self, sctx, eps_app_eng):
        '''
        Get the microplane contributions to the fracture energy
        '''
        e_max_list = self.get_e_max_list(sctx, eps_app_eng)
        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
        return fracture_energy_list

    def get_fracture_energy(self, sctx, eps_app_eng):
        '''
        Get the macroscopic fracture energy as a weighted sum of all mircoplane contributions
        '''
        e_max_list = self.get_e_max_list(sctx, eps_app_eng)
        fracture_energy_list = self.polar_fn.get_fracture_energy_list( e_max_list )
        return array( [dot( self._MPW, fracture_energy_list )], float )
    
    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return { 'eps_app'           :  self.get_eps_app,
                 'sig_app'           :  self.get_sig_app,
                 'sig_norm'          :  self.get_sig_norm,
                 'microplane_damage' :  self.get_microplane_damage,
                 'fracture_energy_list'   :  self.get_fracture_energy_list,
                 'fracture_energy'   :  self.get_fracture_energy }

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


def map3d_tns2_to_tns4(tns2):
    '''map a matrix to a fourth order tensor assuming minor and major symmetry,
    e.g. D_mtx (6x6) in engineering notation to D_tns(3,3,3,3)).
    '''
    n_d = 3
    n_eng = 6
    tns4 = zeros([n_d,n_d,n_d,n_d])
    for i in range(0,n_d):
        for j in range(0,n_d):
            for k in range(0,n_d):
                for l in range(0,n_d):
                    tns4[i,j,k,l] = tns2[map3d_ijkl2mn(i,j,k,l)]
    return tns4
               
                
def map3d_tns4_to_tns2(tns4):
    '''map a fourth order tensor to a matrix assuming minor and major symmetry,
    e.g. D_tns(3,3,3,3) to D_mtx (6x6) in engineering notation.
    '''
    n_d = 3
    n_eng = 6
    tns2 = zeros([n_eng,n_eng])
    for i in range(0,n_d):
        for j in range(0,n_d):
            for k in range(0,n_d):
                for l in range(0,n_d):
                    tns2[map3d_ijkl2mn(i,j,k,l)] = tns4[i,j,k,l]
    return tns2


def map2d_tns2_to_tns4(tns2):
    '''map a matrix to a fourth order tensor assuming minor and major symmetry,
    e.g. D_mtx (3x3) in engineering notation to D_tns(2,2,2,2)).
    '''
    n_d = 2
    n_eng = 3
    tns4 = zeros([n_d,n_d,n_d,n_d])
    for i in range(0,n_d):
        for j in range(0,n_d):
            for k in range(0,n_d):
                for l in range(0,n_d):
                    tns4[i,j,k,l] = tns2[map2d_ijkl2mn(i,j,k,l)]
    return tns4
               
                
def map2d_tns4_to_tns2(tns4):
    '''map a fourth order tensor to a matrix assuming minor and major symmetry,
    e.g. D_tns(2,2,2,2) to D_mtx (3x3) in engineering notation.
    '''
    n_d = 2
    n_eng = 3
    tns2 = zeros([n_eng,n_eng])
    for i in range(0,n_d):
        for j in range(0,n_d):
            for k in range(0,n_d):
                for l in range(0,n_d):
                    tns2[map2d_ijkl2mn(i,j,k,l)] = tns4[i,j,k,l]
    return tns2


# @todo - temporary alias rename the class and test it all
MA2DCompositeMicroplaneDamage = MATS2DMicroplaneDamage

#--------------------------------------------------------------------------------
# Example 
#--------------------------------------------------------------------------------

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

from core.time_loop import TimeLoop as TL, DACRangeValue as TRange
from core.bc import BCDof
from core.ibvp_solve import IBVPSolve as IS

class MATS2DMicroplaneDamageExplore( HasTraits ):

    alpha_degree = Range( 0., 360., 0., 
                          label = 'Loading angle',
                          auto_set = False)

    ### 
    bc_alpha = Instance( BCDof )
    def _bc_alpha_default( self ):
        alpha =  Pi * self.alpha_degree / 180
        value, coeff = get_value_and_coeff( 1., alpha )
        return  BCDof(var='u', dof = 0, value = value,
                      link_dofs = [1],
                      link_coeffs = [coeff],
                      time_function = lambda t: t )

    @on_trait_change('alpha_degree')
    def update_bc_alpha( self ):
        alpha =  Pi * self.alpha_degree / 180
        value, coeff = get_value_and_coeff( 1., alpha )
        self.bc_alpha.value = value
        self.bc_alpha.link_coeffs[0] = coeff

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
                 bc_list = [ self.bc_alpha 
                             ],
                 rv_list = [ RTOGraph(name = 'strain 0 - stress 0',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 0,
                                      record_on = 'update' ),
                             RTOGraph(name = 'strain 1 - stress 1',
                                      var_x = 'eps_app', idx_x = 1,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTOGraph(name = 'strain 0 - stress 1',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTOGraph(name = 'strain 1 - stress 0',
                                      var_x = 'eps_app', idx_x = 1,
                                      var_y = 'sig_app', idx_y = 0,
                                      record_on = 'update' ),
                             RTOGraph(name = 'strain 0 - strain 1',
                                      var_x = 'eps_app', idx_x = 0,
                                      var_y = 'eps_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTOGraph(name = 'stress 0 - stress 1',
                                      var_x = 'sig_app', idx_x = 0,
                                      var_y = 'sig_app', idx_y = 1,
                                      record_on = 'update' ),
                             RTOGraph(name = 'time - sig_norm',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'sig_norm', idx_y = 0,
                                      record_on = 'update' ),
                             RTOGraph(name = 'time - G_f',
                                      var_x = 'time', idx_x = 0,
                                      var_y = 'fracture_energy', idx_y = 0,
                                      record_on = 'update' ),
                             RTOArraySnapshot(name = 'fracture energy contributions',
                                      var = 'fracture_energy_list',
                                      record_on = 'update' ),
                             RTOArraySnapshot(name = 'microplane damage',
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
            tmax = 0.0001
            # tmax = 0.0006
            n_steps = 60

        tl = TL( ts = ts,
                 DT=tmax/n_steps, KMAX = 100, RESETMAX = 0,
                 T = TRange( min = 0.0,  max = tmax ) )

        return IS( time_loop = tl )

    traits_view = View( Item('alpha_degree@', label = 'angle'),
                        Item('bc_alpha@'),
                        Item('sim@'),
                        resizable = True,
                        width = 1.0,
                        height = 1.0
                        )

    traits_view_begehung = View( Item('alpha_degree@', label = 'angle'),
                                 Item('sim'),
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

    bc_alpha = BCDof(var='u', dof = 0, value = value,
                     link_dofs = [1],
                     link_coeffs = [coeff],
                     time_function = lambda t: t )
    
    ts = TS( tse = tseval,
             bc_list = [ bc_alpha 
                         ],
             rv_list = [ RTOGraph(name = 'strain 0 - stress 0',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 0,
                                  record_on = 'update' ),
                         RTOGraph(name = 'strain 1 - stress 1',
                                  var_x = 'eps_app', idx_x = 1,
                                  var_y = 'sig_app', idx_y = 1,
                                  record_on = 'update' ),
                         RTOGraph(name = 'strain 0 - stress 1',
                                  var_x = 'eps_app', idx_x = 0,
                                  var_y = 'sig_app', idx_y = 1,
                                  record_on = 'update' ),
                         RTOGraph(name = 'strain 1 - stress 0',
                                  var_x = 'eps_app', idx_x = 1,
                                  var_y = 'sig_app', idx_y = 0,
                                  record_on = 'update' ),
                         RTOGraph(name = 'strain 0 - strain 1',
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
        tmax = 0.0001
        # tmax = 0.0006
        n_steps = 60

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
        bc_alpha.value = value
        bc_alpha.link_coeffs[0] = coeff

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
