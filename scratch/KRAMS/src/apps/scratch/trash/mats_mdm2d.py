
from enthought.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, \
     Instance, Int, Trait, Range, HasTraits, on_trait_change, Event, \
     implements, Dict, Property, cached_property, Delegate

from enthought.traits.ui.api import \
     Item, View, HSplit, VSplit, VGroup, Group, Spring

# Chaco imports
from enthought.chaco.chaco_plot_editor import \
     ChacoPlotEditor, \
     ChacoPlotItem
from enthought.enable.component_editor import \
     ComponentEditor
from enthought.chaco.tools.api import \
     PanTool, SimpleZoom
from enthought.chaco.api import \
     Plot, AbstractPlotData, ArrayPlotData

#from dacwt import DAC

from numpy import \
     array, ones, zeros, outer, inner, transpose, dot, frompyfunc, \
     fabs, sqrt, linspace, vdot, identity, tensordot, \
     sin as nsin, meshgrid, float_, ix_, \
     vstack, hstack, sqrt as arr_sqrt

from math import pi as Pi, cos, sin, exp, sqrt as scalar_sqrt

from scipy.linalg import eig, inv

from ibvpy.core.tstepper import \
     TStepper as TS

from ibvpy.mats.mats_eval import IMATSEval, MATSEval

from api import RTrace, RTraceGraph, RTraceArraySnapshot
from mats_mdm2d_phi import MPArrayMDM, MPArrayCMDM, IMPArray
#---------------------------------------------------------------------------
# Material time-step-evaluator for Microplane-Damage-Model
#---------------------------------------------------------------------------

class MACMDM( MATSEval ):
    '''
    Microplane Damage Model.
    '''

    implements( IMATSEval )

    #---------------------------------------------------------------------------
    # Parameters of the numerical algorithm (integration)
    #---------------------------------------------------------------------------

    n_mp = Delegate('mparray')
    
    model_version = Enum("compliance","stiffness")
    stress_state  = Enum("plane_stress","plane_strain")

    n_d = Trait( '2D', {'2D' : 2, '3D' : 3 },
                 label = 'Dimensionality',
                 desc = 'This value is inactive yet',
                 auto_set = False)
    #n_d = Int(2)

    @on_trait_change('n_mp')
    def update_arrays( self ):
        self._init_arrays()
        self.changed = True
        
    #---------------------------------------------------------------------------
    # Material parameters 
    #---------------------------------------------------------------------------

    E   = Float( 34e+3,
                 label = "E",
                 desc = "Young's Modulus",
                 auto_set = False )
    nu  = Float( 0.2,
                 label = 'nu',
                 desc = "Poison's ratio",
                 auto_set = False )
    c_T = Float( 0.0,
                 desc = 'fraction of tangential stress accounted on each microplane',
                 auto_set = False)

    damage_function = Trait('IsotropicQuasiBrittle',
                            {'IsotropicQuasiBrittle' : MPArrayMDM,
                             'AnisotropicQuasiDuctile' : MPArrayCMDM } )

    mparray = Property( Instance( IMPArray ), depends_on = 'damage_function' )
    @cached_property
    def _get_mparray(self):
        return self.damage_function_()

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

    view_traits = View( VSplit( Group(Item('E'),
                                      'nu','c_T',
                                      Item('damage_function@'),
                                      Item('mparray', style='custom', show_label = False),
                                      label='Material parameters',
                                      show_border=True),
                                Group( Item('model_version', style = 'custom' ),
                                       Item('stress_state', style = 'custom' ),
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
        super( MACMDM, self ).__init__( **kwtraits )
        self._init_arrays()

    alpha_list = Delegate( 'mparray' )

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
#        sctx.state_array[:] = 0.0
        ## ?TODO change by alex
        state_arr_size = self.get_state_array_size()
        sctx.state_array = zeros(state_arr_size, 'float_')
        ##
        self.update_state_on = False
        self._setup_interim_arrays()

    def _setup_interim_arrays( self ):
        '''
        Intialize intermediate arrays
        '''
        n_d = 3
        # rank-four tensor including damage effect
        self.D4s_mdm = zeros([n_d,n_d,n_d,n_d])
        self.C4c_mdm = zeros([n_d,n_d,n_d,n_d])
        self.D4s_e = zeros([n_d,n_d,n_d,n_d])
        self.C4c_e = zeros([n_d,n_d,n_d,n_d])
        # rank-two tensor (matrix) including damage effect
        self.D2s_mdm = zeros([6,6])
        self.C2c_mdm = zeros([6,6])
        self.D2s_e = zeros([6,6])
        self.C2c_e = zeros([6,6])

    def new_cntl_var(self):
        return zeros( 3, float_ )

    def new_resp_var(self):
        return zeros( 3, float_ )

    def _get_phi_list(self, sctx, eps_app_eng ):

        # Prepare the microplane coefficients needed for the 
        # projection and integration
        #
        MPN  = self._MPN
        MPW  = self._MPW
        MPNN = self._MPNN

        # tangent strain ratio
        #
        c_T = self.c_T
        n_mp = self.n_mp

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
        
        return self.mparray.get_phi_list( e_equiv_list, sctx )
        
    #-----------------------------------------------------------------------------------------------
    # Evaluation - get the corrector and predictor
    #-----------------------------------------------------------------------------------------------

    def get_corr_pred( self, sctx, eps_app_eng, tn, tn1 ):
        '''
        Corrector predictor computation.
        @param eps_app_eng input variable - engineering strain
        '''
        #n_d = self.n_d
        n_mp = self.n_mp

        E   = self.E
        nu  = self.nu

        phi_list, e_max_list = self._get_phi_list( sctx, eps_app_eng )
        
        MPN  = self._MPN
        MPW  = self._MPW
        MPNN = self._MPNN

        # TODO Put it into the parameters - this is a temporary hack
        #
        if self.update_state_on:
            sctx.state_array[:] = e_max_list
            self.update_state_on = False
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
                        # Damaged material
                        #--------------------------------------------------------------------------------
                        D4s_mdm[i,j,k,l] = la * phi_mtx[i,j] * phi_mtx[k,l] + \
                                           mu * ( phi_mtx[i,k] * phi_mtx[j,l] + phi_mtx[i,l] * phi_mtx[j,k] )
                        C4c_mdm[i,j,k,l] = (1+nu)/(2*E) * \
                                           ( psi_mtx[i,k] * psi_mtx[j,l] + psi_mtx[i,l]* psi_mtx[j,k] ) - \
                                           nu / E * psi_mtx[i,j] * psi_mtx[k,l]
                        D2s_mdm[map3d_ijkl2mn(i,j,k,l)] =  D4s_mdm[i,j,k,l]
                        C2c_mdm[map3d_ijkl2mn(i,j,k,l)] =  C4c_mdm[i,j,k,l]
                        #--------------------------------------------------------------------------------
                        # Elastic material
                        #--------------------------------------------------------------------------------
                        D4s_e[i,j,k,l] = la * delta[i,j] * delta[k,l] + \
                                        mu * ( delta[i,k] * delta[j,l] + delta[i,l] * delta[j,k] )
                        C4c_e[i,j,k,l] = (1+nu)/(2*E) * \
                                         ( delta[i,k] * delta[j,l] + delta[i,l]* delta[j,k] ) - \
                                         nu / E * delta[i,j] * delta[k,l]
                        D2s_e[map3d_ijkl2mn(i,j,k,l)] =  D4s_e[i,j,k,l]
                        C2c_e[map3d_ijkl2mn(i,j,k,l)] =  C4c_e[i,j,k,l]

        if self.elastic_debug:
            C2cc_e = self._compliance_mapping( C2c_e )
            D2c_e = inv( C2cc_e )
            if self.stress_state == 'plane_stress':
                D2r_e = self._get_D_plane_stress( D2s_e )
            else: # plane strain
                D2r_e = self._get_D_plane_strain( D2s_e )
            sig_eng = tensordot( D2r_e, eps_app_eng, [[1],[0]])
            return sig_eng, D2r_e

        C2cc_mdm = self._compliance_mapping( C2c_mdm )
        D2c_mdm = inv( C2cc_mdm )

        #---------------------------------------------------------------------------------------------
        # Product symmetrization using explicitly expressed beta tensor
        #---------------------------------------------------------------------------------------------
        phi_eig_val, phi_eig_mtx = eig( phi_mtx )
        phi_eig = array([ pe.real for pe in phi_eig_val] )

        # verify the transformation
        #phi_pdc_mtx = tensordot( tensordot( phi_eig_mtx, phi_mtx, [[0],[0]] ), phi_eig_mtx, [[1],[0]] )
        phi_pdc_mtx = tensordot( tensordot( phi_mtx, phi_eig_mtx, [[0],[0]] ), phi_eig_mtx, [[1],[0]] )
        w_mtx = arr_sqrt( phi_pdc_mtx )

        beta4_tns = zeros([n_d,n_d,n_d,n_d])
        for i in range(0,n_d):
            for j in range(0,n_d):
                for k in range(0,n_d):
                    for l in range(0,n_d):
                        beta4_tns[i,j,k,l] = w_mtx[i,k] * w_mtx[j,l]

        D4_bb_mdm = tensordot( tensordot( D4s_e, beta4_tns, [[0,1],[2,3]] ), beta4_tns, [[2,3],[2,3]] )

        # Postprocess the material stiffness according to the
        # specified configuration
        #
        if self.model_version == 'stiffness' :
            D2_mdm = D2s_mdm
        else: # compliance
            D2_mdm = D2c_mdm

        if self.stress_state == 'plane_stress':
            D2r_mdm = self._get_D_plane_stress( D2_mdm )
        else: # plane strain
            D2r_mdm = self._get_D_plane_strain( D2_mdm )

        sig_eng = tensordot( D2r_mdm, eps_app_eng, [[1],[0]])
        return sig_eng, D2r_mdm

    #---------------------------------------------------------------------------------------------
    # Subsidiary methods realizing configurable features
    #---------------------------------------------------------------------------------------------

    def _compliance_mapping( self, C2_mtx ):
        '''
        Compliance matrix must reflect the symmetry assumption 
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

    #---------------------------------------------------------------------------------------------
    # Update state method called upon an accepted time-step
    #---------------------------------------------------------------------------------------------

    def update_state(self, sctx, eps_app_eng ):
        '''
        Here just set the flag on to make the update afterwords in the method itself.
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
        # cached.
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0 )
        return sig_eng

    def get_sig_norm( self, sctx, eps_app_eng ):
        sig_eng, D_mtx = self.get_corr_pred( sctx, eps_app_eng, 0, 0 )
        return array( [ scalar_sqrt( sig_eng[0]**2 + sig_eng[1]**2 ) ] )
    
    def get_microplane_damage(self, sctx, eps_app_eng ):
        phi_list, e_max_list = self._get_phi_list(sctx, eps_app_eng)
        return phi_list

    # Declare and fill-in the rte_dict - it is used by the clients to
    # assemble all the available time-steppers.
    #
    rte_dict = Trait( Dict )
    def _rte_dict_default(self):
        return { 'eps_app' : self.get_eps_app,
                 'sig_app' : self.get_sig_app,
                 'sig_norm' : self.get_sig_norm,
                 'microplane_damage' : self.get_microplane_damage }

#---------------------------------------------------------------------------------------------
# Subsidiary index mapping functions for rank-four to rank-two tensors
#---------------------------------------------------------------------------------------------

def map_ijkl2mn(i,j,k,l):
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

# @todo - temporary alias rename the class and test it all
MA2DCompositeMicroplaneDamage = MACMDM

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

from ibvpy.core.tloop import TLoop, TLine
from ibvpy.api import BCDof
from ibvpy.core.ibvp_solve import IBVPSolve as IS

class MACMDMExplore( HasTraits ):

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
        elastic_debug = False
        # tseval for a material model
        #
        tseval  = MACMDM( elastic_debug = elastic_debug )
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
    tseval  = MACMDM( elastic_debug = elastic_debug )


    value, coeff = get_value_and_coeff( 1., 0.0 )

    bcond_alpha = BCDof(var='u', dof = 0, value = value,
                     link_dofs = [1],
                     link_coeffs = [coeff],
                     time_function = lambda t: t )
    
    ts = TS( tse = tseval,
             bcond_list = [ bcond_alpha 
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
        n_steps = 60

    tl = TLoop( tstepper = ts,
             DT=tmax/n_steps, KMAX = 100, RESETMAX = 0,
             tline = TLine( min = 0.0,  max = tmax ) )

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

        eps0_sig0 = tl.rtrace_mngr.rtrace_list[0]
        eps1_sig1 = tl.rtrace_mngr.rtrace_list[1]

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

    from mfn_line import MFnLineArray
    
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
    mme = MACMDMExplore()
    #mme.configure_traits( view = 'traits_view_begehung' )
    mme.configure_traits()
