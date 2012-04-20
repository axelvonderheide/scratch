
from enthought.traits.api import \
     Array, Bool, Enum, Float, HasTraits, \
     HasStrictTraits, \
     Instance, Int, Trait, Str, Enum, \
     Callable, List, TraitDict, Any, Range, \
     Delegate, Event, on_trait_change, Button, \
     Interface, Property, cached_property, WeakRef, Dict

from enthought.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, \
     HSplit, Group, Handler, VSplit

from enthought.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, \
     Action

from numpy import zeros, float_

from ibvpy.api import\
    ISDomain, SDomain, SContext, RTraceMngr, ITStepperEval, TStepperEval, TStepper
#from sdomain import SDomain
#from scontext import SContext

from mgrid_domain import MeshManager
from dots_eval_enr import DOTSManager

from ibvpy.core.bc_mngr import BCondMngr
#from rtrace_mngr import RTraceMngr
from ibvpy.core.ibv_resource import IBVResource
#from tstepper_eval import ITStepperEval, TStepperEval

class TStepper_enr(TStepper):
    """
    The TStepper is a spatially bounded TStepperEval.

    @TODO update this for the new class names
    
    The binding is done by associating the time-stepper with a spatial
    object. In fact, any TStepper, not only the top-level one must
    be associated with spatial object. For the chained
    sub-time-steppers, this association is done using a parameter sctx
    (see the specification above). This parameters stands for spatial
    context. Note, the distinction between the spatial object and
    spatial context. While the spatial object represents a single
    spatial entity (one of domain, element, layer or material point)
    the spatial context represents a complex reference within the
    spatial object In particular, spatial context of a particular
    material point is represented as tuple containing tuple of
    references to [domain, element, layer, integration cell, material
    point].
    """

    # service specifiers - used to link the service to this object
    service_class = 'ibvpy.plugins.tstepper_service.TStepperService'
    service_attrib = 'tstepper'

    # Sub-time-stepper or integrator.
    #
    tse = Instance( DOTSManager )

    # Spatial domain to bind the time-stepper to.
    #
    sdomain = Instance( MeshManager )
#    def _sdomain_default( self ):
#        return SDomain()

    state_array = Array

    sctx = Instance( SContext )
    
    # Boundary condition manager
    #
    bcond_mngr = Instance( BCondMngr )
    def _bcond_mngr_default( self ):
        return BCondMngr()

    # Convenience constructor 
    #
    # This property provides the possibility to write
    # tstepper.bcond_list = [BCDof(var='u',dof=5,value=0, ... ]
    # The result gets propageted to the BCondMngr
    #
    bcond_list = Property( List )
    def _set_bcond_list( self, bcond_list ):
        self.bcond_mngr.bcond_list = bcond_list

    # Response variable manager
    #
    rtrace_mngr = Instance( RTraceMngr )
    def _rtrace_mngr_default( self ):
        return RTraceMngr( tstepper = self )

    # Convenience constructor 
    #
    # This property provides the possibility to write
    # tstepper.bcond_list = [RVDof(var='u',dof=5,value=0, ... ]
    # The result gets propageted to the RTraceMngr
    #
    rtrace_list = Property( List )
    def _set_rtrace_list( self, rtrace_list ):
        self.rtrace_mngr.rtrace_list = rtrace_list

    # Backward reference to the time-loop in order to accommadate the
    # response-trace-evaluators from the top level. These include
    # bookkeeping data like memory usage or solving time per selected
    # type of operation.
    #
    tloop = WeakRef

    rte_dict = Property( Dict, depends_on = 'tse:rte_dict')
    def _get_rte_dict( self ):
        '''
        Gather all the currently applicable evaluators from the sub-ts
        and from the time-loop.

        Note the control data (time-loop data) is available within the
        model to construct flexible views (tracers) on the model.
        '''
        _rte_dict = {}
        _rte_dict.update( self.tse.rte_dict )
        if self.tloop:
            _rte_dict.update( self.tloop.rte_dict )
        return _rte_dict

    new_cntl_var = Delegate( 'tse' )

    new_resp_var = Delegate( 'tse' )

    # Vector of external forces
    #
    F_ext = Array

    # missing - setup of the time-stepper itself. reacting to changes
    # in the sub time-steppers. bcond_list and rtrace_list must be reset once
    # a change has been performed either in a spatial domain or in
    # tse.
    #
    def setup( self ):

        # Put the spatial domain into the spatial context
        #
        self.sctx = sctx = self.sdomain.new_scontext()

        # Let the boundary conditions setup themselves within the
        # spatial context
        #
        # TODO - the setup needs the link to the algorithm and to the
        # time-steppers as well.!
        #
        self.bcond_mngr.setup( sctx )

        # Let the response variables setup themselves within the
        # spatial context
        #
        self.rtrace_mngr.setup( sctx )

        # Get the size of the state state and set it up
        #
        sarr_size = self.tse.get_state_array_size( )
        self.state_array = zeros( sarr_size, float_ )
        sctx.state_array = self.state_array

        # Let the time-step-evaluators setup themselves within spatial context
        #
        # This steps includes initialization of state arrays and
        # eventual DOF mappings.
        #
        self.tse.setup( sctx )

        # Return the response variable to be used when assembling the
        # boundary conditions. Should the bcond_mngr take care of this? 
        # That's the source object, isn't it? BCondMngr is the bounded
        # version of the conditions, it could supply the matrix
        # autonomously.
        #
        self.F_ext = self.new_resp_var()
        
        # Prepare the global update flag
        sctx.update_state_on = False

    def eval( self, step_flag, U_k, d_U, t_n, t_n1 ):

        # Put the spatial domain into the spatial context
        #
        sctx = self.sctx

        # Let the time sub-stepper evaluate its contribution.
        #
        F_int, K = self.tse.get_corr_pred( sctx, U_k, d_U, t_n, t_n1 )
        
        #Switch off the global update flag
        sctx.update_state_on = False

        self.F_ext[:] = 0.

        # Apply the boundary conditions
        #
        self.bcond_mngr.apply( step_flag, K, F_int, self.F_ext, t_n, t_n1 )

        return K, F_int, self.F_ext
    
    def update_state(self, U):
        '''
        spatial context represents a stack with the top object
         representing the current level.
        @param U: 
        '''
        #sctx = ( self.sdomain, )
        self.sctx.update_state_on = True
        #self.tse.update_state( sctx, U )

    def register_mv_pipelines(self,e):
        '''Register the visualization pipelines in mayavi engine
        '''
        self.tse.register_mv_pipelines(e)
        scene = e.new_scene()
        scene.name = 'Spatial domain'
        self.sdomain.register_mv_pipelines(e)
        self.rtrace_mngr.register_mv_pipelines(e)

    traits_view = View( Group( Item('sdomain', style = 'custom', show_label = False),
                               label = 'Discretization'),
                        Group( Item( 'tse', style = 'custom', show_label = False),
                               label = 'Integrator'),
                        Group( Item('bcond_mngr', style = 'custom', show_label = False),
                               label = 'Boundary conditions'),
                        resizable = True,
                        height = 0.8,
                        width = 0.8,
                        buttons = [OKButton,CancelButton],
                        kind = 'subpanel',
                        )
