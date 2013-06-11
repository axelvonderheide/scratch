'''
Generic implementation of the time loop.
'''

from enthought.traits.api import Array, Bool, Enum, Float, HasTraits, \
                                 Instance, Int, Trait, Str, Enum, \
                                 Callable, List, TraitDict, Any, Range, \
                                 Delegate, Event, on_trait_change, Button, Property, \
                                 cached_property
from enthought.traits.ui.api import Item, View, HGroup, ListEditor, VGroup, \
     HSplit, Group, Handler, VSplit
from enthought.traits.ui.menu import NoButtons, OKButton, CancelButton, \
     Action
     
from core.tloop import TLoop

from weakref import ref

from math import pow, fabs
from numpy import copy, array
from scipy.linalg import norm, solve

from core.rv import RTraceMngr
from core.astrategy import AStrategyBase

from threading import Thread
import time
from scipy import linsolve, sparse

#import pymem
from core.tstepper import TStepper as TS

LOGGING_ON = False

if LOGGING_ON:
    import logging
    log = logging.getLogger('debug')

class DACRangeValue(HasTraits):
    min = Float( 0.0 )
    max = Float( 1.0 )
    val = Float( 0.0 )

class CompTimer(object):
    def __init__(self):
        self.last_duration  = 0.0
        self.duration  = 0.0
        self.reset()

    def reset(self):
        self.start     = time.clock() # time of this process
        self.start_all = time.time()  # wall clock time

    def record(self):
        self.last_duration = self.duration
        self.duration = time.clock() - self.start
        self.overall_duration = time.time() - self.start_all

    def report(self):
        print "Computation Time : %8.2f sec" % self.overall_duration 
        print "Elapsed time     : %8.2f sec" % self.duration

class TLoopHandler( Handler ):

    computation_thread = Instance(Thread)
    
    def setattr( self, info, object, name, value ):
        Handler.setattr( self, info, object, name, value )
        info.object._updated += 1

    def object__updated_changed( self, info ):
        if info.initialized:
            info.ui.title += '*'

    def recalculate( self, info ):
        '''
        Eval in a thread.
        '''
        if self.computation_thread and self.computation_thread.isAlive():
            info.object.user_wants_abort = True
        else:
            info.object.user_wants_abort = False
            self.computation_thread = Thread(target=info.object.eval, name="computation")
            self.computation_thread.start()

RecalcAction = Action( name = 'Recalculate', action = 'recalculate' )

from enthought.traits.ui.api import TreeNodeObject

class STLoop(TLoop):

    user_wants_abort = False

    name     = Str('<unnamed>')

    DT       = Float(1.0)
    tline        = Trait( DACRangeValue() )
    KMAX     = Int(40)
    RESETMAX = Int(10)
    tolerance     = Float( 1e-8 )

    # 
    #
    adap     = Trait( AStrategyBase() )
    ts       = Trait( TS )
    rtrace_mngr  = Delegate('ts')

    # Supply the global resp-trace-evals. They can be accessed from
    # the time-stepper and included in views to the model.
    #
    rte_dict = Property( List, depends_on = 'ts:rte_dict' )
    @cached_property
    def _get_rte_dict(self):
        return { 'time'        : lambda sctx, U_k: array([self.t_n1]) ,
                 'mem_counter' : lambda sctx, U_k: self.mem_counter,
                 'U_k'         : lambda sctx, U_k: self.U_k,
                 'F_int'       : lambda sctx, U_k: self.F_int,
                 'F_ext'       : lambda sctx, U_k: self.F_ext }

    sync_resp_tracing = Bool( False )
    
    def __init__(self,*args,**kwtraits):
        super(TLoop,self).__init__(**kwtraits)
        self._updated = 0

        # make algorithm accessible from the components
        #
        self.adap.alg = ref(self)
        self.ts.tloop = self

        # Computation timer for individual steps.
        #
        self.eval_timer = CompTimer()
        self.iter_timer = CompTimer()
        self.crpr_timer = CompTimer()
        self.solv_timer = CompTimer()
        self.rtrace_mngr_timer = CompTimer()
        #self.mem_counter = pymem.MemCounter()

    def setup(self):
        self.t_n   = self.T.min
        self.t_n1  = self.t_n
        self.d_t   = self.DT
        self.k     = 0
        self.n_reset = 0
        self.ttot_k = 0

        self.ts.setup()

        # Initialize the variables
        #
        self.R_k  = self.ts.new_resp_var()
        self.U_k  = self.ts.new_cntl_var()
        self.U_n  = self.ts.new_cntl_var()

        self.U_k[:] = 0.0
        self.tot_k  = 0.0
        self.ttot_k = 0.0
        self.t_n   = self.T.min
        self.t_n1  = self.t_n

    def eval( self, e = None ):

        if not self.sync_resp_tracing:
            self.rtrace_mngr.start_timer()

        adap = self.adap
        tstepper = self.ts
        self.setup()

        self.ls_counter = 0
        
        # Measure computation time
        
        self.eval_timer.reset()

        # Register the list of variables requrired by the adaptive
        # strategy in the time stepper to include their evaluation in
        # the corrector-predictor evaluation ts.extend_output_list(
        # adap.get_indicator_list() )
        #
        while self.t_n1 <= self.T.max and not self.user_wants_abort:

            self.iter_timer.reset()

            if LOGGING_ON:
                log.info( '===================================')
                log.info('time t_n1 = %f',self.t_n1)

            self.ls_counter += 1
            self.report_load_step()
            self.k = 0
            step_flag = 'predictor'
            while self.k < self.KMAX:
                
                self.report_iteration()

                # Extract the matrix representation of the time t_n1
                # stepper for the given time and control variable u
                #
                self.crpr_timer.reset()
                self.K, self.F_int, self.F_ext = tstepper.eval( step_flag,
                                                                    self.U_k,
                                                                    self.t_n,
                                                                    self.t_n1 )
                
                #self.K = sparse.csr_matrix(K_full)
                self.crpr_timer.record()

                self.rtrace_mngr_timer.reset()
                self.rtrace_mngr.record_iter( self.U_k )
                self.rtrace_mngr_timer.record()

                if adap.ihandler_needed():
                    adap.ihandler_invoke()
                    self.d_t *= adap.ihandler_get_scale()
                    self.t_n1 = self.t_n + self.d_t
                    step_flag = 'predictor'
                    self.k = 0

                    if LOGGING_ON:
                        log.info( "redo time step" )
                    continue

                self.R_k = self.F_ext - self.F_int

                if norm(self.R_k) < self.tolerance:  # convergence satisfied
                    self.n_reset = 0
                    if LOGGING_ON:
                        log.info( "time step equilibrated in %d step(s)", self.k )
                    break                        # update_switch -> on

                self.solv_timer.reset()
                #self.d_U = solve(self.K,self.R_k)
                #s_solve = linsolve.factorized(self.K)
                #self.d_U = s_solve(self.R_k)
                self.d_U = linsolve.spsolve(self.K,self.R_k) # DG_k * d_U = r
                self.solv_timer.record()

                self.U_k += self.d_U  

                self.k      += 1
                self.tot_k  += 1
                self.ttot_k += 1
                step_flag = 'corrector'

            else:                                # handling nonconverged step

                if LOGGING_ON:
                    log.info("no convergence reached - refinement in time")
                self.n_reset = self.n_reset + 1  # adaptive strategy halving d_t
                if self.n_reset > self.RESETMAX: # max number of resets achieved

                    # Handle this situation with an exception with an
                    # associated dialog box and concise report of the
                    # iteration process
                    #
                    # @TODO raise the NoConvergence exception
                    if not self.sync_resp_tracing:
                        self.rtrace_mngr.stop_timer()
                        self.rtrace_mngr.timer_tick()
                    return
                
                self.d_t *= adap.fhandler_get_scale()
                self.t_n1 = self.t_n + self.d_t

                if LOGGING_ON:
                    log.info("redo time step with dt = %s", )
                continue

            self.iter_timer.record()
            #self.mem_counter.record()

            if self.sync_resp_tracing:
                self.rtrace_mngr.timer_tick()

            if adap.ehandler_needed(): # explicit adaptations 
                if adap.ehandler_accept():  # accept equilibrium?
                    self.accept_time_step() # register the state and response
                adap.ehandler_invoke()
            else:
                self.accept_time_step()
                self.d_t *= adap.ehandler_get_scale()
                self.t_n  = self.t_n1
                self.t_n1 = self.t_n + self.d_t
                self.T.val = self.t_n1
                
        # Report computation time
        #
        self.eval_timer.record()
        self.eval_timer.report()

        if not self.sync_resp_tracing:
            self.rtrace_mngr.stop_timer()
            self.rtrace_mngr.timer_tick()

        return

    def accept_time_step(self):
        self.ts.update_state( self.U_k )
        self.U_n = copy( self.U_k )
        if LOGGING_ON:
            log.debug("u_n = \n%s", str(self.U_n))

        self.rtrace_mngr_timer.reset()
        self.rtrace_mngr.record_equilibrium( self.U_k )
        self.rtrace_mngr_timer.record()
        
    def report_load_step(self):
        print '--- LS: %3d ------------------------' % self.ls_counter
        
    def report_iteration(self):
        print '--- LS: %3d ---- IT: %3d -----------' % (self.ls_counter, self.k)
 


    calculate = Button
    computation_thread = Instance(Thread)

    @on_trait_change('calculate')
    def teval( self ):
        '''
        Eval in a thread.
        '''
        if self.computation_thread and self.computation_thread.isAlive():
            self.user_wants_abort = True
        else:
            self.user_wants_abort = False
            self.computation_thread = Thread(target=self.eval, name="computation")
            self.computation_thread.start()

    view = View( Group( Item('calculate', show_label = False ),
                        'sync_resp_tracing',
                        Item('DT',label = 'Time step size'),
                        Item('T', style='custom'),
                        HGroup( Item('KMAX',label = 'Max. number of iterations'),
                                Item('RESETMAX',label = 'Max. number of resets' ),
                                Item('tolerance',label = 'Tolerance on residual norm')
                                ),
                        ),
                 #Item('rmgr', style="custom"),
                 #handler = TCHandler(),
                 resizable=True,
                 height=0.75, width=0.75,
                 handler = TLoopHandler(),
                 buttons= [OKButton,CancelButton,RecalcAction])

if LOGGING_ON:
    import logging.config
    logging.config.fileConfig("logging.conf")
