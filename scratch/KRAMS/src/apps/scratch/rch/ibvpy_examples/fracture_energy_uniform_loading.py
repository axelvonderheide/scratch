#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 10, 2009 by: rch

from enthought.traits.api import \
    HasTraits, Property, Float, Trait, cached_property, Instance, Int

from enthought.traits.ui.api import \
    View, Item, Group, HGroup, VGroup
     
from ibvpy.util.simgrid import \
    simgrid

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.mats.mats3D.mats3D_cmdm.mats3D_cmdm import \
    MATS3DMicroplaneDamage, PhiFnStrainSoftening

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H

from ibvpy.fets.fets2D5.fets2D58h import \
    FETS2D58H

from ibvpy.fets.fets2D.fets2Drotsym import \
    FETS2Drotsym

from ibvpy.fets.fets2D.fets2D4q import \
    FETS2D4Q

from ibvpy.api import \
    RTraceGraph
    
from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray
    
from mathkit.mfn.mfn_line.mfn_matplotlib_editor import \
    MFnMatplotlibEditor
    
from mathkit.mfn.mfn_line.mfn_chaco_editor import \
    MFnChacoEditor
    
from numpy import \
    column_stack, sum
    
from math import \
    pi as Pi, sqrt

class CMDMFractureEnergyTester( HasTraits ):

    # dimension-independent attributes
    
    mats = Instance( MATS3DMicroplaneDamage )
    def _mats_default(self):
        return MATS3DMicroplaneDamage( E = 34000, nu = 0.25,
                                       model_version = 'stiffness',
                                       phi_fn = PhiFnStrainSoftening(
                                                                     G_f = 0.001117,
                                                                     f_t = 2.8968,
                                                                     md = 0.0 ) )

    u_max   = Float( 0.0005, enter_set = True, auto_set = False,
                     desc = 'maximum control displacement' )
    
    n_steps =   Int( 10, enter_set = True, auto_set = False,
                     desc = 'number of load steps to reach maximum displacement' )

    _results_3d = Property( depends_on = 'n_steps, u_max, mats.+' )
    @cached_property
    def _get__results_3d( self ):
        
        fets_eval = FETS3D8H( mats_eval  = self.mats )        

        support_slice =  [ 
                (0   ,slice(None),slice(None),0   ,slice(None),slice(None)), # yz plane  0
                (0   ,0   ,slice(None),0   ,0   ,slice(None)), #  z-axis   1
                (0   ,0   ,   0,0   ,0   ,0   )  #  origin   2
                          ]
        
        support_dirs = [[0],[1],[2]]
        
        loading_slice = (-1  ,slice(None),slice(None),-1  ,slice(None),slice(None))  # loading in x dir

        tl, u1, fields, integs, xydata = simgrid( fets_eval, (1,1,1), (1,1,1), 
                                  support_slice, support_dirs,
                                  loading_slice, 
                                  0, 
                                  self.u_max, 
                                  self.n_steps,
                                  vars  = ['fracture_energy'],
                                  ivars = ['fracture_energy'] )

        return integs[0][0], xydata

    _results_2d = Property( depends_on = 'n_steps, u_max, mats.+' )
    @cached_property
    def _get__results_2d( self ):

        fets_eval = FETS2Drotsym( prototype_fets = FETS2D4Q(),
                                  mats_eval  = self.mats )        

        support_slice = [
                          (0   ,slice(None),0   ,slice(None)), #  y-axis   1
                          (0   ,     0,    0   ,0   )  #  origin   2
                          ]
        support_dirs = [[0],[1]]
    
        loading_slice = (-1  ,slice(None),-1  ,slice(None))  # loading in x dir

        # Get the radius for a circle with the surface equal to 1
        # A = Pi * r^2 --> r = sqrt( A / Pi )
        R = sqrt( 1. / Pi )
        tl, u1, fields, integs, xydata = simgrid( fets_eval, (1,R), (1,1), 
                              support_slice, support_dirs,
                              loading_slice, 
                              0, 
                              self.u_max, 
                              self.n_steps,
                              vars  = ['fracture_energy'],
                              ivars = ['fracture_energy'] )

        return integs[0][0], xydata        

    fracture_energy3d = Property(
                                 depends_on = 'n_steps, u_max, mats.+' )
    @cached_property
    def _get_fracture_energy3d(self):
        return self._results_3d[0]

    graph3d = Property( Instance( MFnLineArray ),
                        depends_on = 'n_steps, u_max, mats.+' )
    @cached_property
    def _get_graph3d(self):
        graph = self._results_3d[1]
        return MFNLineArray( xdata = graph[0], ydata = graph[1] ) 
    
    fracture_energy2d = Property(
                                 depends_on = 'n_steps, u_max, mats.+' )
    @cached_property
    def _get_fracture_energy2d(self):
        return self._results_2d[0]

    graph2d = Property( Instance( MFnLineArray ),
                        depends_on = 'n_steps, u_max, mats.+' )
    @cached_property
    def _get_graph2d(self):
        graph = self._results_2d[1]
        return MFNLineArray( xdata = graph[0], ydata = graph[1] ) 

    traits_view = View( Item('mats', show_label = False ),
                        HGroup( Item('u_max', label = 'maximum displacement'),
                                Item('n_steps', label = 'number of load steps' ),
                        ),
                        HGroup( 
                        VGroup( Item('fracture_energy3d', label = 'Gf_3D'),
                               Item('graph3d', editor = MFnMatplotlibEditor(), 
                                    show_label = False ),
                               label = 'Three-dimensional'                                
                               ),
                        VGroup( Item('fracture_energy2d', label = 'Gf_2D'),
                               Item('graph2d', editor = MFnMatplotlibEditor(), 
                                    show_label = False ),
                               label = 'Rotational symmetry'
                               ),                                
                        ),
                        resizable = True,
                        height = 0.9, width = 0.9,
                        buttons = ['OK','Cancel']
                        )

    mfn = Property( Instance( MFnLineArray ),
                    depends_on = 'n_steps, u_max' )
    @cached_property
    def _get_mfn( self ):
        graph2d = self._results_2d[1]
        graph3d = self._results_3d[1]
        
        ydata = column_stack( [sum( graph2d[1], axis = 1),
                               sum( graph3d[1], axis = 1) ] )
        
        mfn = MFnLineArray( xdata = graph2d[0], ydata = ydata )
        return mfn

    straits_view = View( Item('mats', show_label = False ),
                        HGroup( Item('u_max', label = 'maximum displacement'),
                                Item('n_steps', label = 'number of load steps' ),
                        ),
                        HGroup( Item('fracture_energy3d', label = 'Gf_3D'),
                                Item('fracture_energy2d', label = 'Gf_2D'),
                        ),
                        Item('mfn', editor = MFnMatplotlibEditor(), 
                                    show_label = False ),
                        resizable = True,
                        height = 0.9, width = 0.9,
                        buttons = ['OK','Cancel']
                        )

if __name__ == '__main__':
    cmdm = CMDMFractureEnergyTester()
    cmdm.configure_traits( view = 'straits_view' )    
