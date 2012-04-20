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
# Created on Jul 29, 2009 by: rchx

#----------------------- example --------------------

if __name__ == '__main__':

    from math import fabs
    from numpy import loadtxt
    from os.path import join
        
#   from ibvpy.fets.fets2D5.fets2D58h import FETS2D58H
    #from ibvpy.bcond.bc_slice import BCSlice    
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDofGroup, IBVPSolve as IS, DOTSEval, BCSlice   
    from ibvpy.fets.fets3D.fets3D8h import FETS3D8H
        
    from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import MATS3DElastic
          
    dim = '3D'
    if dim == '3D':
        fets_eval = FETS3D8H(mats_eval = MATS3DElastic( E = 35000, nu = 0.2)) # [N/mm], [-]
    fets_eval.vtk_r *= 1.0

                                            
    from ibvpy.mesh.fe_grid import FEGrid
    
   
    # Discretization
    domain = FEGrid( coord_max = (1.,1.,0.5 ), 
                     shape   = (1, 1, 1),
                     fets_eval = fets_eval)
                                 
    ts = TS(
            sdomain = domain,
            bcond_list = [BCSlice( var = 'u', value = 0., dims = [0,1,2], slice = domain[0,:,:,0,:,:] ),
                          BCSlice( var = 'f', value = -1.0, dims = [2], slice = domain[-1,:,:,-1,:,-1] )],
             rtrace_list = [ 
                            
                         RTraceDomainListField(name = 'Displacement' ,
                                        var = 'u', idx = 0, warp = True),
                         RTraceDomainListField(name = 'Stress' ,
                                        var = 'sig_app', idx = 0, warp = True, 
                                        record_on = 'update'),
                        ]             
            )
    
    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
                   #tolerance = 1e-3,
                   tline  = TLine( min = 0.0,  step = 1., max = 1.0 ))


    tloop.eval()

  

    
    # Put the whole stuff into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp( ibv_resource = tloop )
    app.main()
