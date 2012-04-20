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
# Created on Jun 15, 2010 by: rch

# Construct a single domain

from ibvpy.fets.fets2D.fets2D4q import \
    FETS2D4Q

from ibvpy.api import \
    TStepper as TS, RTraceGraph, RTraceDomainListField, \
    RTraceDomainListInteg, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import \
    MATS2DElastic

from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import \
    MATS2DScalarDamage

from numpy import \
    array

from ibvpy.api import \
    BCDofGroup

from ibvpy.mesh.fe_grid import \
    FEGrid

from ibvpy.mesh.fe_refinement_grid import \
    FERefinementGrid

from ibvpy.mesh.fe_domain import \
    FEDomain

from mathkit.mfn.mfn_line.mfn_line import \
    MFnLineArray

def run():
    fets_eval = FETS2D4Q( mats_eval = MATS2DElastic() )

    # Discretization
    fe_grid = FEGrid( coord_max = ( 10., 4., 0. ),
                      shape = ( 10, 3 ),
                      fets_eval = fets_eval )

    bcg = BCDofGroup( var = 'u', value = 0., dims = [0],
                   get_dof_method = fe_grid.get_left_dofs )
    bcg.setup( None )
    print 'labels', bcg._get_labels()
    print 'points', bcg._get_mvpoints()

    mf = MFnLineArray( #xdata = arange(10),
                       ydata = array( [0, 1, 2, 3] ) )

    right_dof = 2
    tstepper = TS( sdomain = fe_grid,
                   bcond_list = [ BCDofGroup( var = 'u', value = 0., dims = [0, 1],
                                               get_dof_method = fe_grid.get_left_dofs ),
                                 BCDofGroup( var = 'u', value = .005, dims = [1],
                                          time_function = mf.get_value,
                                          get_dof_method = fe_grid.get_right_dofs ) ],
            )

    # Add the time-loop control
    tloop = TLoop( tstepper = tstepper, KMAX = 300, tolerance = 1e-4,
                   tline = TLine( min = 0.0, step = 1.0, max = 1.0 ) )

    U_k = tloop.eval()
    print 'dir', tloop.rtrace_mngr.dir

    # RTrace should not contain backward link to RTraceMngr
    # The definition should be forward. 
    #

    rt1 = RTraceGraph( sd = tstepper.sdomain,
                       rmgr = tstepper.rtrace_mngr,
                       name = 'Fi,right over u_right (iteration)' ,
                       var_y = 'F_int', idx_y = right_dof,
                       var_x = 'U_k', idx_x = right_dof,
                       record_on = 'update' )
    print 'dir', rt1.dir

if __name__ == '__main__':
    run()
