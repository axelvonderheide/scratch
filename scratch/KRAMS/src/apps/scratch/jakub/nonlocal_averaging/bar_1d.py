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

def run():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats1D.mats1D_elastic.mats1D_elastic import MATS1DElastic
    from ibvpy.mats.mats1D.mats1D_damage.mats1D_damage import MATS1DDamage
    from ibvpy.fets.fets1D.fets1D2l import FETS1D2L

    #fets_eval = FETS1D2L( mats_eval = MATS1DElastic( E = 10. ) )
    fets_eval = FETS1D2L( mats_eval = MATS1DDamage( E = 10. ) )
    from ibvpy.mesh.fe_grid import FEGrid

    # Discretization
    domain = FEGrid( coord_max = ( 1., 0., 0. ),
                           shape = ( 10, ),
                           fets_eval = fets_eval )

    right_dof = domain[-1, -1].dofs[0, 0, 0]

    ts = TS( nonlocal_avg = True,
             dof_resultants = True,
             sdomain = domain,
         # conversion to list (square brackets) is only necessary for slicing of 
         # single dofs, e.g "get_left_dofs()[0,1]"
#         bcond_list =  [ BCDof(var='u', dof = 0, value = 0.)     ] +  
#                    [ BCDof(var='u', dof = 2, value = 0.001 ) ]+
#                    [ )     ],
         bcond_list = [BCDof( var = 'u', dof = 0, value = 0. ),
#                        BCDof(var='u', dof = 1, link_dofs = [2], link_coeffs = [0.5],
#                              value = 0. ),
#                        BCDof(var='u', dof = 2, link_dofs = [3], link_coeffs = [1.],
#                              value = 0. ),
                        BCDof( var = 'f', dof = right_dof, value = 1,
                                  #link_dofs = [2], link_coeffs = [2]
                                   ) ],
         rtrace_list = [ RTraceGraph( name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = 0,
                               var_x = 'U_k', idx_x = 1 ),
                    RTraceDomainListField( name = 'Stress' ,
                         var = 'sig_app', idx = 0 ),
                     RTraceDomainListField( name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    warp = True ),
                             RTraceDomainListField( name = 'N0' ,
                                          var = 'N_mtx', idx = 0,
                                          record_on = 'update' )

                ]
            )

    # Add the time-loop control
    tloop = TLoop( tstepper = ts,
                   tline = TLine( min = 0.0, step = 0.5, max = 1.0 ) )

    print '---- result ----'
    print tloop.eval()
    print ts.F_int
    print ts.rtrace_list[0].trace.ydata

if __name__ == '__main__':
    run()
